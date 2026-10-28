/**
 * _fused_eval.c -- the fused evaluate of the S2 key FusedCost (docs/todo_s2_fused_c_eval.md):
 * one C call for VoxelCostFunction.evaluate in hard mode without a ray rule. Compiled into the
 * _rasterize module (setup.py: a second source, with ICENINE_FUSED defined, which adds the method
 * table entries and the HAS_FUSED constant in _rasterize.c).
 *
 *   fused_context(...) -> capsule   the per-(cost function, phase) context (cost_functions.py,
 *                                   _FusedContext, builds it and rebuilds it when it is stale)
 *   fused_evaluate(ctx, R, verts[, debug]) -> the 7 OverlapInfo values, or the list of frames the
 *                                   FrameTable must fill first (Python calls ensure and repeats the
 *                                   call); debug=1: a dict of the intermediates
 *
 * What it computes (cost_functions.py VoxelCostFunction.evaluate and
 * calculate_diffraction_overlap_batched; diffraction_core.py scattering_omegas_fixed_beam): g = R
 * g_hkl and the Bragg solve; the observable rows (omega1 of every observable reflection, then
 * omega2, as torch.cat); the lab normal, the reflected ray, eta and the min_sin_eta filter; the
 * frame lookup (frame_of_omegas); stages B and C (lab vertices, ray-plane intersection, pixel
 * coordinate with the CppEdgeRules -1 rule); stage D through icenine_stage_d_peak, the per-peak
 * body that stage_d_overlap uses too (_rasterize.c), so stage D exists once.
 *
 * Float32 throughout, in the operation order of the torch code. Where torch's kernels fuse a
 * multiply and an add (bmm, the Rz @ base product, the einsum of the vertices, torch.norm) the code
 * calls fmaf explicitly in the order that matches torch bit for bit (S2 Step 0 probe,
 * benchmarks/s2_step0/fmaf_probe.json); FP contraction is off in this file, so the compiler fuses
 * nothing else. Not bit for bit: the 3 x K product R @ g_hkl^T (torch: BLAS sgemm; here the plain
 * order by default, glab_order 1, which gave fewer flips than the prototype's fmaf chain, S2 Step
 * 1) and sin, cos, atan2, asin (torch: vectorised kernels; here: libm). The differences are float32
 * rounding, so a pixel coordinate within a few ulp of an integer can floor to the other pixel (the
 * parity rule, D-20261027-I2: such pixel-edge flips are accepted; frame, eta, observability and
 * t > 0 flips are not).
 *
 * Validated on macOS clang (arm64) only; GCC and other platforms are untested (the float order is
 * pinned by the explicit fmaf calls and -ffp-contract=off, but libm's sin, cos, atan2 and asin
 * differ between platforms).
 *
 * Not thread-safe: the scratch buffers live in the context. The GIL is held for the whole call,
 * and in the non-debug path no Python code runs between filling the scratch buffers and reading
 * them (rows and cols are attribute reads; Python objects are created only after stage D), so
 * calls on one context are serialised. A cost function shared between threads of a free-threaded
 * Python would need per-call scratch.
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NO_IMPORT_ARRAY /* imported by the module init in _rasterize.c */
#include "_stage_d.h" /* numpy/arrayobject.h, MAX_DETECTORS, IceStageDAcc, icenine_stage_d_peak */
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* no multiply-add contraction: every fused multiply-add is an explicit fmaf. setup.py also passes
   -ffp-contract=off for the whole extension (D-20261027-I8), since GCC ignores STDC FP_CONTRACT
   and documents its optimize pragma as unreliable; the pragmas stay for builds without the flag */
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif
#pragma STDC FP_CONTRACT OFF

#define FUSED_MAX_DET MAX_DETECTORS /* _stage_d.h */
#define FUSED_CAPSULE "icenine._rasterize.fused_ctx"

typedef struct {
    /* phase: K reflections */
    int K;
    const float *g_hkl;     /* (K, 3) */
    const float *sin_theta; /* (K,) */
    const float *g_mag_eps; /* (K,) */
    float sin_chi_laue, cos_chi_laue, eps;
    /* beam and sample */
    float beam[3];
    float base[9]; /* sample_to_lab[:3, :3] row-major */
    float trans[3];
    /* filters */
    float eta_limit;
    float min_sin_eta;
    int use_min_sin_eta;
    /* modes */
    int pixel_radius;
    int center_screen;
    int rule_bits; /* 1 C++ CountQualifiedPeaks rule, 2 CppEdgeRules */
    int glab_order; /* g = R g_hkl: 0 the fmaf chain, 1 the plain order (a*b + c*d) + e*f */
    /* frame lookup: f = (omega - low) / width in the dtypes numpy uses (cost_functions.py,
       _numpy_frame_mode) */
    double low, width;
    float low_f, width_f;
    int frame_mode; /* 0: sub and div in float64; 1: sub float32, div float64; 2: all float32 */
    const long long *frame_index;
    long long n_index;
    /* detectors */
    int n_det;
    float det_n[FUSED_MAX_DET][3], det_d[FUSED_MAX_DET];
    float det_pos[FUSED_MAX_DET][3], det_org[FUSED_MAX_DET][3];
    float det_bj[FUSED_MAX_DET][3], det_bk[FUSED_MAX_DET][3];
    float det_hw[FUSED_MAX_DET], det_pw[FUSED_MAX_DET], det_hh[FUSED_MAX_DET],
        det_ph[FUSED_MAX_DET];
    /* frame table (FrameTable): the images list, the filled flags, the table object (rows, cols) */
    PyObject *table;
    PyObject *keep[5]; /* g_hkl, sin_theta, g_mag_eps, frame_index, filled */
    PyObject *images;
    const npy_bool *filled;
    long long n_frames;
    /* scratch, sized for N = 2K rows */
    float *g_lab; /* (K, 3) */
    unsigned char *obs;
    float *om1, *om2;                                 /* (K,) */
    float *row_om, *row_n, *row_eta, *row_c, *row_s; /* (N,), (N, 3), (N,), (N,), (N,) */
    float *row_refl;                                  /* (N, 3) */
    unsigned char *row_valid;
    int *row_frame;     /* (N,) -1: none */
    float *pix;         /* (N, n_det, 3, 2) [col, row] */
    unsigned char *hit; /* (N, n_det) all three vertices hit */
    float *tval;        /* (N, n_det, 3) */
} FusedCtx;

static void fused_free(PyObject *cap) {
    FusedCtx *c = (FusedCtx *)PyCapsule_GetPointer(cap, FUSED_CAPSULE);
    if (!c) return;
    Py_XDECREF(c->table);
    for (int i = 0; i < 5; i++) Py_XDECREF(c->keep[i]);
    Py_XDECREF(c->images);
    free(c->g_lab);
    free(c->obs);
    free(c->om1);
    free(c->om2);
    free(c->row_om);
    free(c->row_n);
    free(c->row_eta);
    free(c->row_c);
    free(c->row_s);
    free(c->row_refl);
    free(c->row_valid);
    free(c->row_frame);
    free(c->pix);
    free(c->hit);
    free(c->tval);
    free(c);
}

static int fused_f32(PyObject *o, int n, float *out, const char *name) {
    PyArrayObject *a = (PyArrayObject *)PyArray_FROM_OTF(o, NPY_FLOAT32, NPY_ARRAY_IN_ARRAY);
    if (!a) return -1;
    if (PyArray_SIZE(a) != n) {
        PyErr_Format(PyExc_ValueError, "%s: expected %d float32 values", name, n);
        Py_DECREF(a);
        return -1;
    }
    memcpy(out, PyArray_DATA(a), n * sizeof(float));
    Py_DECREF(a);
    return 0;
}

/**
 * fused_context(g_hkl, sin_theta, g_mag_eps, scal, beam, base, trans, ints, frame_f, frame_index,
 *               det_geom, table, images, filled) -> capsule
 * g_hkl (K, 3) float32, sin_theta and g_mag_eps (K,) float32. scal float32 [sin_chi_laue,
 * cos_chi_laue, eps, eta_limit, min_sin_eta]; ints [pixel_radius, center_screen, rule_bits,
 * use_min_sin_eta, frame_mode, n_det, glab_order]; frame_f float64 [low, width]; frame_index int64;
 * det_geom (n_det, 20) float32 n[3] d pos[3] org[3] bj[3] bk[3] hw pw hh ph; table: the
 * FrameTable (rows and cols read on each call), images: its list, filled: its bool array. The
 * context keeps a reference to every array it points into.
 */
PyObject *icenine_fused_context(PyObject *self, PyObject *args) {
    PyArrayObject *g_hkl, *sin_theta, *g_mag_eps, *frame_index, *filled;
    PyObject *scal, *beam, *base, *trans, *ints_o, *frame_f, *geom, *table, *images;
    (void)self;
    if (!PyArg_ParseTuple(args, "O!O!O!OOOOOOO!OOO!O!", &PyArray_Type, &g_hkl, &PyArray_Type,
                          &sin_theta, &PyArray_Type, &g_mag_eps, &scal, &beam, &base, &trans,
                          &ints_o, &frame_f, &PyArray_Type, &frame_index, &geom, &table,
                          &PyList_Type, &images, &PyArray_Type, &filled))
        return NULL;
    if (PyArray_NDIM(g_hkl) != 2) {
        PyErr_SetString(PyExc_TypeError, "fused_context: g_hkl must be (K, 3)");
        return NULL;
    }
    int K = (int)PyArray_DIM(g_hkl, 0);
    if (PyArray_TYPE(g_hkl) != NPY_FLOAT32 || !PyArray_IS_C_CONTIGUOUS(g_hkl) ||
        PyArray_DIM(g_hkl, 1) != 3 || PyArray_TYPE(sin_theta) != NPY_FLOAT32 ||
        !PyArray_IS_C_CONTIGUOUS(sin_theta) || PyArray_SIZE(sin_theta) != K ||
        PyArray_TYPE(g_mag_eps) != NPY_FLOAT32 || !PyArray_IS_C_CONTIGUOUS(g_mag_eps) ||
        PyArray_SIZE(g_mag_eps) != K || PyArray_TYPE(frame_index) != NPY_INT64 ||
        !PyArray_IS_C_CONTIGUOUS(frame_index) || PyArray_SIZE(frame_index) < 1 ||
        PyArray_TYPE(filled) != NPY_BOOL || !PyArray_IS_C_CONTIGUOUS(filled)) {
        PyErr_SetString(PyExc_TypeError, "fused_context: bad array type or shape");
        return NULL;
    }
    int iv[7];
    double ff[2];
    {
        PyArrayObject *a = (PyArrayObject *)PyArray_FROM_OTF(ints_o, NPY_INT32, NPY_ARRAY_IN_ARRAY);
        PyArrayObject *b =
            (PyArrayObject *)PyArray_FROM_OTF(frame_f, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
        if (!a || !b || PyArray_SIZE(a) != 7 || PyArray_SIZE(b) != 2) {
            Py_XDECREF(a);
            Py_XDECREF(b);
            if (!PyErr_Occurred())
                PyErr_SetString(PyExc_ValueError, "fused_context: ints (7) / frame_f (2)");
            return NULL;
        }
        memcpy(iv, PyArray_DATA(a), sizeof(iv));
        memcpy(ff, PyArray_DATA(b), sizeof(ff));
        Py_DECREF(a);
        Py_DECREF(b);
    }
    int n_det = iv[5];
    if (n_det < 1 || n_det > FUSED_MAX_DET) {
        PyErr_SetString(PyExc_ValueError, "fused_context: n_det out of range (1-8)");
        return NULL;
    }
    float sc[5];
    float g[FUSED_MAX_DET * 20];
    FusedCtx *c = (FusedCtx *)calloc(1, sizeof(FusedCtx));
    if (!c) return PyErr_NoMemory();
    if (fused_f32(scal, 5, sc, "scal") || fused_f32(beam, 3, c->beam, "beam") ||
        fused_f32(base, 9, c->base, "base") || fused_f32(trans, 3, c->trans, "trans") ||
        fused_f32(geom, n_det * 20, g, "det_geom")) {
        free(c);
        return NULL;
    }
    c->K = K;
    c->g_hkl = (const float *)PyArray_DATA(g_hkl);
    c->sin_theta = (const float *)PyArray_DATA(sin_theta);
    c->g_mag_eps = (const float *)PyArray_DATA(g_mag_eps);
    c->sin_chi_laue = sc[0];
    c->cos_chi_laue = sc[1];
    c->eps = sc[2];
    c->eta_limit = sc[3];
    c->min_sin_eta = sc[4];
    c->pixel_radius = iv[0];
    c->center_screen = iv[1];
    c->rule_bits = iv[2];
    c->use_min_sin_eta = iv[3];
    c->frame_mode = iv[4];
    c->n_det = n_det;
    c->glab_order = iv[6];
    c->low = ff[0];
    c->width = ff[1];
    c->low_f = (float)ff[0];
    c->width_f = (float)ff[1];
    c->frame_index = (const long long *)PyArray_DATA(frame_index);
    c->n_index = PyArray_SIZE(frame_index);
    for (int d = 0; d < n_det; d++) {
        const float *q = g + 20 * d;
        for (int i = 0; i < 3; i++) {
            c->det_n[d][i] = q[i];
            c->det_pos[d][i] = q[4 + i];
            c->det_org[d][i] = q[7 + i];
            c->det_bj[d][i] = q[10 + i];
            c->det_bk[d][i] = q[13 + i];
        }
        c->det_d[d] = q[3];
        c->det_hw[d] = q[16];
        c->det_pw[d] = q[17];
        c->det_hh[d] = q[18];
        c->det_ph[d] = q[19];
    }
    Py_INCREF(table);
    c->table = table;
    Py_INCREF(images);
    c->images = images;
    c->filled = (const npy_bool *)PyArray_DATA(filled);
    c->n_frames = PyArray_SIZE(filled);
    c->keep[0] = (PyObject *)g_hkl;
    c->keep[1] = (PyObject *)sin_theta;
    c->keep[2] = (PyObject *)g_mag_eps;
    c->keep[3] = (PyObject *)frame_index;
    c->keep[4] = (PyObject *)filled;
    for (int i = 0; i < 5; i++) Py_INCREF(c->keep[i]);
    /* from here on fused_free releases everything, through the capsule */
    PyObject *cap = PyCapsule_New(c, FUSED_CAPSULE, fused_free);
    if (!cap) {
        Py_DECREF(table);
        Py_DECREF(images);
        for (int i = 0; i < 5; i++) Py_DECREF(c->keep[i]);
        free(c);
        return NULL;
    }
    size_t N = 2 * (size_t)K + 1;
    c->g_lab = (float *)malloc(3 * (K + 1) * sizeof(float));
    c->obs = (unsigned char *)malloc(K + 1);
    c->om1 = (float *)malloc((K + 1) * sizeof(float));
    c->om2 = (float *)malloc((K + 1) * sizeof(float));
    c->row_om = (float *)malloc(N * sizeof(float));
    c->row_n = (float *)malloc(3 * N * sizeof(float));
    c->row_eta = (float *)malloc(N * sizeof(float));
    c->row_c = (float *)malloc(N * sizeof(float));
    c->row_s = (float *)malloc(N * sizeof(float));
    c->row_refl = (float *)malloc(3 * N * sizeof(float));
    c->row_valid = (unsigned char *)malloc(N);
    c->row_frame = (int *)malloc(N * sizeof(int));
    c->pix = (float *)malloc(N * n_det * 6 * sizeof(float));
    c->hit = (unsigned char *)malloc(N * n_det);
    c->tval = (float *)malloc(N * n_det * 3 * sizeof(float));
    if (!c->g_lab || !c->obs || !c->om1 || !c->om2 || !c->row_om || !c->row_n || !c->row_eta ||
        !c->row_c || !c->row_s || !c->row_refl || !c->row_valid || !c->row_frame || !c->pix ||
        !c->hit || !c->tval) {
        Py_DECREF(cap);
        return PyErr_NoMemory();
    }
    return cap;
}

/* frame of an omega: _omega_bins + frame_of_omegas (cost_functions.py) */
static int fused_frame(const FusedCtx *c, float om) {
    double f;
    if (c->frame_mode == 0) {
        f = ((double)om - c->low) / c->width;
    } else if (c->frame_mode == 1) {
        f = (double)(om - c->low_f) / c->width;
    } else {
        f = (double)((om - c->low_f) / c->width_f);
    }
    long long b;
    if (c->rule_bits & 2) f = floor(f); /* CppEdgeRules: floor; default: truncation (astype) */
    /* astype(int) truncates toward zero; the rows that reach here have a finite omega in
       [-pi, pi] (a NaN omega fails the eta filter first), the guards are for safety only */
    if (!(f > -1e15))
        b = -1000000000000000LL;
    else if (f > 1e15)
        b = 1000000000000000LL;
    else
        b = (long long)f;
    if (b < -1) b = -1;
    if (b > c->n_index - 1) b = c->n_index - 1;
    long long idx = b < 0 ? c->n_index + b : b; /* index -1 is the last (-1) entry */
    return (int)c->frame_index[idx];
}

static inline float fused_sum3(float a, float b, float d) { return (a + b) + d; }

/**
 * fused_evaluate(ctx, R, verts, debug=0) -> (po, pon, pko, pkon, quality, nq, det_ovlp)
 *   or a list of frame indices the FrameTable must fill first (call again after ensure);
 *   with debug=1 a dict of the intermediates instead ("result" is None and "missing" lists the
 *   frames while some are not filled).
 * R: (3, 3) array (any float dtype; cast to float64, then to float32 as torch .float()).
 * verts: (3, 3) float32 C-contiguous.
 */
PyObject *icenine_fused_evaluate(PyObject *self, PyObject *args) {
    PyObject *cap, *R_o;
    PyArrayObject *V_a;
    int debug = 0;
    (void)self;
    if (!PyArg_ParseTuple(args, "OOO!|i", &cap, &R_o, &PyArray_Type, &V_a, &debug)) return NULL;
    FusedCtx *c = (FusedCtx *)PyCapsule_GetPointer(cap, FUSED_CAPSULE);
    if (!c) return NULL;
    if (PyArray_TYPE(V_a) != NPY_FLOAT32 || !PyArray_IS_C_CONTIGUOUS(V_a) ||
        PyArray_SIZE(V_a) != 9) {
        PyErr_SetString(PyExc_TypeError, "verts must be (3, 3) C-contiguous float32");
        return NULL;
    }
    float Rf[9];
    if (PyArray_Check(R_o) && PyArray_TYPE((PyArrayObject *)R_o) == NPY_FLOAT64 &&
        PyArray_IS_C_CONTIGUOUS((PyArrayObject *)R_o) && PyArray_SIZE((PyArrayObject *)R_o) == 9) {
        const double *Rd = (const double *)PyArray_DATA((PyArrayObject *)R_o);
        for (int i = 0; i < 9; i++) Rf[i] = (float)Rd[i];
    } else {
        PyArrayObject *a = (PyArrayObject *)PyArray_FROM_OTF(R_o, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
        if (!a) return NULL;
        if (PyArray_SIZE(a) != 9) {
            Py_DECREF(a);
            PyErr_SetString(PyExc_ValueError, "R must have 9 elements");
            return NULL;
        }
        const double *Rd = (const double *)PyArray_DATA(a);
        for (int i = 0; i < 9; i++) Rf[i] = (float)Rd[i];
        Py_DECREF(a);
    }
    const float *V0 = (const float *)PyArray_DATA(V_a);
    float V[9];
    if (c->center_screen && c->pixel_radius > 0) {
        /* verts.mean(dim=0, keepdim=True).repeat(3, 1): ((a + b) + c) / 3 */
        for (int i = 0; i < 3; i++) {
            float m = fused_sum3(V0[i], V0[3 + i], V0[6 + i]) / 3.0f;
            V[i] = V[3 + i] = V[6 + i] = m;
        }
    } else {
        memcpy(V, V0, sizeof(V));
    }

    const int K = c->K, n_det = c->n_det;
    const float eps = c->eps;
    const float PI_F = (float)M_PI, TWO_PI_F = (float)(2.0 * M_PI);
    /* 1. g_lab = (R @ g^T)^T, and the Bragg solve (scattering_omegas_fixed_beam) */
    for (int k = 0; k < K; k++) {
        float *gl = c->g_lab + 3 * k;
        const float *g = c->g_hkl + 3 * k;
        if (c->glab_order == 1) {
            for (int i = 0; i < 3; i++)
                gl[i] = (Rf[3 * i] * g[0] + Rf[3 * i + 1] * g[1]) + Rf[3 * i + 2] * g[2];
        } else {
            for (int i = 0; i < 3; i++)
                gl[i] = fmaf(Rf[3 * i + 2], g[2], fmaf(Rf[3 * i + 1], g[1], Rf[3 * i] * g[0]));
        }
        float cos_chi = gl[2] / c->g_mag_eps[k];
        cos_chi = cos_chi < -1.0f ? -1.0f : (cos_chi > 1.0f ? 1.0f : cos_chi);
        float sin_chi = sqrtf((1.0f - cos_chi * cos_chi) + eps);
        float num = c->sin_theta[k] + cos_chi * c->sin_chi_laue;
        float den = sin_chi * c->cos_chi_laue + eps;
        unsigned char ob = fabsf(num) <= fabsf(den);
        c->obs[k] = ob;
        if (!ob && !debug) continue;
        float d0 = atan2f(gl[0], gl[1]);
        float arg = num / den;
        arg = arg < -1.0f ? -1.0f : (arg > 1.0f ? 1.0f : arg);
        float b1 = asinf(arg);
        float b2 = PI_F - b1;
        float o1 = b1 + d0, o2 = b2 + d0;
        if (o1 > PI_F) o1 = o1 - TWO_PI_F;
        if (o1 < -PI_F) o1 = o1 + TWO_PI_F;
        if (o2 > PI_F) o2 = o2 - TWO_PI_F;
        if (o2 < -PI_F) o2 = o2 + TWO_PI_F;
        c->om1[k] = o1;
        c->om2[k] = o2;
    }

    /* 2. rows: omega1 of the observable reflections, then omega2 (torch.cat order); normals,
       Rz @ base, lab normal, reflection, eta filter, frame. In debug mode every (solution,
       reflection) row is computed, row r = s * K + k, and rows of unobservable reflections are
       marked invalid */
    const float *B = c->base;
    int n_rows = 0;
    for (int s = 0; s < 2; s++) {
        for (int k = 0; k < K; k++) {
            if (!c->obs[k] && !debug) continue;
            int r = n_rows++;
            const float *gl = c->g_lab + 3 * k;
            float mag = sqrtf(fmaf(gl[2], gl[2], fmaf(gl[1], gl[1], gl[0] * gl[0])));
            float n[3] = {gl[0] / mag, gl[1] / mag, gl[2] / mag};
            float om = s == 0 ? c->om1[k] : c->om2[k];
            float cw = cosf(om), sw = sinf(om);
            float Rz[9] = {cw, -sw, 0.0f, sw, cw, 0.0f, 0.0f, 0.0f, 1.0f};
            float F[9];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    F[3 * i + j] = fmaf(Rz[3 * i + 2], B[6 + j],
                                        fmaf(Rz[3 * i + 1], B[3 + j], Rz[3 * i] * B[j]));
            float ln[3];
            for (int i = 0; i < 3; i++)
                ln[i] = fmaf(F[3 * i + 2], n[2], fmaf(F[3 * i + 1], n[1], F[3 * i] * n[0]));
            float dot = fused_sum3(c->beam[0] * ln[0], c->beam[1] * ln[1], c->beam[2] * ln[2]);
            float two_dot = 2.0f * dot;
            float rf[3];
            for (int i = 0; i < 3; i++) rf[i] = c->beam[i] - two_dot * ln[i];
            float rn = sqrtf(fmaf(rf[2], rf[2], fmaf(rf[1], rf[1], rf[0] * rf[0])));
            float safe = rn > 0.0f ? rn : 1.0f;
            float ry = fabsf(rf[1]) / safe, rzv = fabsf(rf[2]) / safe;
            float eta = atan2f(ry, rzv);
            unsigned char valid = (eta < c->eta_limit) && (rn > 0.0f);
            if (c->use_min_sin_eta) valid = valid && (sinf(eta) >= c->min_sin_eta);
            if (!c->obs[k]) valid = 0;
            int fr = valid ? fused_frame(c, om) : -1;
            c->row_om[r] = om;
            c->row_c[r] = cw;
            c->row_s[r] = sw;
            c->row_eta[r] = eta;
            c->row_valid[r] = valid;
            c->row_frame[r] = fr;
            for (int i = 0; i < 3; i++) {
                c->row_n[3 * r + i] = ln[i];
                c->row_refl[3 * r + i] = rf[i];
            }

            /* 3. stages B and C for this row (torch skips rows without a frame; computed for
               them only in debug mode) */
            if (fr < 0 && !debug) continue;
            float LV[3][3]; /* lab vertices: einsum('mij,vj->mvi', full_4x4, verts_4d) */
            for (int v = 0; v < 3; v++)
                for (int i = 0; i < 3; i++)
                    LV[v][i] = fmaf(c->trans[i], 1.0f,
                                    fmaf(F[3 * i + 2], V[3 * v + 2],
                                         fmaf(F[3 * i + 1], V[3 * v + 1], F[3 * i] * V[3 * v])));
            for (int d = 0; d < n_det; d++) {
                const float *pn = c->det_n[d];
                float denom = fused_sum3(rf[0] * pn[0], rf[1] * pn[1], rf[2] * pn[2]);
                int parallel = fabsf(denom) < 1e-8f;
                int all_hit = 1;
                for (int v = 0; v < 3; v++) {
                    const float *o = LV[v];
                    float numer =
                        -(fused_sum3(o[0] * pn[0], o[1] * pn[1], o[2] * pn[2]) + c->det_d[d]);
                    float t = parallel ? 0.0f : numer / denom;
                    int h = !parallel && (t > 0.0f);
                    all_hit &= h;
                    float pl[3];
                    for (int i = 0; i < 3; i++) {
                        float pt = o[i] + t * rf[i];
                        float rel = pt - c->det_pos[d][i];
                        pl[i] = rel - c->det_org[d][i];
                    }
                    float j = fused_sum3(pl[0] * c->det_bj[d][0], pl[1] * c->det_bj[d][1],
                                         pl[2] * c->det_bj[d][2]);
                    float kk = fused_sum3(pl[0] * c->det_bk[d][0], pl[1] * c->det_bk[d][1],
                                          pl[2] * c->det_bk[d][2]);
                    float col = (j + c->det_hw[d]) / c->det_pw[d];
                    float row = (kk + c->det_hh[d]) / c->det_ph[d];
                    if (c->rule_bits & 2) { /* CppEdgeRules: d < 0 is pixel -1 */
                        if (j < 0.0f) col = -1.0f;
                        if (kk < 0.0f) row = -1.0f;
                    }
                    size_t bi = (((size_t)r * n_det + d) * 3 + v);
                    c->pix[2 * bi] = col;
                    c->pix[2 * bi + 1] = row;
                    c->tval[bi] = t;
                }
                c->hit[(size_t)r * n_det + d] = (unsigned char)all_hit;
            }
        }
    }

    /* frames to fill first */
    PyObject *missing = NULL;
    for (int r = 0; r < n_rows; r++) {
        int fr = c->row_frame[r];
        if (fr < 0) continue;
        if (fr >= c->n_frames) {
            PyErr_SetString(PyExc_IndexError, "frame index beyond the FrameTable");
            Py_XDECREF(missing);
            return NULL;
        }
        if (!c->filled[fr]) {
            if (!missing && !(missing = PyList_New(0))) return NULL;
            PyObject *v = PyLong_FromLong(fr);
            if (!v || PyList_Append(missing, v) < 0) {
                Py_XDECREF(v);
                Py_DECREF(missing);
                return NULL;
            }
            Py_DECREF(v);
        }
    }
    if (missing && !debug) return missing;

    /* 4. stage D (icenine_stage_d_peak, as stage_d_overlap), rows in order */
    IceStageDAcc acc = {0, 0, 0, 0, 0, 0, 0.0};
    if (!missing) {
        PyObject *ro = PyObject_GetAttrString(c->table, "rows");
        PyObject *co = PyObject_GetAttrString(c->table, "cols");
        if (!ro || !co) {
            Py_XDECREF(ro);
            Py_XDECREF(co);
            return NULL;
        }
        long image_rows = PyLong_AsLong(ro), image_cols = PyLong_AsLong(co);
        Py_DECREF(ro);
        Py_DECREF(co);
        if (PyErr_Occurred()) return NULL;
        /* stage D reads every frame as image_rows x image_cols; FrameTable checks the size of
           detector 0 only, so each frame's shape is checked here */
        const npy_intp want_rows = (npy_intp)image_rows, want_cols = (npy_intp)image_cols;
        Py_ssize_t n_images = PyList_GET_SIZE(c->images);
        const unsigned char *peak_images[FUSED_MAX_DET];
        npy_int32 centers[2 * FUSED_MAX_DET];
        for (int r = 0; r < n_rows; r++) {
            int wedge_idx = c->row_frame[r];
            if (wedge_idx < 0) continue;
            const unsigned char *hit = c->hit + (size_t)r * n_det;
            const float *px = c->pix + 6 * (size_t)r * n_det;
            for (int d = 0; d < n_det; d++) {
                peak_images[d] = NULL;
                Py_ssize_t img_idx = (Py_ssize_t)wedge_idx * n_det + d;
                if (!hit[d] || img_idx < 0 || img_idx >= n_images) continue;
                PyObject *item = PyList_GET_ITEM(c->images, img_idx);
                if (!PyArray_Check(item) || PyArray_TYPE((PyArrayObject *)item) != NPY_UINT8 ||
                    PyArray_NDIM((PyArrayObject *)item) != 2 ||
                    !PyArray_IS_C_CONTIGUOUS((PyArrayObject *)item)) {
                    PyErr_SetString(PyExc_TypeError,
                                    "frame images must be C-contiguous 2D uint8 arrays");
                    return NULL;
                }
                if (PyArray_DIM((PyArrayObject *)item, 0) != want_rows ||
                    PyArray_DIM((PyArrayObject *)item, 1) != want_cols) {
                    PyErr_Format(PyExc_ValueError,
                                 "frame %zd (detector %d) has shape (%zd, %zd), expected the "
                                 "FrameTable size (%zd, %zd)",
                                 img_idx, d, (Py_ssize_t)PyArray_DIM((PyArrayObject *)item, 0),
                                 (Py_ssize_t)PyArray_DIM((PyArrayObject *)item, 1),
                                 (Py_ssize_t)want_rows, (Py_ssize_t)want_cols);
                    return NULL;
                }
                peak_images[d] = (const unsigned char *)PyArray_DATA((PyArrayObject *)item);
                /* pixel radius: vertex 0 truncated to int, as .astype(np.int32) */
                centers[2 * d] = (npy_int32)px[6 * d];
                centers[2 * d + 1] = (npy_int32)px[6 * d + 1];
            }
            icenine_stage_d_peak(peak_images, image_rows, image_cols, hit,
                                 c->pixel_radius > 0 ? centers : NULL,
                                 c->pixel_radius > 0 ? NULL : px, n_det, c->pixel_radius,
                                 c->rule_bits, NULL, NULL, NULL, &acc);
        }
    }
    if (!debug) {
        return Py_BuildValue("(iiiidii)", acc.pixel_overlap, acc.pixel_on_det, acc.peak_overlap,
                             acc.peak_on_det, acc.quality, acc.n_quality_points,
                             acc.last_n_det_ovlp);
    }

    /* debug: the intermediates of every (solution, reflection) row, row r = s * K + k */
    npy_intp dK[2] = {K, 3}, dN1[1] = {n_rows}, dpix[4] = {n_rows, n_det, 3, 2},
             dhit[2] = {n_rows, n_det}, dt[3] = {n_rows, n_det, 3}, dN3[2] = {n_rows, 3},
             dKo[1] = {K};
    PyObject *out = PyDict_New();
    if (!out) {
        Py_XDECREF(missing);
        return NULL;
    }
    int ok = 1;
#define FUSED_PUT(name, nd, dims, typ, src, bytes)                                                 \
    do {                                                                                           \
        PyObject *a_ = ok ? PyArray_SimpleNew(nd, dims, typ) : NULL;                               \
        if (!a_) {                                                                                 \
            ok = 0;                                                                                \
            break;                                                                                 \
        }                                                                                          \
        memcpy(PyArray_DATA((PyArrayObject *)a_), src, bytes);                                     \
        if (PyDict_SetItemString(out, name, a_) < 0) ok = 0;                                       \
        Py_DECREF(a_);                                                                             \
    } while (0)
    FUSED_PUT("g_lab", 2, dK, NPY_FLOAT32, c->g_lab, 3 * (size_t)K * sizeof(float));
    FUSED_PUT("obs", 1, dKo, NPY_UINT8, c->obs, (size_t)K);
    FUSED_PUT("omega1", 1, dKo, NPY_FLOAT32, c->om1, (size_t)K * sizeof(float));
    FUSED_PUT("omega2", 1, dKo, NPY_FLOAT32, c->om2, (size_t)K * sizeof(float));
    FUSED_PUT("omega", 1, dN1, NPY_FLOAT32, c->row_om, (size_t)n_rows * sizeof(float));
    FUSED_PUT("cos", 1, dN1, NPY_FLOAT32, c->row_c, (size_t)n_rows * sizeof(float));
    FUSED_PUT("sin", 1, dN1, NPY_FLOAT32, c->row_s, (size_t)n_rows * sizeof(float));
    FUSED_PUT("eta", 1, dN1, NPY_FLOAT32, c->row_eta, (size_t)n_rows * sizeof(float));
    FUSED_PUT("valid", 1, dN1, NPY_UINT8, c->row_valid, (size_t)n_rows);
    FUSED_PUT("frame", 1, dN1, NPY_INT32, c->row_frame, (size_t)n_rows * sizeof(int));
    FUSED_PUT("lab_normal", 2, dN3, NPY_FLOAT32, c->row_n, 3 * (size_t)n_rows * sizeof(float));
    FUSED_PUT("reflected", 2, dN3, NPY_FLOAT32, c->row_refl, 3 * (size_t)n_rows * sizeof(float));
    FUSED_PUT("pix", 4, dpix, NPY_FLOAT32, c->pix, (size_t)n_rows * n_det * 6 * sizeof(float));
    FUSED_PUT("hit", 2, dhit, NPY_UINT8, c->hit, (size_t)n_rows * n_det);
    FUSED_PUT("t", 3, dt, NPY_FLOAT32, c->tval, (size_t)n_rows * n_det * 3 * sizeof(float));
#undef FUSED_PUT
    if (ok) {
        if (missing) {
            ok = PyDict_SetItemString(out, "result", Py_None) == 0 &&
                 PyDict_SetItemString(out, "missing", missing) == 0;
        } else {
            PyObject *res = Py_BuildValue("(iiiidii)", acc.pixel_overlap, acc.pixel_on_det,
                                          acc.peak_overlap, acc.peak_on_det, acc.quality,
                                          acc.n_quality_points, acc.last_n_det_ovlp);
            ok = res && PyDict_SetItemString(out, "result", res) == 0;
            Py_XDECREF(res);
        }
    }
    Py_XDECREF(missing);
    if (!ok) {
        Py_DECREF(out);
        return NULL;
    }
    return out;
}
