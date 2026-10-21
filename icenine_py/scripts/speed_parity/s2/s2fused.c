/**
 * s2fused.c -- S2 Step 0 scratch extension (docs/todo_s2_fused_c_eval.md): one C call for
 * VoxelCostFunction.evaluate in hard mode without a ray rule (Bragg solve, observable rows, eta
 * filter, frame lookup, stages B and C, stage D).
 *
 * Scripts only: built by s2fused.py with the system cc into scripts/speed_parity/cache/s2/. Stage
 * D is not copied: the library's _rasterize.c is included, so triangle_overlap_uint8,
 * pixel_radius_overlap_uint8 and count_qualified_peaks_c are the library's own functions (its
 * module init is compiled in but unused).
 *
 * Float32 throughout, in the operation order of the torch code (cost_functions.py evaluate and
 * calculate_diffraction_overlap_batched; diffraction_core.py scattering_omegas_fixed_beam). Where
 * torch's kernels fuse a multiply and an add (bmm, the Rz @ base product, the einsum of the
 * vertices, torch.norm) the code calls fmaf explicitly, in the order found to match torch bit for
 * bit (scripts/speed_parity/s2/README: probe); FP_CONTRACT is off (the pragma in _rasterize.c), so
 * the compiler fuses nothing else. Not reproduced bit for bit: the 3 x K product R @ g_hkl^T
 * (Accelerate sgemm, no simple formula found) and sin, cos, atan2, asin (torch: SLEEF vector
 * kernels; here: libm). Differences are measured by s2_parity.py.
 *
 * Not thread-safe: the scratch buffers live in the context.
 */

#include "../../../icenine/_rasterize.c"

#define S2_MAX_DET MAX_DETECTORS

#ifdef S2_ACCELERATE
/* Accelerate's CBLAS (linked with -framework Accelerate; the enums are ints in the ABI) */
void cblas_sgemm(int order, int transa, int transb, int m, int n, int k, float alpha,
                 const float *a, int lda, const float *b, int ldb, float beta, float *c, int ldc);
#define S2_MODNAME "_s2fused_acc"
#define S2_INIT PyInit__s2fused_acc
#else
#define S2_MODNAME "_s2fused"
#define S2_INIT PyInit__s2fused
#endif

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
    int rule_bits; /* 1 cpp rule, 2 cpp edge rules */
    /* frame lookup: f = (omega - low) / width in the dtypes numpy uses (see s2fused.py) */
    double low, width;
    float low_f, width_f;
    int frame_mode; /* 0: sub and div in float64; 1: sub float32, div float64; 2: all float32 */
    const long long *frame_index;
    long long n_index;
    /* detectors */
    int n_det;
    float det_n[S2_MAX_DET][3], det_d[S2_MAX_DET];
    float det_pos[S2_MAX_DET][3], det_org[S2_MAX_DET][3];
    float det_bj[S2_MAX_DET][3], det_bk[S2_MAX_DET][3];
    float det_hw[S2_MAX_DET], det_pw[S2_MAX_DET], det_hh[S2_MAX_DET], det_ph[S2_MAX_DET];
    /* frame table (FrameTable): images list, filled flags, the table object (rows, cols) */
    PyObject *table;
    PyObject *keep[5]; /* g_hkl, sin_theta, g_mag_eps, frame_index, filled */
    PyObject *images;
    const npy_bool *filled;
    long long n_frames;
    /* scratch, sized for N = 2K rows */
    float *g_lab;  /* (K, 3) */
    float *g_tmp;  /* (3, K): the sgemm output (S2_ACCELERATE) */
    unsigned char *obs;
    float *om1, *om2; /* (K,) */
    float *row_om, *row_n, *row_eta, *row_c, *row_s; /* (N,), (N, 3), (N,), (N,), (N,) */
    float *row_refl;                                  /* (N, 3) */
    int *row_k;                                       /* (N,) reflection, -1 if row not used */
    unsigned char *row_valid;
    int *row_frame; /* (N,) -1: none */
    float *pix;     /* (N, n_det, 3, 2) */
    unsigned char *hit; /* (N, n_det) all three vertices hit */
    float *tval;        /* (N, n_det, 3) */
} S2Ctx;

static void s2_free(PyObject *cap) {
    S2Ctx *c = (S2Ctx *)PyCapsule_GetPointer(cap, "s2fused.ctx");
    if (!c) return;
    Py_XDECREF(c->table);
    for (int i = 0; i < 5; i++) Py_XDECREF(c->keep[i]);
    Py_XDECREF(c->images);
    free(c->g_lab);
    free(c->g_tmp);
    free(c->obs);
    free(c->om1);
    free(c->om2);
    free(c->row_om);
    free(c->row_n);
    free(c->row_eta);
    free(c->row_c);
    free(c->row_s);
    free(c->row_refl);
    free(c->row_k);
    free(c->row_valid);
    free(c->row_frame);
    free(c->pix);
    free(c->hit);
    free(c->tval);
    free(c);
}

static int s2_f32(PyObject *o, int n, float *out, const char *name) {
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
 * make_context(g_hkl, sin_theta, g_mag_eps, scal, beam, base, trans, ints, frame_f, frame_index,
 *              det_geom, table, images, filled) -> capsule
 * g_hkl (K,3) float32, sin_theta and g_mag_eps (K,) float32: kept by reference (the caller keeps
 * them alive). scal float32 [sin_chi_laue, cos_chi_laue, eps, eta_limit, min_sin_eta]; ints
 * [pixel_radius, center_screen, rule_bits, use_min_sin_eta, frame_mode, n_det]; frame_f float64
 * [low, width]; frame_index int64; det_geom (n_det, 20) float32 n[3] d pos[3] org[3] bj[3] bk[3]
 * hw pw hh ph; filled: the FrameTable's bool array (kept by reference).
 */
static PyObject *make_context(PyObject *self, PyObject *args) {
    PyArrayObject *g_hkl, *sin_theta, *g_mag_eps, *frame_index, *filled;
    PyObject *scal, *beam, *base, *trans, *ints_o, *frame_f, *geom, *table, *images;
    if (!PyArg_ParseTuple(args, "O!O!O!OOOOOOO!OOO!O!", &PyArray_Type, &g_hkl, &PyArray_Type,
                          &sin_theta, &PyArray_Type, &g_mag_eps, &scal, &beam, &base, &trans,
                          &ints_o, &frame_f, &PyArray_Type, &frame_index, &geom, &table,
                          &PyList_Type, &images, &PyArray_Type, &filled))
        return NULL;
    int K = (int)PyArray_DIM(g_hkl, 0);
    if (PyArray_TYPE(g_hkl) != NPY_FLOAT32 || !PyArray_IS_C_CONTIGUOUS(g_hkl) ||
        PyArray_NDIM(g_hkl) != 2 || PyArray_DIM(g_hkl, 1) != 3 ||
        PyArray_TYPE(sin_theta) != NPY_FLOAT32 || !PyArray_IS_C_CONTIGUOUS(sin_theta) ||
        PyArray_SIZE(sin_theta) != K || PyArray_TYPE(g_mag_eps) != NPY_FLOAT32 ||
        !PyArray_IS_C_CONTIGUOUS(g_mag_eps) || PyArray_SIZE(g_mag_eps) != K ||
        PyArray_TYPE(frame_index) != NPY_INT64 || !PyArray_IS_C_CONTIGUOUS(frame_index) ||
        PyArray_TYPE(filled) != NPY_BOOL || !PyArray_IS_C_CONTIGUOUS(filled)) {
        PyErr_SetString(PyExc_TypeError, "make_context: bad array type or shape");
        return NULL;
    }
    S2Ctx *c = (S2Ctx *)calloc(1, sizeof(S2Ctx));
    if (!c) return PyErr_NoMemory();
    c->K = K;
    c->g_hkl = (const float *)PyArray_DATA(g_hkl);
    c->sin_theta = (const float *)PyArray_DATA(sin_theta);
    c->g_mag_eps = (const float *)PyArray_DATA(g_mag_eps);
    float sc[5];
    int iv[6];
    double ff[2];
    {
        PyArrayObject *a = (PyArrayObject *)PyArray_FROM_OTF(ints_o, NPY_INT32, NPY_ARRAY_IN_ARRAY);
        PyArrayObject *b = (PyArrayObject *)PyArray_FROM_OTF(frame_f, NPY_FLOAT64,
                                                             NPY_ARRAY_IN_ARRAY);
        if (!a || !b || PyArray_SIZE(a) != 6 || PyArray_SIZE(b) != 2) {
            Py_XDECREF(a);
            Py_XDECREF(b);
            free(c);
            if (!PyErr_Occurred()) PyErr_SetString(PyExc_ValueError, "ints (6) / frame_f (2)");
            return NULL;
        }
        memcpy(iv, PyArray_DATA(a), sizeof(iv));
        memcpy(ff, PyArray_DATA(b), sizeof(ff));
        Py_DECREF(a);
        Py_DECREF(b);
    }
    int n_det = iv[5];
    if (n_det < 1 || n_det > S2_MAX_DET) {
        free(c);
        PyErr_SetString(PyExc_ValueError, "n_det out of range");
        return NULL;
    }
    float g[S2_MAX_DET * 20];
    if (s2_f32(scal, 5, sc, "scal") || s2_f32(beam, 3, c->beam, "beam") ||
        s2_f32(base, 9, c->base, "base") || s2_f32(trans, 3, c->trans, "trans") ||
        s2_f32(geom, n_det * 20, g, "det_geom")) {
        free(c);
        return NULL;
    }
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
    /* the context owns a reference to every array it points into */
    c->keep[0] = (PyObject *)g_hkl;
    c->keep[1] = (PyObject *)sin_theta;
    c->keep[2] = (PyObject *)g_mag_eps;
    c->keep[3] = (PyObject *)frame_index;
    c->keep[4] = (PyObject *)filled;
    for (int i = 0; i < 5; i++) Py_INCREF(c->keep[i]);
    int N = 2 * K + 1;
    c->g_lab = (float *)malloc(3 * (K + 1) * sizeof(float));
    c->g_tmp = (float *)malloc(3 * (K + 1) * sizeof(float));
    c->obs = (unsigned char *)malloc(K + 1);
    c->om1 = (float *)malloc((K + 1) * sizeof(float));
    c->om2 = (float *)malloc((K + 1) * sizeof(float));
    c->row_om = (float *)malloc(N * sizeof(float));
    c->row_n = (float *)malloc(3 * N * sizeof(float));
    c->row_eta = (float *)malloc(N * sizeof(float));
    c->row_c = (float *)malloc(N * sizeof(float));
    c->row_s = (float *)malloc(N * sizeof(float));
    c->row_refl = (float *)malloc(3 * N * sizeof(float));
    c->row_k = (int *)malloc(N * sizeof(int));
    c->row_valid = (unsigned char *)malloc(N);
    c->row_frame = (int *)malloc(N * sizeof(int));
    c->pix = (float *)malloc((size_t)N * n_det * 6 * sizeof(float));
    c->hit = (unsigned char *)malloc((size_t)N * n_det);
    c->tval = (float *)malloc((size_t)N * n_det * 3 * sizeof(float));
    PyObject *cap = PyCapsule_New(c, "s2fused.ctx", s2_free);
    if (!cap) {
        Py_DECREF(table);
        Py_DECREF(images);
        for (int i = 0; i < 5; i++) Py_DECREF(c->keep[i]);
        free(c);
        return NULL;
    }
    if (!c->g_lab || !c->g_tmp || !c->obs || !c->om1 || !c->om2 || !c->row_om || !c->row_n ||
        !c->row_eta ||
        !c->row_c || !c->row_s || !c->row_refl || !c->row_k || !c->row_valid || !c->row_frame ||
        !c->pix || !c->hit || !c->tval) {
        Py_DECREF(cap);
        return PyErr_NoMemory();
    }
    return cap;
}

/* frame of an omega: _omega_bins + frame_of_omegas (cost_functions.py) */
static int s2_frame(const S2Ctx *c, float om) {
    double f;
    if (c->frame_mode == 0) {
        f = ((double)om - c->low) / c->width;
    } else if (c->frame_mode == 1) {
        f = (double)(om - c->low_f) / c->width;
    } else {
        f = (double)((om - c->low_f) / c->width_f);
    }
    long long b;
    if (c->rule_bits & 2) {
        f = floor(f);
    }
    /* astype(int): truncation toward zero; values are small (|omega| <= pi) */
    if (!(f > -1e15)) b = -1000000000000000LL;
    else if (f > 1e15) b = 1000000000000000LL;
    else b = (long long)f;
    if (b < -1) b = -1;
    if (b > c->n_index - 1) b = c->n_index - 1;
    long long idx = b < 0 ? c->n_index + b : b; /* index -1 is the last (-1) entry */
    return (int)c->frame_index[idx];
}

static inline float s2_sum3(float a, float b, float d) { return (a + b) + d; }

/**
 * evaluate(ctx, R, verts, debug) -> (po, pon, pko, pkon, quality, nq, det_ovlp)
 *   or a list of frame indices that the FrameTable must fill first (call again after ensure);
 *   with debug=1 a dict of the intermediates instead (no stage D unless all frames are filled).
 * R: (3, 3) array (any float dtype; cast to float64 then to float32 as torch .float()).
 * verts: (3, 3) float32 C-contiguous.
 */
static PyObject *s2_evaluate(PyObject *self, PyObject *args) {
    PyObject *cap, *R_o;
    PyArrayObject *V_a;
    int debug = 0;
    if (!PyArg_ParseTuple(args, "OOO!|i", &cap, &R_o, &PyArray_Type, &V_a, &debug)) return NULL;
    S2Ctx *c = (S2Ctx *)PyCapsule_GetPointer(cap, "s2fused.ctx");
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
            float m = s2_sum3(V0[i], V0[3 + i], V0[6 + i]) / 3.0f;
            V[i] = V[3 + i] = V[6 + i] = m;
        }
    } else {
        memcpy(V, V0, sizeof(V));
    }

    const int K = c->K, n_det = c->n_det;
    const float eps = c->eps;
    const float PI_F = (float)M_PI, TWO_PI_F = (float)(2.0 * M_PI);
    /* 1. g_lab = (R @ g^T)^T, and the Bragg solve (scattering_omegas_fixed_beam) */
#ifdef S2_ACCELERATE
    /* the call torch makes on macOS (BLAS_INFO=accelerate): bit-identical g_lab (probe) */
    cblas_sgemm(101 /* row major */, 111 /* no trans */, 112 /* trans */, 3, K, 3, 1.0f, Rf, 3,
                c->g_hkl, 3, 0.0f, c->g_tmp, K);
    for (int k = 0; k < K; k++)
        for (int i = 0; i < 3; i++) c->g_lab[3 * k + i] = c->g_tmp[i * K + k];
#endif
    for (int k = 0; k < K; k++) {
        float *gl = c->g_lab + 3 * k;
#ifndef S2_ACCELERATE
        const float *g = c->g_hkl + 3 * k;
        for (int i = 0; i < 3; i++)
            gl[i] = fmaf(Rf[3 * i + 2], g[2], fmaf(Rf[3 * i + 1], g[1], Rf[3 * i] * g[0]));
#endif
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
       Rz @ base, lab normal, reflection, eta filter, frame */
    const float *B = c->base;
    int n_rows = 0;
    for (int s = 0; s < 2; s++) {
        for (int k = 0; k < K; k++) {
            if (!c->obs[k] && !debug) continue;
            int r = n_rows++;
            c->row_k[r] = c->obs[k] ? k : -1;
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
            float dot = s2_sum3(c->beam[0] * ln[0], c->beam[1] * ln[1], c->beam[2] * ln[2]);
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
            int fr = valid ? s2_frame(c, om) : -1;
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

            /* 3. stages B and C for this row (rows without a frame are skipped by torch; computed
               here only in debug mode) */
            if (fr < 0 && !debug) continue;
            float LV[3][3]; /* lab vertices: einsum('mij,vj->mvi', full_4x4, verts_4d) */
            for (int v = 0; v < 3; v++)
                for (int i = 0; i < 3; i++)
                    LV[v][i] = fmaf(c->trans[i], 1.0f,
                                    fmaf(F[3 * i + 2], V[3 * v + 2],
                                         fmaf(F[3 * i + 1], V[3 * v + 1], F[3 * i] * V[3 * v])));
            for (int d = 0; d < n_det; d++) {
                const float *pn = c->det_n[d];
                float denom = s2_sum3(rf[0] * pn[0], rf[1] * pn[1], rf[2] * pn[2]);
                int parallel = fabsf(denom) < 1e-8f;
                int all_hit = 1;
                for (int v = 0; v < 3; v++) {
                    const float *o = LV[v];
                    float numer = -(s2_sum3(o[0] * pn[0], o[1] * pn[1], o[2] * pn[2]) + c->det_d[d]);
                    float t = parallel ? 0.0f : numer / denom;
                    int h = !parallel && (t > 0.0f);
                    all_hit &= h;
                    float pl[3];
                    for (int i = 0; i < 3; i++) {
                        float pt = o[i] + t * rf[i];
                        float rel = pt - c->det_pos[d][i];
                        pl[i] = rel - c->det_org[d][i];
                    }
                    float j = s2_sum3(pl[0] * c->det_bj[d][0], pl[1] * c->det_bj[d][1],
                                      pl[2] * c->det_bj[d][2]);
                    float kk = s2_sum3(pl[0] * c->det_bk[d][0], pl[1] * c->det_bk[d][1],
                                       pl[2] * c->det_bk[d][2]);
                    float col = (j + c->det_hw[d]) / c->det_pw[d];
                    float row = (kk + c->det_hh[d]) / c->det_ph[d];
                    if (c->rule_bits & 2) {
                        if (j < 0.0f) col = -1.0f;
                        if (kk < 0.0f) row = -1.0f;
                    }
                    size_t base_i = (((size_t)r * n_det + d) * 3 + v);
                    c->pix[2 * base_i] = col;
                    c->pix[2 * base_i + 1] = row;
                    c->tval[base_i] = t;
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
            if (!missing) missing = PyList_New(0);
            PyObject *v = PyLong_FromLong(fr);
            PyList_Append(missing, v);
            Py_DECREF(v);
        }
    }
    if (missing && !debug) return missing;

    /* 4. stage D: the per-peak body of stage_d_overlap (_rasterize.c), rows in order */
    int total_pixel_overlap = 0, total_pixel_on_det = 0, total_peak_overlap = 0,
        total_peak_on_det = 0, n_quality_points = 0, last_n_det_ovlp = 0;
    double quality = 0.0;
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
        Py_ssize_t n_images = PyList_GET_SIZE(c->images);
        const int cpp_vertex_lit = (c->rule_bits & 2) != 0;
        for (int r = 0; r < n_rows; r++) {
            int wedge_idx = c->row_frame[r];
            if (wedge_idx < 0) continue;
            int detector_lit[MAX_DETECTORS];
            int spot_overlap[MAX_DETECTORS];
            int peak_pixel_overlap = 0, peak_pixel_on_det = 0;
            memset(detector_lit, 0, n_det * sizeof(int));
            memset(spot_overlap, 0, n_det * sizeof(int));
            for (int d = 0; d < n_det; d++) {
                if (!c->hit[(size_t)r * n_det + d]) continue;
                Py_ssize_t img_idx = (Py_ssize_t)wedge_idx * n_det + d;
                if (img_idx < 0 || img_idx >= n_images) continue;
                PyObject *item = PyList_GET_ITEM(c->images, img_idx);
                if (!PyArray_Check(item) || PyArray_TYPE((PyArrayObject *)item) != NPY_UINT8 ||
                    !PyArray_IS_C_CONTIGUOUS((PyArrayObject *)item)) {
                    PyErr_SetString(PyExc_TypeError, "frame images must be C-contiguous uint8");
                    return NULL;
                }
                const unsigned char *image_data =
                    (const unsigned char *)PyArray_DATA((PyArrayObject *)item);
                const float *px = c->pix + 6 * ((size_t)r * n_det + d);
                if (c->pixel_radius > 0) {
                    int cx = (int)px[0];
                    int cy = (int)px[1];
                    int in_bounds, bright;
                    pixel_radius_overlap_uint8(image_data, image_rows, image_cols, cx, cy,
                                               c->pixel_radius, &in_bounds, &bright);
                    if (in_bounds) {
                        detector_lit[d] = 1;
                        peak_pixel_on_det++;
                    }
                    if (bright) {
                        peak_pixel_overlap++;
                        spot_overlap[d] = 1;
                    }
                } else {
                    double v0x = (double)px[0], v0y = (double)px[1];
                    double v1x = (double)px[2], v1y = (double)px[3];
                    double v2x = (double)px[4], v2y = (double)px[5];
                    int n_overlap, n_lit;
                    triangle_overlap_uint8(image_data, image_rows, image_cols, v0x, v0y, v1x,
                                           v1y, v2x, v2y, &n_overlap, &n_lit);
                    peak_pixel_overlap += n_overlap;
                    peak_pixel_on_det += n_lit;
                    if (n_overlap > 0) {
                        detector_lit[d] = 1;
                        spot_overlap[d] = 1;
                    } else if (cpp_vertex_lit) {
                        double vx[3] = {v0x, v1x, v2x}, vy[3] = {v0y, v1y, v2y};
                        for (int q = 0; q < 3; q++) {
                            double qx = (vx[q] < 0) ? -1.0 : floor(vx[q]);
                            double qy = (vy[q] < 0) ? -1.0 : floor(vy[q]);
                            if (qx >= 0 && qx < image_cols && qy >= 0 && qy < image_rows) {
                                detector_lit[d] = 1;
                                break;
                            }
                        }
                    } else if (n_lit > 0) {
                        detector_lit[d] = 1;
                    }
                }
            }
            int peak_on_det, peak_ovlp, n_det_ovlp;
            count_qualified_peaks_c(detector_lit, spot_overlap, n_det, &peak_on_det, &peak_ovlp,
                                    &n_det_ovlp, c->rule_bits & 1);
            last_n_det_ovlp = n_det_ovlp;
            total_pixel_overlap += peak_pixel_overlap;
            total_pixel_on_det += peak_pixel_on_det;
            total_peak_overlap += peak_ovlp;
            total_peak_on_det += peak_on_det;
            if (peak_pixel_on_det > 0 && n_det > 0) {
                double pixel_ratio = (double)peak_pixel_overlap / peak_pixel_on_det;
                double cur_quality;
                if (n_det_ovlp > 0) {
                    double det_ratio = (double)n_det_ovlp / n_det;
                    cur_quality = pixel_ratio * det_ratio;
                } else {
                    cur_quality = 0.0;
                }
                quality += (cur_quality - quality) / (n_quality_points + 1);
                n_quality_points++;
            }
        }
    }
    if (!debug) {
        return Py_BuildValue("(iiiidii)", total_pixel_overlap, total_pixel_on_det,
                             total_peak_overlap, total_peak_on_det, quality, n_quality_points,
                             last_n_det_ovlp);
    }

    /* debug: the intermediates of every (solution, reflection) row, row r = s * K + k */
    npy_intp dK[2] = {K, 3}, dN1[1] = {n_rows}, dpix[4] = {n_rows, n_det, 3, 2},
             dhit[2] = {n_rows, n_det}, dt[3] = {n_rows, n_det, 3}, dN3[2] = {n_rows, 3},
             dKo[1] = {K};
    PyObject *out = PyDict_New();
#define S2_PUT(name, nd, dims, typ, src, bytes)                                \
    do {                                                                       \
        PyObject *a_ = PyArray_SimpleNew(nd, dims, typ);                       \
        memcpy(PyArray_DATA((PyArrayObject *)a_), src, bytes);                 \
        PyDict_SetItemString(out, name, a_);                                   \
        Py_DECREF(a_);                                                         \
    } while (0)
    S2_PUT("g_lab", 2, dK, NPY_FLOAT32, c->g_lab, 3 * K * sizeof(float));
    S2_PUT("obs", 1, dKo, NPY_UINT8, c->obs, K);
    S2_PUT("omega1", 1, dKo, NPY_FLOAT32, c->om1, K * sizeof(float));
    S2_PUT("omega2", 1, dKo, NPY_FLOAT32, c->om2, K * sizeof(float));
    S2_PUT("omega", 1, dN1, NPY_FLOAT32, c->row_om, n_rows * sizeof(float));
    S2_PUT("cos", 1, dN1, NPY_FLOAT32, c->row_c, n_rows * sizeof(float));
    S2_PUT("sin", 1, dN1, NPY_FLOAT32, c->row_s, n_rows * sizeof(float));
    S2_PUT("eta", 1, dN1, NPY_FLOAT32, c->row_eta, n_rows * sizeof(float));
    S2_PUT("valid", 1, dN1, NPY_UINT8, c->row_valid, n_rows);
    S2_PUT("frame", 1, dN1, NPY_INT32, c->row_frame, n_rows * sizeof(int));
    S2_PUT("lab_normal", 2, dN3, NPY_FLOAT32, c->row_n, 3 * n_rows * sizeof(float));
    S2_PUT("reflected", 2, dN3, NPY_FLOAT32, c->row_refl, 3 * n_rows * sizeof(float));
    S2_PUT("pix", 4, dpix, NPY_FLOAT32, c->pix, (size_t)n_rows * n_det * 6 * sizeof(float));
    S2_PUT("hit", 2, dhit, NPY_UINT8, c->hit, (size_t)n_rows * n_det);
    S2_PUT("t", 3, dt, NPY_FLOAT32, c->tval, (size_t)n_rows * n_det * 3 * sizeof(float));
#undef S2_PUT
    PyObject *res = Py_BuildValue("(iiiidii)", total_pixel_overlap, total_pixel_on_det,
                                  total_peak_overlap, total_peak_on_det, quality,
                                  n_quality_points, last_n_det_ovlp);
    PyDict_SetItemString(out, "result", missing ? Py_None : res);
    Py_DECREF(res);
    if (missing) {
        PyDict_SetItemString(out, "missing", missing);
        Py_DECREF(missing);
    }
    return out;
}

static PyMethodDef s2_methods[] = {
    {"make_context", make_context, METH_VARARGS, "build the per-cost-function context"},
    {"evaluate", s2_evaluate, METH_VARARGS, "fused evaluate (see s2fused.c)"},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef s2_module = {PyModuleDef_HEAD_INIT, S2_MODNAME,
                                       "S2 Step 0 scratch fused evaluate", -1, s2_methods};

PyMODINIT_FUNC S2_INIT(void) {
    import_array();
    return PyModule_Create(&s2_module);
}
