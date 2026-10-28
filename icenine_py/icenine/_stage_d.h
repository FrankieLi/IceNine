/**
 * _stage_d.h -- internal declarations shared by the two translation units of the _rasterize module
 * (_rasterize.c and _fused_eval.c, S2): the shared numpy C-API table, the detector limit, the
 * stage-D accumulator, the per-peak stage-D function and the fused entry points. One copy, so the
 * two files cannot disagree on a struct layout. Listed in setup.py `depends=` so that an edit
 * rebuilds both files.
 *
 * Include after <Python.h> and in place of <numpy/arrayobject.h> (this header includes it, after
 * declaring the shared API table hidden). _fused_eval.c defines NO_IMPORT_ARRAY first: the table
 * is defined and imported by the module init in _rasterize.c.
 *
 * The functions and the API table are not static (each is defined in one file and used from the
 * other) but hidden: the .so exports only PyInit__rasterize.
 */
#ifndef ICENINE_STAGE_D_H
#define ICENINE_STAGE_D_H

#if defined(__GNUC__) || defined(__clang__)
#define ICENINE_HIDDEN __attribute__((visibility("hidden")))
#else
#define ICENINE_HIDDEN
#endif

/* one numpy C-API table for both files; numpy 2 adds the RUNTIME_VERSION companion (an unused
   extern declaration under numpy 1) */
#define PY_ARRAY_UNIQUE_SYMBOL icenine_rasterize_ARRAY_API
extern ICENINE_HIDDEN void **icenine_rasterize_ARRAY_API;
extern ICENINE_HIDDEN int icenine_rasterize_ARRAY_APIPyArray_RUNTIME_VERSION;
#include <numpy/arrayobject.h>

/* Maximum detectors of stage_d_overlap, ray_flags and the fused evaluate
   (cost_functions._FUSED_MAX_DETECTORS) */
#define MAX_DETECTORS 8

/* Stage-D accumulator of one evaluation (the OverlapInfo fields) */
typedef struct {
    int pixel_overlap, pixel_on_det, peak_overlap, peak_on_det, n_quality_points, last_n_det_ovlp;
    double quality;
} IceStageDAcc;

/* _rasterize.c: the per-peak body of stage D (see its definition for the arguments) */
ICENINE_HIDDEN void icenine_stage_d_peak(const unsigned char *const *images, npy_intp image_rows,
                                         npy_intp image_cols, const npy_uint8 *det_hit,
                                         const npy_int32 *centers, const float *vertices,
                                         int n_detectors, int pixel_radius, int rule_bits,
                                         const npy_uint8 *obs_in, const npy_uint8 *obs_bright,
                                         long *cases, IceStageDAcc *acc);

#ifdef ICENINE_FUSED
/* _fused_eval.c: the fused evaluate (S2), in the method table of _rasterize.c */
ICENINE_HIDDEN PyObject *icenine_fused_context(PyObject *self, PyObject *args);
ICENINE_HIDDEN PyObject *icenine_fused_evaluate(PyObject *self, PyObject *args);
#endif

#endif /* ICENINE_STAGE_D_H */
