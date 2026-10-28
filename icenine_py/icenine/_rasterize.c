/**
 * _rasterize.c — CPython C extension for fast triangle rasterization + overlap.
 *
 * Implements the hot path of cost function evaluation:
 *   triangle_overlap(image, v0x, v0y, v1x, v1y, v2x, v2y) -> (overlap, lit)
 *   pixel_radius_overlap(image, cx, cy, radius) -> (in_bounds, bright)
 *   stage_d_overlap(...) -> (pixel_overlap, pixel_on_det, peak_overlap,
 *                            peak_on_det, quality, n_quality_points)
 *
 * Ports the Python algorithms from image_data.py:
 *   _sutherland_hodgman_clip + _scanline_fill + _bresenham_edge + overlap count
 *
 * C++ Reference: Src/XDMRaster.cpp GetTriangleOverlapProperty,
 *                Src/Raster.tmpl.cpp GeneralRasterizePolygon + CalculateScanline
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
/* the module has a second translation unit (_fused_eval.c, built with ICENINE_FUSED): both share
   one numpy C-API table, imported here in the module init */
#define PY_ARRAY_UNIQUE_SYMBOL icenine_rasterize_ARRAY_API
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Maximum polygon vertices after clipping (triangle clipped by 4 edges -> max 7 vertices) */
#define MAX_POLY_VERTS 16
/* Maximum scanlines for edge table */
#define MAX_SCANLINES 4096
/* Maximum detectors */
#define MAX_DETECTORS 8

/* Stage-D accumulator of one evaluation (the OverlapInfo fields). The same typedef is in
   _fused_eval.c; keep the two identical. */
typedef struct {
    int pixel_overlap, pixel_on_det, peak_overlap, peak_on_det, n_quality_points, last_n_det_ovlp;
    double quality;
} IceStageDAcc;

/* ---- Sutherland-Hodgman polygon clipping ---- */

typedef struct {
    double x, y;
} Point2D;

static int clip_edge(
    const Point2D *in, int n_in,
    Point2D *out,
    int (*is_inside)(double, double, double),
    void (*intersect)(const Point2D*, const Point2D*, double, Point2D*),
    double boundary
) {
    if (n_in == 0) return 0;
    int n_out = 0;
    Point2D prev = in[n_in - 1];
    int prev_inside = is_inside(prev.x, prev.y, boundary);

    for (int i = 0; i < n_in; i++) {
        Point2D curr = in[i];
        int curr_inside = is_inside(curr.x, curr.y, boundary);
        if (curr_inside) {
            if (!prev_inside) {
                intersect(&prev, &curr, boundary, &out[n_out++]);
            }
            out[n_out++] = curr;
        } else if (prev_inside) {
            intersect(&prev, &curr, boundary, &out[n_out++]);
        }
        prev = curr;
        prev_inside = curr_inside;
        if (n_out >= MAX_POLY_VERTS - 1) break;
    }
    return n_out;
}

/* Boundary tests — matching C++ SutherlandHodgman.h conventions */
static int inside_left(double x, double y, double b) { (void)y; return x >= b; }
static int inside_right(double x, double y, double b) { (void)y; return x < b; }
static int inside_top(double x, double y, double b) { (void)x; return y >= b; }
static int inside_bottom(double x, double y, double b) { (void)x; return y < b; }

static void intersect_left(const Point2D *p0, const Point2D *p1, double b, Point2D *out) {
    double dx = p1->x - p0->x;
    if (fabs(dx) < 0.01) { out->x = b; out->y = p0->y; return; }
    double slope = (p1->y - p0->y) / dx;
    out->x = b;
    out->y = p0->y + slope * (b - p0->x);
}

static void intersect_right(const Point2D *p0, const Point2D *p1, double b, Point2D *out) {
    double dx = p1->x - p0->x;
    if (fabs(dx) < 0.01) { out->x = b; out->y = p0->y; return; }
    double slope = (p1->y - p0->y) / dx;
    out->x = b;
    out->y = p0->y + slope * (b - p0->x);
}

static void intersect_top(const Point2D *p0, const Point2D *p1, double b, Point2D *out) {
    double dy = p1->y - p0->y;
    if (fabs(dy) < 0.01) { out->x = p0->x; out->y = b; return; }
    double slope = (p1->x - p0->x) / dy;
    out->x = p0->x + slope * (b - p0->y);
    out->y = b;
}

static void intersect_bottom(const Point2D *p0, const Point2D *p1, double b, Point2D *out) {
    double dy = p1->y - p0->y;
    if (fabs(dy) < 0.01) { out->x = p0->x; out->y = b; return; }
    double slope = (p1->x - p0->x) / dy;
    out->x = p0->x + slope * (b - p0->y);
    out->y = b;
}

static int sutherland_hodgman_clip(
    Point2D *polygon, int n,
    double x_min, double x_max, double y_min, double y_max
) {
    Point2D buf1[MAX_POLY_VERTS], buf2[MAX_POLY_VERTS];

    /* C++ clips in order: Right, Top, Left, Bottom */
    memcpy(buf1, polygon, n * sizeof(Point2D));
    n = clip_edge(buf1, n, buf2, inside_right, intersect_right, x_max);
    if (n < 3) return 0;
    n = clip_edge(buf2, n, buf1, inside_top, intersect_top, y_min);
    if (n < 3) return 0;
    n = clip_edge(buf1, n, buf2, inside_left, intersect_left, x_min);
    if (n < 3) return 0;
    n = clip_edge(buf2, n, buf1, inside_bottom, intersect_bottom, y_max);

    memcpy(polygon, buf1, n * sizeof(Point2D));
    return n;
}

/* ---- Bresenham edge tracing + scanline fill ---- */

typedef struct {
    int left;
    int right;
} ScanlineEntry;

static void bresenham_edge(ScanlineEntry *table, int y_offset, int table_size,
                           int v0x, int v0y, int v1x, int v1y) {
    int steep = abs(v1y - v0y) > abs(v1x - v0x);
    if (steep) {
        int tmp;
        tmp = v0x; v0x = v0y; v0y = tmp;
        tmp = v1x; v1x = v1y; v1y = tmp;
    }
    if (v0x > v1x) {
        int tmp;
        tmp = v0x; v0x = v1x; v1x = tmp;
        tmp = v0y; v0y = v1y; v1y = tmp;
    }

    int delta_x = v1x - v0x;
    int delta_y = abs(v1y - v0y);
    int error = delta_x;
    int y_step = (v0y < v1y) ? 1 : -1;
    int y = v0y;

    for (int x = v0x; x <= v1x; x++) {
        int row, col;
        if (steep) { row = x; col = y; }
        else { row = y; col = x; }

        int idx = row - y_offset;
        if (idx >= 0 && idx < table_size) {
            if (table[idx].left < 0) {
                table[idx].left = col;
            } else if (table[idx].left > col) {
                table[idx].left = col;
            }
            if (table[idx].right < col) {
                table[idx].right = col;
            }
        }

        error -= 2 * delta_y;
        if (error < 0) {
            y += y_step;
            error += 2 * delta_x;
        }
    }
}

/* ---- Internal helpers for batch operations ---- */

/**
 * Triangle overlap on uint8 binary image. Returns (overlap, lit).
 * image_data: C-contiguous uint8 array, row-major
 */
static void triangle_overlap_uint8(
    const unsigned char *image_data, npy_intp num_rows, npy_intp num_cols,
    double v0x, double v0y, double v1x, double v1y, double v2x, double v2y,
    int *out_overlap, int *out_lit
) {
    *out_overlap = 0;
    *out_lit = 0;

    /* Truncate pixel coordinates */
    Point2D polygon[MAX_POLY_VERTS];
    polygon[0].x = (v0x < 0) ? -1.0 : floor(v0x);
    polygon[0].y = (v0y < 0) ? -1.0 : floor(v0y);
    polygon[1].x = (v1x < 0) ? -1.0 : floor(v1x);
    polygon[1].y = (v1y < 0) ? -1.0 : floor(v1y);
    polygon[2].x = (v2x < 0) ? -1.0 : floor(v2x);
    polygon[2].y = (v2y < 0) ? -1.0 : floor(v2y);

    int n_verts = sutherland_hodgman_clip(
        polygon, 3,
        0.0, (double)(num_cols - 1), 0.0, (double)(num_rows - 1)
    );

    if (n_verts < 3) return;

    int int_verts[MAX_POLY_VERTS][2];
    int y_min = INT_MAX, y_max = INT_MIN;
    for (int i = 0; i < n_verts; i++) {
        int_verts[i][0] = (int)round(polygon[i].x);
        int_verts[i][1] = (int)round(polygon[i].y);
        if (int_verts[i][1] < y_min) y_min = int_verts[i][1];
        if (int_verts[i][1] > y_max) y_max = int_verts[i][1];
    }

    int n_scanlines = y_max - y_min + 1;
    if (n_scanlines <= 0 || n_scanlines > MAX_SCANLINES) return;

    ScanlineEntry *table = (ScanlineEntry*)malloc(n_scanlines * sizeof(ScanlineEntry));
    if (!table) return;
    for (int i = 0; i < n_scanlines; i++) {
        table[i].left = -1;
        table[i].right = -1;
    }

    for (int i = 0; i < n_verts; i++) {
        int prev = (i == 0) ? n_verts - 1 : i - 1;
        bresenham_edge(table, y_min, n_scanlines,
                       int_verts[prev][0], int_verts[prev][1],
                       int_verts[i][0], int_verts[i][1]);
    }

    /* Degenerate case: horizontal line */
    if (y_min == y_max) {
        int x_min_h = INT_MAX, x_max_h = INT_MIN;
        for (int i = 0; i < n_verts; i++) {
            if (int_verts[i][0] < x_min_h) x_min_h = int_verts[i][0];
            if (int_verts[i][0] > x_max_h) x_max_h = int_verts[i][0];
        }
        for (int x = x_min_h; x <= x_max_h; x++) {
            if (x >= 0 && x < num_cols && y_min >= 0 && y_min < num_rows) {
                (*out_lit)++;
                if (image_data[y_min * num_cols + x]) (*out_overlap)++;
            }
        }
        free(table);
        return;
    }

    for (int i = 0; i < n_scanlines; i++) {
        int left = table[i].left;
        int right = table[i].right;
        int y = y_min + i;

        if (left < 0 && right < 0) continue;
        if (left < 0) left = right;
        if (right < 0) right = left;
        if (left > right) { int tmp = left; left = right; right = tmp; }

        for (int x = left; x <= right; x++) {
            if (x >= 0 && x < num_cols && y >= 0 && y < num_rows) {
                (*out_lit)++;
                if (image_data[y * num_cols + x]) (*out_overlap)++;
            }
        }
    }

    free(table);
}

/**
 * Pixel-radius overlap on uint8 binary image. Returns (in_bounds, bright).
 */
static void pixel_radius_overlap_uint8(
    const unsigned char *image_data, npy_intp num_rows, npy_intp num_cols,
    int cx, int cy, int radius,
    int *out_in_bounds, int *out_bright
) {
    *out_in_bounds = 0;
    *out_bright = 0;

    for (int dy = -radius; dy <= radius && !(*out_bright); dy++) {
        for (int dx = -radius; dx <= radius; dx++) {
            int px = cx + dx;
            int py = cy + dy;
            if (px >= 0 && px < num_cols && py >= 0 && py < num_rows) {
                *out_in_bounds = 1;
                if (image_data[py * num_cols + px]) {
                    *out_bright = 1;
                    break;
                }
            }
        }
    }
}

/**
 * count_qualified_peaks — C port of Python count_qualified_peaks().
 *
 * Check if a peak has valid (contiguous) detector coverage.
 * A peak is "qualified" as on-detector if it lights up a contiguous set
 * of detectors starting from detector 0 (allowing trailing gap).
 *
 * C++ Reference: OverlapInfo.tmpl.cpp CountQualifiedPeaks
 */
static void count_qualified_peaks_c(
    const int *detector_lit, const int *spot_overlap, int n_detectors,
    int *out_peak_on_det, int *out_peak_overlap, int *out_n_det_overlap,
    int cpp_rule
) {
    *out_peak_on_det = 0;
    *out_peak_overlap = 0;
    *out_n_det_overlap = 0;

    if (n_detectors == 0) return;

    if (cpp_rule) {
        /* Literal port of C++ CountQualifiedPeaks (OverlapInfo.tmpl.cpp:78-139) */
        int bValidSimPeak = detector_lit[0];
        int bOnPreviousDetector = detector_lit[0];
        int nDetectorsOverlap = 0;
        for (int i = 1; i < n_detectors && bValidSimPeak; i++) {
            if (detector_lit[i]) {
                bValidSimPeak = 1;
            } else if (bOnPreviousDetector) {
                bValidSimPeak = 1;
                bOnPreviousDetector = 0;
            } else {
                bValidSimPeak = 0;
            }
        }
        int bPreviousOverlap = spot_overlap[0];
        int bValidOverlap = spot_overlap[0];
        if (spot_overlap[0]) nDetectorsOverlap++;
        for (int i = 1; i < n_detectors && bValidOverlap; i++) {
            if (spot_overlap[i]) {
                nDetectorsOverlap++;
                bValidOverlap = 1;
            } else if (bPreviousOverlap && !detector_lit[i]) {
                bValidOverlap = 1;
                bPreviousOverlap = 0;
            } else {
                bValidOverlap = 0;
            }
        }
        if (bValidSimPeak) *out_peak_on_det = 1;
        if (bValidOverlap) *out_peak_overlap = 1;
        if (bValidSimPeak && bValidOverlap) *out_n_det_overlap = nDetectorsOverlap;
        return;
    }

    /* Check contiguous lit pattern (must start at detector 0) */
    int valid_sim_peak = 0;
    if (detector_lit[0]) {
        valid_sim_peak = 1;
        int gap_seen = 0;
        for (int i = 1; i < n_detectors; i++) {
            if (!detector_lit[i]) {
                gap_seen = 1;
            } else if (gap_seen) {
                valid_sim_peak = 0;
                break;
            }
        }
    }

    if (!valid_sim_peak) return;
    *out_peak_on_det = 1;

    /* Check contiguous overlap pattern */
    int valid_overlap = 0;
    int n_det_overlap = 0;

    if (spot_overlap[0]) {
        valid_overlap = 1;
        n_det_overlap = 1;
        int gap_seen = 0;
        for (int i = 1; i < n_detectors; i++) {
            if (spot_overlap[i]) {
                if (gap_seen && detector_lit[i]) {
                    valid_overlap = 0;
                    break;
                }
                n_det_overlap++;
            } else if (detector_lit[i]) {
                gap_seen = 1;
            }
        }
    }

    if (valid_overlap) {
        *out_peak_overlap = 1;
        *out_n_det_overlap = n_det_overlap;
    }
}

/* ---- Main triangle_overlap function (Python-facing, float32/float64) ---- */

static PyObject* triangle_overlap(PyObject *self, PyObject *args) {
    PyArrayObject *image_array;
    double v0x, v0y, v1x, v1y, v2x, v2y;

    if (!PyArg_ParseTuple(args, "O!dddddd",
            &PyArray_Type, &image_array,
            &v0x, &v0y, &v1x, &v1y, &v2x, &v2y))
        return NULL;

    /* Validate image array */
    if (PyArray_NDIM(image_array) != 2) {
        PyErr_SetString(PyExc_ValueError, "image must be 2D array");
        return NULL;
    }
    if (!PyArray_IS_C_CONTIGUOUS(image_array)) {
        PyErr_SetString(PyExc_ValueError, "image must be C-contiguous");
        return NULL;
    }

    npy_intp num_rows = PyArray_DIM(image_array, 0);
    npy_intp num_cols = PyArray_DIM(image_array, 1);
    int dtype = PyArray_TYPE(image_array);
    if (dtype != NPY_FLOAT32 && dtype != NPY_FLOAT64) {
        PyErr_SetString(PyExc_TypeError, "image must be float32 or float64");
        return NULL;
    }

    /* Truncate pixel coordinates (matching Python truncate_pixel) */
    #define TRUNC_PIXEL(val) ((val) < 0 ? -1.0 : floor(val))

    Point2D polygon[MAX_POLY_VERTS];
    polygon[0].x = TRUNC_PIXEL(v0x); polygon[0].y = TRUNC_PIXEL(v0y);
    polygon[1].x = TRUNC_PIXEL(v1x); polygon[1].y = TRUNC_PIXEL(v1y);
    polygon[2].x = TRUNC_PIXEL(v2x); polygon[2].y = TRUNC_PIXEL(v2y);
    #undef TRUNC_PIXEL

    /* Sutherland-Hodgman clip */
    int n_verts = sutherland_hodgman_clip(
        polygon, 3,
        0.0, (double)(num_cols - 1), 0.0, (double)(num_rows - 1)
    );

    if (n_verts < 3) {
        return Py_BuildValue("(ii)", 0, 0);
    }

    /* Round to integer vertices */
    int int_verts[MAX_POLY_VERTS][2];
    int y_min = INT_MAX, y_max = INT_MIN;
    for (int i = 0; i < n_verts; i++) {
        int_verts[i][0] = (int)round(polygon[i].x);
        int_verts[i][1] = (int)round(polygon[i].y);
        if (int_verts[i][1] < y_min) y_min = int_verts[i][1];
        if (int_verts[i][1] > y_max) y_max = int_verts[i][1];
    }

    int n_scanlines = y_max - y_min + 1;
    if (n_scanlines <= 0 || n_scanlines > MAX_SCANLINES) {
        return Py_BuildValue("(ii)", 0, 0);
    }

    /* Build edge table via Bresenham */
    ScanlineEntry *table = (ScanlineEntry*)malloc(n_scanlines * sizeof(ScanlineEntry));
    if (!table) {
        PyErr_NoMemory();
        return NULL;
    }
    for (int i = 0; i < n_scanlines; i++) {
        table[i].left = -1;
        table[i].right = -1;
    }

    for (int i = 0; i < n_verts; i++) {
        int prev = (i == 0) ? n_verts - 1 : i - 1;
        bresenham_edge(table, y_min, n_scanlines,
                       int_verts[prev][0], int_verts[prev][1],
                       int_verts[i][0], int_verts[i][1]);
    }

    /* Degenerate case: horizontal line */
    if (y_min == y_max) {
        int x_min_h = INT_MAX, x_max_h = INT_MIN;
        for (int i = 0; i < n_verts; i++) {
            if (int_verts[i][0] < x_min_h) x_min_h = int_verts[i][0];
            if (int_verts[i][0] > x_max_h) x_max_h = int_verts[i][0];
        }
        int num_lit = 0, num_overlap = 0;
        for (int x = x_min_h; x <= x_max_h; x++) {
            if (x >= 0 && x < num_cols && y_min >= 0 && y_min < num_rows) {
                num_lit++;
                /* Check pixel value */
                if (dtype == NPY_FLOAT32) {
                    float val = *(float*)PyArray_GETPTR2(image_array, y_min, x);
                    if (val > 0) num_overlap++;
                } else if (dtype == NPY_FLOAT64) {
                    double val = *(double*)PyArray_GETPTR2(image_array, y_min, x);
                    if (val > 0) num_overlap++;
                }
            }
        }
        free(table);
        return Py_BuildValue("(ii)", num_overlap, num_lit);
    }

    /* Fill scanlines and count overlap */
    int num_lit = 0, num_overlap = 0;

    for (int i = 0; i < n_scanlines; i++) {
        int left = table[i].left;
        int right = table[i].right;
        int y = y_min + i;

        if (left < 0 && right < 0) continue;
        if (left < 0) left = right;
        if (right < 0) right = left;
        if (left > right) { int tmp = left; left = right; right = tmp; }

        for (int x = left; x <= right; x++) {
            if (x >= 0 && x < num_cols && y >= 0 && y < num_rows) {
                num_lit++;
                if (dtype == NPY_FLOAT32) {
                    float val = *(float*)PyArray_GETPTR2(image_array, y, x);
                    if (val > 0) num_overlap++;
                } else if (dtype == NPY_FLOAT64) {
                    double val = *(double*)PyArray_GETPTR2(image_array, y, x);
                    if (val > 0) num_overlap++;
                }
            }
        }
    }

    free(table);
    return Py_BuildValue("(ii)", num_overlap, num_lit);
}

/* ---- pixel_radius_overlap function (Python-facing, float32/float64) ---- */

static PyObject* pixel_radius_overlap(PyObject *self, PyObject *args) {
    PyArrayObject *image_array;
    int cx, cy, radius;

    if (!PyArg_ParseTuple(args, "O!iii",
            &PyArray_Type, &image_array,
            &cx, &cy, &radius))
        return NULL;

    if (PyArray_NDIM(image_array) != 2) {
        PyErr_SetString(PyExc_ValueError, "image must be 2D array");
        return NULL;
    }

    npy_intp num_rows = PyArray_DIM(image_array, 0);
    npy_intp num_cols = PyArray_DIM(image_array, 1);
    int dtype = PyArray_TYPE(image_array);
    if (dtype != NPY_FLOAT32 && dtype != NPY_FLOAT64) {
        PyErr_SetString(PyExc_TypeError, "image must be float32 or float64");
        return NULL;
    }

    int found_in_bounds = 0;
    int found_bright = 0;

    for (int dy = -radius; dy <= radius && !found_bright; dy++) {
        for (int dx = -radius; dx <= radius; dx++) {
            int px = cx + dx;
            int py = cy + dy;
            if (px >= 0 && px < num_cols && py >= 0 && py < num_rows) {
                found_in_bounds = 1;
                if (dtype == NPY_FLOAT32) {
                    float val = *(float*)PyArray_GETPTR2(image_array, py, px);
                    if (val > 0) { found_bright = 1; break; }
                } else if (dtype == NPY_FLOAT64) {
                    double val = *(double*)PyArray_GETPTR2(image_array, py, px);
                    if (val > 0) { found_bright = 1; break; }
                }
            }
        }
    }

    return Py_BuildValue("(ii)", found_in_bounds, found_bright);
}

/* ---- Stage D: one peak (shared by stage_d_overlap and the fused evaluate, _fused_eval.c) ---- */

/**
 * icenine_stage_d_peak — the per-peak body of stage D: overlap of one peak on every detector, the
 * optional ray-consistency rule, the qualified-peak count and the accumulation into ``acc``
 * (OverlapInfo.update_counts and the Welford mean of OverlapInfo.update_quality).
 *
 *   images[d]: the binary frame (uint8, image_rows x image_cols) of this peak on detector d, or
 *              NULL (no image: the detector is skipped, as an index outside the image list).
 *   det_hit[d]: 1 if all three vertices hit detector d's plane (t > 0).
 *   centers: (n_det, 2) int [col, row] when pixel_radius > 0; vertices: (n_det, 3, 2) float32
 *            [col, row] when pixel_radius == 0 (converted to double, as stage_d_overlap always did).
 *   rule_bits: 1 = C++ CountQualifiedPeaks rule, 2 = CppEdgeRules detector-lit rule.
 *   obs_in, obs_bright: (n_det,) uint8 of the ray rule, or both NULL; cases: its 8 counters.
 *
 * Not static: _fused_eval.c (a second translation unit of this module) calls it, so stage D exists
 * once.
 */
void icenine_stage_d_peak(
    const unsigned char *const *images, npy_intp image_rows, npy_intp image_cols,
    const npy_uint8 *det_hit, const npy_int32 *centers, const float *vertices,
    int n_detectors, int pixel_radius, int rule_bits,
    const npy_uint8 *obs_in, const npy_uint8 *obs_bright, long *cases,
    IceStageDAcc *acc
) {
    /* rule_bits is a bit field: 1 = C++ CountQualifiedPeaks rule, 2 = CppEdgeRules detector-lit
       rule, exactly C++ ShapePixelOverlapCounter (OverlapInfo.h:501-514): lit iff the pixel overlap
       is > 0 or a vertex is in bounds, so a triangle with pixels but no vertex in bounds stays unlit
       (without bit 2: lit on any produced pixel). Bit 2 acts only in the triangle (pixel_radius <= 0)
       branch; 1 is what older callers pass */
    const int cpp_vertex_lit = (rule_bits & 2) != 0;

    int detector_lit[MAX_DETECTORS];
    int spot_overlap[MAX_DETECTORS];
    int peak_pixel_overlap = 0;
    int peak_pixel_on_det = 0;

    memset(detector_lit, 0, n_detectors * sizeof(int));
    memset(spot_overlap, 0, n_detectors * sizeof(int));

    for (int d = 0; d < n_detectors; d++) {
        if (!det_hit[d]) continue;
        const unsigned char *image_data = images[d];
        if (image_data == NULL) continue;

        if (pixel_radius > 0) {
            int cx = centers[2 * d];
            int cy = centers[2 * d + 1];

            int in_bounds, bright;
            pixel_radius_overlap_uint8(
                image_data, image_rows, image_cols,
                cx, cy, pixel_radius,
                &in_bounds, &bright
            );

            if (in_bounds) {
                detector_lit[d] = 1;
                peak_pixel_on_det++;
            }
            if (bright) {
                peak_pixel_overlap++;
                spot_overlap[d] = 1;
            }
        } else {
            const float *v = vertices + 6 * d;
            double v0x = (double)v[0];
            double v0y = (double)v[1];
            double v1x = (double)v[2];
            double v1y = (double)v[3];
            double v2x = (double)v[4];
            double v2y = (double)v[5];

            int n_overlap, n_lit;
            triangle_overlap_uint8(
                image_data, image_rows, image_cols,
                v0x, v0y, v1x, v1y, v2x, v2y,
                &n_overlap, &n_lit
            );

            peak_pixel_overlap += n_overlap;
            peak_pixel_on_det += n_lit;

            if (n_overlap > 0) {
                detector_lit[d] = 1;
                spot_overlap[d] = 1;
            } else if (cpp_vertex_lit) {
                /* CppEdgeRules, exact C++ ShapePixelOverlapCounter (OverlapInfo.h:501-514): a
                   detector is lit iff the pixel overlap is > 0 (branch above) or one of the
                   three vertices is in bounds (IsInBound). There is no "pixels produced"
                   clause, so a triangle that covers image pixels with all three vertices out of
                   bounds (near a corner) stays unlit. Without CppEdgeRules the default rule
                   (branch below: lit on any produced pixel) is unchanged. */
                double vx[3] = {v0x, v1x, v2x}, vy[3] = {v0y, v1y, v2y};
                for (int q = 0; q < 3; q++) {
                    double px = (vx[q] < 0) ? -1.0 : floor(vx[q]);
                    double py = (vy[q] < 0) ? -1.0 : floor(vy[q]);
                    if (px >= 0 && px < image_cols && py >= 0 && py < image_rows) {
                        detector_lit[d] = 1;
                        break;
                    }
                }
            } else if (n_lit > 0) {
                detector_lit[d] = 1;
            }
        }
    }

    /* Ray-consistency rule (see ray_overlap.py): once the predicted spot is bright on
       detector 0, a later detector d whose predicted spot is not bright is decided by the
       observed ray. a: predicted bright (nothing). b: predicted lit and dark, observed ray
       bright -> overlap. c: both in bounds and dark -> unchanged (invalid under the C++
       rule). d: both out of bounds -> unchanged (unlit tolerated). e: predicted lit and
       dark, observed out of bounds -> treated as unlit. f: predicted out of bounds,
       observed in bounds and dark -> unchanged (always forgiven, as in the C++ rule; the
       invalidating variant f1 was dropped 2026-10-12). Predicted out of bounds with the observed ray bright is
       counted separately (f_bright) and left unchanged.
       Pixel accounting in triangle mode (pixel_radius 0): the detector-d triangle already
       added n_lit pixels to peak_pixel_on_det. Case b adds ONE overlap pixel (a partial
       credit against n_lit pixels) and case e removes ONE on-detector pixel (an "unlit"
       detector keeps n_lit-1 pixels in the denominator). Owner decision 2026-10-12: kept as
       is, to be revisited; see open question 9 in
       docs/findings_finisher_mc_study_2026-10.md. */
    if (obs_in && spot_overlap[0]) {
        cases[7]++;
        for (int d = 1; d < n_detectors; d++) {
            int oi = obs_in[d];
            int ob = obs_bright[d];
            if (spot_overlap[d]) { cases[0]++; continue; }
            if (detector_lit[d]) {
                if (ob) { spot_overlap[d] = 1; peak_pixel_overlap++; cases[1]++; }
                else if (!oi) { detector_lit[d] = 0; peak_pixel_on_det--; cases[4]++; }
                else cases[2]++;
            } else if (!oi) {
                cases[3]++;
            } else if (ob) {
                cases[6]++;
            } else {
                cases[5]++;
            }
        }
    }

    /* Qualified peak counting */
    int peak_on_det, peak_ovlp, n_det_ovlp;
    count_qualified_peaks_c(
        detector_lit, spot_overlap, n_detectors,
        &peak_on_det, &peak_ovlp, &n_det_ovlp, rule_bits & 1
    );

    /* Update counts (matches OverlapInfo.update_counts) */
    acc->last_n_det_ovlp = n_det_ovlp;
    acc->pixel_overlap += peak_pixel_overlap;
    acc->pixel_on_det += peak_pixel_on_det;
    acc->peak_overlap += peak_ovlp;
    acc->peak_on_det += peak_on_det;

    /* Update quality (Welford mean, matches OverlapInfo.update_quality) */
    if (peak_pixel_on_det > 0 && n_detectors > 0) {
        double pixel_ratio = (double)peak_pixel_overlap / peak_pixel_on_det;
        double cur_quality;
        if (n_det_ovlp > 0) {
            double det_ratio = (double)n_det_ovlp / n_detectors;
            cur_quality = pixel_ratio * det_ratio;
        } else {
            cur_quality = 0.0;
        }
        acc->quality += (cur_quality - acc->quality) / (acc->n_quality_points + 1);
        acc->n_quality_points++;
    }
}

/* ---- Batch Stage D overlap function ---- */

/**
 * stage_d_overlap(image_list, wedge_indices, det_hit_mask,
 *                 pixel_centers, pixel_vertices,
 *                 n_detectors, pixel_radius, n_peaks,
 *                 image_rows, image_cols)
 *
 * Replaces the entire Stage D Python loop with a single C call.
 * Processes all M peaks × N detectors, computing overlap against
 * pre-cached uint8 binary images.
 *
 * Args:
 *   image_list: Python list of uint8 2D numpy arrays (binary images).
 *               Indexed as image_list[wedge_local * n_detectors + det].
 *   wedge_indices: int32 array (M,) — local wedge index for each peak
 *   det_hit_mask: uint8 array (M, N_det) — 1 if peak hits detector
 *   pixel_centers: int32 array (M, N_det, 2) — (cx, cy) for pixel_radius mode
 *                  (may be None if pixel_radius == 0)
 *   pixel_vertices: float32 array (M, N_det, 3, 2) — triangle vertices
 *                   (may be None if pixel_radius > 0)
 *   n_detectors: int
 *   pixel_radius: int (0 for triangle mode)
 *   n_peaks: int (M)
 *   image_rows: int (rows per image)
 *   image_cols: int (cols per image)
 *
 * Returns:
 *   (pixel_overlap, pixel_on_detector, peak_overlap, peak_on_detector,
 *    quality, n_quality_points)
 */
static PyObject* stage_d_overlap(PyObject *self, PyObject *args) {
    PyObject *image_list;
    PyArrayObject *wedge_indices_arr, *det_hit_mask_arr;
    PyObject *pixel_centers_obj, *pixel_vertices_obj;
    int n_detectors, pixel_radius, n_peaks, image_rows, image_cols;
    int cpp_rule = 0;
    PyObject *obs_in_obj = Py_None, *obs_bright_obj = Py_None;

    if (!PyArg_ParseTuple(args, "O!O!O!OOiiiii|iOO",
            &PyList_Type, &image_list,
            &PyArray_Type, &wedge_indices_arr,
            &PyArray_Type, &det_hit_mask_arr,
            &pixel_centers_obj,
            &pixel_vertices_obj,
            &n_detectors, &pixel_radius, &n_peaks,
            &image_rows, &image_cols, &cpp_rule,
            &obs_in_obj, &obs_bright_obj))
        return NULL;

    /* Optional ray-consistency inputs (icenine/ray_overlap.py): (M, N_det) uint8, observed ray
       in bounds / bright on detector d. Both given -> the combined rule is applied. */
    const npy_uint8 *obs_in = NULL, *obs_bright = NULL;
    if (obs_in_obj != Py_None || obs_bright_obj != Py_None) {
        if (obs_in_obj == Py_None || obs_bright_obj == Py_None ||
            !PyArray_Check(obs_in_obj) || !PyArray_Check(obs_bright_obj) ||
            PyArray_TYPE((PyArrayObject *)obs_in_obj) != NPY_UINT8 ||
            PyArray_TYPE((PyArrayObject *)obs_bright_obj) != NPY_UINT8 ||
            !PyArray_IS_C_CONTIGUOUS((PyArrayObject *)obs_in_obj) ||
            !PyArray_IS_C_CONTIGUOUS((PyArrayObject *)obs_bright_obj) ||
            PyArray_SIZE((PyArrayObject *)obs_in_obj) < (npy_intp)n_peaks * n_detectors ||
            PyArray_SIZE((PyArrayObject *)obs_bright_obj) < (npy_intp)n_peaks * n_detectors) {
            PyErr_SetString(PyExc_TypeError,
                "obs_in/obs_bright must both be C-contiguous (M, N_det) uint8 arrays");
            return NULL;
        }
        obs_in = (const npy_uint8 *)PyArray_DATA((PyArrayObject *)obs_in_obj);
        obs_bright = (const npy_uint8 *)PyArray_DATA((PyArrayObject *)obs_bright_obj);
    }
    /* case counters: a b c d e f_dark f_bright det0_bright */
    long cases[8] = {0, 0, 0, 0, 0, 0, 0, 0};

    if (n_detectors > MAX_DETECTORS) {
        PyErr_SetString(PyExc_ValueError, "n_detectors exceeds MAX_DETECTORS");
        return NULL;
    }

    /* Validate wedge_indices */
    if (PyArray_NDIM(wedge_indices_arr) != 1 ||
        PyArray_TYPE(wedge_indices_arr) != NPY_INT32 ||
        !PyArray_IS_C_CONTIGUOUS(wedge_indices_arr)) {
        PyErr_SetString(PyExc_TypeError,
            "wedge_indices must be C-contiguous 1D int32 array");
        return NULL;
    }

    /* Validate det_hit_mask */
    if (PyArray_NDIM(det_hit_mask_arr) != 2 ||
        PyArray_TYPE(det_hit_mask_arr) != NPY_UINT8 ||
        !PyArray_IS_C_CONTIGUOUS(det_hit_mask_arr)) {
        PyErr_SetString(PyExc_TypeError,
            "det_hit_mask must be C-contiguous 2D uint8 array");
        return NULL;
    }

    const npy_int32 *wedge_indices = (const npy_int32 *)PyArray_DATA(wedge_indices_arr);
    const npy_uint8 *det_hit_mask = (const npy_uint8 *)PyArray_DATA(det_hit_mask_arr);

    /* Optional arrays */
    const npy_int32 *pixel_centers = NULL;
    const float *pixel_vertices = NULL;

    PyArrayObject *pixel_centers_arr = NULL;
    PyArrayObject *pixel_vertices_arr = NULL;

    if (pixel_radius > 0) {
        if (pixel_centers_obj == Py_None) {
            PyErr_SetString(PyExc_ValueError,
                "pixel_centers required when pixel_radius > 0");
            return NULL;
        }
        pixel_centers_arr = (PyArrayObject *)pixel_centers_obj;
        if (PyArray_TYPE(pixel_centers_arr) != NPY_INT32 ||
            !PyArray_IS_C_CONTIGUOUS(pixel_centers_arr)) {
            PyErr_SetString(PyExc_TypeError,
                "pixel_centers must be C-contiguous int32 array");
            return NULL;
        }
        pixel_centers = (const npy_int32 *)PyArray_DATA(pixel_centers_arr);
    } else {
        if (pixel_vertices_obj == Py_None) {
            PyErr_SetString(PyExc_ValueError,
                "pixel_vertices required when pixel_radius == 0");
            return NULL;
        }
        pixel_vertices_arr = (PyArrayObject *)pixel_vertices_obj;
        if (PyArray_TYPE(pixel_vertices_arr) != NPY_FLOAT32 ||
            !PyArray_IS_C_CONTIGUOUS(pixel_vertices_arr)) {
            PyErr_SetString(PyExc_TypeError,
                "pixel_vertices must be C-contiguous float32 array");
            return NULL;
        }
        pixel_vertices = (const float *)PyArray_DATA(pixel_vertices_arr);
    }

    /* Pre-fetch image data pointers */
    Py_ssize_t n_images = PyList_Size(image_list);
    const unsigned char **image_ptrs = (const unsigned char **)malloc(
        n_images * sizeof(unsigned char *));
    if (!image_ptrs) {
        PyErr_NoMemory();
        return NULL;
    }

    for (Py_ssize_t i = 0; i < n_images; i++) {
        PyObject *item = PyList_GET_ITEM(image_list, i);
        if (!PyArray_Check(item)) {
            free(image_ptrs);
            PyErr_SetString(PyExc_TypeError, "all images must be numpy arrays");
            return NULL;
        }
        PyArrayObject *img = (PyArrayObject *)item;
        if (PyArray_TYPE(img) != NPY_UINT8 || PyArray_NDIM(img) != 2 ||
            !PyArray_IS_C_CONTIGUOUS(img)) {
            free(image_ptrs);
            PyErr_SetString(PyExc_TypeError,
                "all images must be C-contiguous 2D uint8 arrays");
            return NULL;
        }
        image_ptrs[i] = (const unsigned char *)PyArray_DATA(img);
    }

    /* Main loop over peaks: the per-peak body is icenine_stage_d_peak (shared with the fused
       evaluate, _fused_eval.c) */
    IceStageDAcc acc = {0, 0, 0, 0, 0, 0, 0.0};
    const unsigned char *peak_images[MAX_DETECTORS];
    for (int p = 0; p < n_peaks; p++) {
        int wedge_idx = wedge_indices[p];
        /* image_list[wedge_idx * n_detectors + d]; NULL outside the list (skipped) */
        for (int d = 0; d < n_detectors; d++) {
            Py_ssize_t img_idx = (Py_ssize_t)wedge_idx * n_detectors + d;
            peak_images[d] = (img_idx < 0 || img_idx >= n_images) ? NULL : image_ptrs[img_idx];
        }
        icenine_stage_d_peak(
            peak_images, image_rows, image_cols,
            det_hit_mask + (npy_intp)p * n_detectors,
            pixel_centers ? pixel_centers + (npy_intp)p * n_detectors * 2 : NULL,
            pixel_vertices ? pixel_vertices + (npy_intp)p * n_detectors * 6 : NULL,
            n_detectors, pixel_radius, cpp_rule,
            obs_in ? obs_in + (npy_intp)p * n_detectors : NULL,
            obs_bright ? obs_bright + (npy_intp)p * n_detectors : NULL,
            cases, &acc
        );
    }

    free(image_ptrs);

    if (obs_in) {
        return Py_BuildValue("(iiiidii(llllllll))",
            acc.pixel_overlap, acc.pixel_on_det, acc.peak_overlap, acc.peak_on_det,
            acc.quality, acc.n_quality_points, acc.last_n_det_ovlp,
            cases[0], cases[1], cases[2], cases[3], cases[4], cases[5], cases[6], cases[7]);
    }
    return Py_BuildValue("(iiiidii)",
        acc.pixel_overlap,
        acc.pixel_on_det,
        acc.peak_overlap,
        acc.peak_on_det,
        acc.quality,
        acc.n_quality_points,
        acc.last_n_det_ovlp  /* detectors_overlap from last peak, matching Python fallback */
    );
}

/* ---- Observed-ray flags (C port of ray_overlap.ray_flags, which stays the reference) ---- */

#define RAY_MAX_RAYS 64
#define RAY_GEOM_LEN 20 /* n[3] d pos[3] org[3] bj[3] bk[3] pw ph hw hh */

#pragma clang fp contract(off)
#pragma STDC FP_CONTRACT OFF

static inline double dot3(const double *a, const double *b) {
    return (a[0] * b[0] + a[1] * b[1]) + a[2] * b[2];
}

static inline long to_long(double x) {
    if (!(x > -1e15)) return -1000000000000000L; /* also NaN */
    if (x > 1e15) return 1000000000000000L;
    return (long)x;
}

/**
 * ray_flags(image_list, local_wedge, det_hit0, det0_px, origins, geom, n_det,
 *           pixel_radius, arm, radius_px, voxel_size, omega_step, registration)
 *   -> (obs_in, obs_bright), both (M, n_det) uint8.
 * image_list is indexed [local_wedge * n_det + d] like stage_d_overlap; det0_px is (M, 3, 2)
 * float64 [col, row]; origins (M, 3) float64; geom (n_det, 20) float64. arm 0 = R, 1 = A.
 * The arithmetic follows the numpy reference operation by operation (no FMA contraction).
 */
static PyObject* ray_flags_c(PyObject *self, PyObject *args) {
    PyObject *image_list;
    PyArrayObject *wedge_arr, *hit_arr, *px_arr, *org_arr, *geom_arr;
    int n_det, pixel_radius, arm, radius_px;
    double voxel_size, omega_step, registration;
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!iiiiddd",
            &PyList_Type, &image_list,
            &PyArray_Type, &wedge_arr, &PyArray_Type, &hit_arr, &PyArray_Type, &px_arr,
            &PyArray_Type, &org_arr, &PyArray_Type, &geom_arr,
            &n_det, &pixel_radius, &arm, &radius_px, &voxel_size, &omega_step, &registration))
        return NULL;
    if (n_det > MAX_DETECTORS || n_det < 1) {
        PyErr_SetString(PyExc_ValueError, "bad n_det");
        return NULL;
    }
    npy_intp m = PyArray_DIM(wedge_arr, 0);
    if (PyArray_TYPE(wedge_arr) != NPY_INT32 || !PyArray_IS_C_CONTIGUOUS(wedge_arr) ||
        PyArray_TYPE(hit_arr) != NPY_UINT8 || !PyArray_IS_C_CONTIGUOUS(hit_arr) ||
        PyArray_SIZE(hit_arr) < m ||
        PyArray_TYPE(px_arr) != NPY_FLOAT64 || !PyArray_IS_C_CONTIGUOUS(px_arr) ||
        PyArray_SIZE(px_arr) < m * 6 ||
        PyArray_TYPE(org_arr) != NPY_FLOAT64 || !PyArray_IS_C_CONTIGUOUS(org_arr) ||
        PyArray_SIZE(org_arr) < m * 3 ||
        PyArray_TYPE(geom_arr) != NPY_FLOAT64 || !PyArray_IS_C_CONTIGUOUS(geom_arr) ||
        PyArray_SIZE(geom_arr) < (npy_intp)n_det * RAY_GEOM_LEN) {
        PyErr_SetString(PyExc_TypeError, "ray_flags: bad array dtype/shape/contiguity");
        return NULL;
    }
    const int *wedge = (const int *)PyArray_DATA(wedge_arr);
    const unsigned char *hit0 = (const unsigned char *)PyArray_DATA(hit_arr);
    const double *px = (const double *)PyArray_DATA(px_arr);
    const double *orgs = (const double *)PyArray_DATA(org_arr);
    const double *geom = (const double *)PyArray_DATA(geom_arr);

    npy_intp dims[2] = {m, n_det};
    PyArrayObject *oin = (PyArrayObject *)PyArray_ZEROS(2, dims, NPY_UINT8, 0);
    PyArrayObject *obr = (PyArrayObject *)PyArray_ZEROS(2, dims, NPY_UINT8, 0);
    if (!oin || !obr) { Py_XDECREF(oin); Py_XDECREF(obr); return NULL; }
    unsigned char *out_in = (unsigned char *)PyArray_DATA(oin);
    unsigned char *out_br = (unsigned char *)PyArray_DATA(obr);

    Py_ssize_t n_images = PyList_Size(image_list);
    const unsigned char *ip[MAX_DETECTORS];
    long irows[MAX_DETECTORS], icols[MAX_DETECTORS];

    if (n_det >= 2) {
        for (npy_intp p = 0; p < m; p++) {
            if (!hit0[p]) continue;
            long base = (long)wedge[p] * n_det;
            if (base < 0 || base + n_det > n_images) {
                Py_DECREF(oin); Py_DECREF(obr);
                PyErr_SetString(PyExc_IndexError, "ray_flags: wedge index outside image_list");
                return NULL;
            }
            for (int d = 0; d < n_det; d++) {
                PyObject *item = PyList_GET_ITEM(image_list, base + d);
                if (!PyArray_Check(item) || PyArray_TYPE((PyArrayObject *)item) != NPY_UINT8 ||
                    PyArray_NDIM((PyArrayObject *)item) != 2 ||
                    !PyArray_IS_C_CONTIGUOUS((PyArrayObject *)item)) {
                    Py_DECREF(oin); Py_DECREF(obr);
                    PyErr_SetString(PyExc_TypeError, "images must be C-contiguous 2D uint8");
                    return NULL;
                }
                ip[d] = (const unsigned char *)PyArray_DATA((PyArrayObject *)item);
                irows[d] = (long)PyArray_DIM((PyArrayObject *)item, 0);
                icols[d] = (long)PyArray_DIM((PyArrayObject *)item, 1);
            }
            const double *v = px + p * 6; /* (3, 2): [col, row] per vertex */
            long c_lo, c_hi, r_lo, r_hi;
            if (pixel_radius > 0) {
                long cx = to_long(trunc(v[0])), cy = to_long(trunc(v[1]));
                c_lo = cx - pixel_radius; c_hi = cx + pixel_radius;
                r_lo = cy - pixel_radius; r_hi = cy + pixel_radius;
            } else {
                double cmin = v[0], cmax = v[0], rmin = v[1], rmax = v[1];
                for (int k = 1; k < 3; k++) {
                    if (v[2 * k] < cmin) cmin = v[2 * k];
                    if (v[2 * k] > cmax) cmax = v[2 * k];
                    if (v[2 * k + 1] < rmin) rmin = v[2 * k + 1];
                    if (v[2 * k + 1] > rmax) rmax = v[2 * k + 1];
                }
                c_lo = to_long(floor(cmin)); c_hi = to_long(floor(cmax));
                r_lo = to_long(floor(rmin)); r_hi = to_long(floor(rmax));
            }
            if (c_lo < 0) c_lo = 0;
            if (c_hi > icols[0] - 1) c_hi = icols[0] - 1;
            if (r_lo < 0) r_lo = 0;
            if (r_hi > irows[0] - 1) r_hi = irows[0] - 1;
            if (c_lo > c_hi || r_lo > r_hi) continue;

            /* bright detector-0 pixels of the window, row-major, at most RAY_MAX_RAYS */
            int ray_c[RAY_MAX_RAYS], ray_r[RAY_MAX_RAYS];
            int nr = 0;
            for (long r = r_lo; r <= r_hi && nr < RAY_MAX_RAYS; r++)
                for (long c = c_lo; c <= c_hi && nr < RAY_MAX_RAYS; c++)
                    if (ip[0][r * icols[0] + c]) { ray_c[nr] = (int)c; ray_r[nr] = (int)r; nr++; }
            if (nr == 0) continue;

            const double *g0 = geom;
            const double *origin = orgs + p * 3;
            double q_base[3];
            for (int i = 0; i < 3; i++) q_base[i] = g0[4 + i] + g0[7 + i];
            double rad = hypot(origin[0], origin[1]);
            double arc[3] = {-origin[1], origin[0], 0.0};
            double nrm = sqrt(dot3(arc, arc));
            if (nrm > 0) { arc[0] /= nrm; arc[1] /= nrm; arc[2] /= nrm; }

            double q0[RAY_MAX_RAYS][3];
            for (int i = 0; i < nr; i++) {
                double j = ((double)ray_c[i] + 0.5) * g0[16] - g0[18];
                double k = ((double)ray_r[i] + 0.5) * g0[17] - g0[19];
                for (int a = 0; a < 3; a++)
                    q0[i][a] = (q_base[a] + j * g0[10 + a]) + k * g0[13 + a];
            }
            for (int d = 1; d < n_det; d++) {
                const double *g = geom + d * RAY_GEOM_LEN;
                double on = dot3(origin, g); /* origin . n */
                double num = -(on + g[3]);
                int any_in = 0, any_bright = 0;
                for (int i = 0; i < nr && !any_bright; i++) {
                    double dirs[3];
                    for (int a = 0; a < 3; a++) dirs[a] = q0[i][a] - origin[a];
                    double denom = dot3(dirs, g);
                    if (!(fabs(denom) > 1e-12)) continue;
                    double t = num / denom;
                    if (!(t > 0)) continue;
                    double rel[3];
                    for (int a = 0; a < 3; a++) {
                        double hit = origin[a] + t * dirs[a];
                        rel[a] = (hit - g[4 + a]) - g[7 + a];
                    }
                    double col = (dot3(rel, g + 10) + g[18]) / g[16];
                    double row = (dot3(rel, g + 13) + g[19]) / g[17];
                    long lo_c, hi_c, lo_r, hi_r;
                    if (arm == 0) {
                        lo_c = to_long(trunc(col)) - radius_px; hi_c = lo_c + 2 * radius_px;
                        lo_r = to_long(trunc(row)) - radius_px; hi_r = lo_r + 2 * radius_px;
                    } else {
                        /* Arm A window. WINDOW DECISION (owner, 2026-10-12): every term is a
                         * half-width (intended reading); the window stays as set, untuned.
                         * Same formula and note as the ray_overlap.py module docstring. */
                        double ratio = t;
                        double iso = (voxel_size * (ratio - 1.0)
                                      + 0.5 * (g0[16] + g0[17]) * ratio) + registration;
                        double arc_mm = rad * omega_step * (ratio - 1.0);
                        double hc = (iso + arc_mm * fabs(dot3(arc, g + 10))) / g[16];
                        double hr = (iso + arc_mm * fabs(dot3(arc, g + 13))) / g[17];
                        lo_c = to_long(floor(col - hc)); hi_c = to_long(floor(col + hc));
                        lo_r = to_long(floor(row - hr)); hi_r = to_long(floor(row + hr));
                    }
                    /* _box_flags */
                    if (hi_c < 0 || hi_r < 0 || lo_c >= icols[d] || lo_r >= irows[d]) continue;
                    any_in = 1;
                    long a0 = lo_c < 0 ? 0 : lo_c, a1 = hi_c > icols[d] - 1 ? icols[d] - 1 : hi_c;
                    long b0 = lo_r < 0 ? 0 : lo_r, b1 = hi_r > irows[d] - 1 ? irows[d] - 1 : hi_r;
                    for (long r = b0; r <= b1 && !any_bright; r++)
                        for (long c = a0; c <= a1; c++)
                            if (ip[d][r * icols[d] + c]) { any_bright = 1; break; }
                }
                out_in[p * n_det + d] = (unsigned char)any_in;
                out_br[p * n_det + d] = (unsigned char)any_bright;
            }
        }
    }
    return Py_BuildValue("(NN)", (PyObject *)oin, (PyObject *)obr);
}

/* ---- Module definition ---- */

#ifdef ICENINE_FUSED
/* the fused evaluate (S2), defined in _fused_eval.c */
PyObject *icenine_fused_context(PyObject *self, PyObject *args);
PyObject *icenine_fused_evaluate(PyObject *self, PyObject *args);
#endif

static PyMethodDef rasterize_methods[] = {
    {"triangle_overlap", triangle_overlap, METH_VARARGS,
     "triangle_overlap(image, v0x, v0y, v1x, v1y, v2x, v2y) -> (overlap, lit)\n\n"
     "Fast triangle rasterization + overlap counting against experimental image.\n"
     "Implements Sutherland-Hodgman clipping + Bresenham scanline fill."},
    {"pixel_radius_overlap", pixel_radius_overlap, METH_VARARGS,
     "pixel_radius_overlap(image, cx, cy, radius) -> (in_bounds, bright)\n\n"
     "Check if any pixel in +-radius square around (cx, cy) is bright."},
    {"ray_flags", ray_flags_c, METH_VARARGS,
     "ray_flags(image_list, local_wedge, det_hit0, det0_px, origins, geom, n_det, pixel_radius,\n"
     "          arm, radius_px, voxel_size, omega_step, registration) -> (obs_in, obs_bright)"},
    {"stage_d_overlap", stage_d_overlap, METH_VARARGS,
     "stage_d_overlap(image_list, wedge_indices, det_hit_mask, pixel_centers,\n"
     "                pixel_vertices, n_detectors, pixel_radius, n_peaks,\n"
     "                image_rows, image_cols)\n"
     "  -> (pixel_overlap, pixel_on_det, peak_overlap, peak_on_det,\n"
     "      quality, n_quality_points, detectors_overlap)\n\n"
     "Batch Stage D overlap: processes all peaks in a single C call.\n"
     "Uses uint8 binary images from ImageData.get_binary_numpy()."},
#ifdef ICENINE_FUSED
    {"fused_context", icenine_fused_context, METH_VARARGS,
     "fused_context(...) -> capsule: the per-cost-function context of fused_evaluate\n"
     "(see _fused_eval.c and cost_functions._FusedContext)."},
    {"fused_evaluate", icenine_fused_evaluate, METH_VARARGS,
     "fused_evaluate(ctx, R, verts[, debug]) -> (pixel_overlap, pixel_on_det, peak_overlap,\n"
     "    peak_on_det, quality, n_quality_points, detectors_overlap), or the list of frames the\n"
     "    FrameTable must fill first; debug=1: a dict of the intermediates (_fused_eval.c)."},
#endif
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef rasterize_module = {
    PyModuleDef_HEAD_INIT,
    "_rasterize",
    "Fast triangle rasterization C extension for IceNine cost functions.",
    -1,
    rasterize_methods
};

PyMODINIT_FUNC PyInit__rasterize(void) {
    import_array();
    PyObject *m = PyModule_Create(&rasterize_module);
    if (m != NULL) PyModule_AddIntConstant(m, "HAS_CPP_RULE", 1);
    if (m != NULL) PyModule_AddIntConstant(m, "HAS_RAY_RULE", 1);
    /* 1 = vertex-lit rule that also lit on produced pixels (superset of C++); 2 = exact C++ rule */
    if (m != NULL) PyModule_AddIntConstant(m, "HAS_VERTEX_LIT", 2);
#ifdef ICENINE_FUSED
    /* 1 = fused_context / fused_evaluate (S2, the FusedCost key; _fused_eval.c) */
    if (m != NULL) PyModule_AddIntConstant(m, "HAS_FUSED", 1);
#endif
    return m;
}
