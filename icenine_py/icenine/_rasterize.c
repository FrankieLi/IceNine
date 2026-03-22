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
    int *out_peak_on_det, int *out_peak_overlap, int *out_n_det_overlap
) {
    *out_peak_on_det = 0;
    *out_peak_overlap = 0;
    *out_n_det_overlap = 0;

    if (n_detectors == 0) return;

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

    if (!PyArg_ParseTuple(args, "O!O!O!OOiiiii",
            &PyList_Type, &image_list,
            &PyArray_Type, &wedge_indices_arr,
            &PyArray_Type, &det_hit_mask_arr,
            &pixel_centers_obj,
            &pixel_vertices_obj,
            &n_detectors, &pixel_radius, &n_peaks,
            &image_rows, &image_cols))
        return NULL;

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

    /* Accumulation variables */
    int total_pixel_overlap = 0;
    int total_pixel_on_det = 0;
    int total_peak_overlap = 0;
    int total_peak_on_det = 0;
    double quality = 0.0;
    int n_quality_points = 0;
    int last_n_det_ovlp = 0;  /* detectors_overlap from last peak */

    /* Main loop over peaks */
    for (int p = 0; p < n_peaks; p++) {
        int wedge_idx = wedge_indices[p];

        int detector_lit[MAX_DETECTORS];
        int spot_overlap[MAX_DETECTORS];
        int peak_pixel_overlap = 0;
        int peak_pixel_on_det = 0;

        memset(detector_lit, 0, n_detectors * sizeof(int));
        memset(spot_overlap, 0, n_detectors * sizeof(int));

        for (int d = 0; d < n_detectors; d++) {
            /* det_hit_mask is (M, N_det), row-major */
            if (!det_hit_mask[p * n_detectors + d]) continue;

            /* Look up image: image_list[wedge_idx * n_detectors + d] */
            Py_ssize_t img_idx = (Py_ssize_t)wedge_idx * n_detectors + d;
            if (img_idx < 0 || img_idx >= n_images) continue;

            const unsigned char *image_data = image_ptrs[img_idx];

            if (pixel_radius > 0) {
                /* pixel_centers is (M, N_det, 2), row-major */
                int base = (p * n_detectors + d) * 2;
                int cx = pixel_centers[base];
                int cy = pixel_centers[base + 1];

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
                /* pixel_vertices is (M, N_det, 3, 2), row-major */
                int base = ((p * n_detectors + d) * 3) * 2;
                double v0x = (double)pixel_vertices[base + 0];
                double v0y = (double)pixel_vertices[base + 1];
                double v1x = (double)pixel_vertices[base + 2];
                double v1y = (double)pixel_vertices[base + 3];
                double v2x = (double)pixel_vertices[base + 4];
                double v2y = (double)pixel_vertices[base + 5];

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
                } else if (n_lit > 0) {
                    detector_lit[d] = 1;
                }
            }
        }

        /* Qualified peak counting */
        int peak_on_det, peak_ovlp, n_det_ovlp;
        count_qualified_peaks_c(
            detector_lit, spot_overlap, n_detectors,
            &peak_on_det, &peak_ovlp, &n_det_ovlp
        );

        /* Update counts (matches OverlapInfo.update_counts) */
        last_n_det_ovlp = n_det_ovlp;
        total_pixel_overlap += peak_pixel_overlap;
        total_pixel_on_det += peak_pixel_on_det;
        total_peak_overlap += peak_ovlp;
        total_peak_on_det += peak_on_det;

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
            quality += (cur_quality - quality) / (n_quality_points + 1);
            n_quality_points++;
        }
    }

    free(image_ptrs);

    return Py_BuildValue("(iiiidii)",
        total_pixel_overlap,
        total_pixel_on_det,
        total_peak_overlap,
        total_peak_on_det,
        quality,
        n_quality_points,
        last_n_det_ovlp  /* detectors_overlap from last peak, matching Python fallback */
    );
}

/* ---- Module definition ---- */

static PyMethodDef rasterize_methods[] = {
    {"triangle_overlap", triangle_overlap, METH_VARARGS,
     "triangle_overlap(image, v0x, v0y, v1x, v1y, v2x, v2y) -> (overlap, lit)\n\n"
     "Fast triangle rasterization + overlap counting against experimental image.\n"
     "Implements Sutherland-Hodgman clipping + Bresenham scanline fill."},
    {"pixel_radius_overlap", pixel_radius_overlap, METH_VARARGS,
     "pixel_radius_overlap(image, cx, cy, radius) -> (in_bounds, bright)\n\n"
     "Check if any pixel in +-radius square around (cx, cy) is bright."},
    {"stage_d_overlap", stage_d_overlap, METH_VARARGS,
     "stage_d_overlap(image_list, wedge_indices, det_hit_mask, pixel_centers,\n"
     "                pixel_vertices, n_detectors, pixel_radius, n_peaks,\n"
     "                image_rows, image_cols)\n"
     "  -> (pixel_overlap, pixel_on_det, peak_overlap, peak_on_det,\n"
     "      quality, n_quality_points, detectors_overlap)\n\n"
     "Batch Stage D overlap: processes all peaks in a single C call.\n"
     "Uses uint8 binary images from ImageData.get_binary_numpy()."},
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
    return PyModule_Create(&rasterize_module);
}
