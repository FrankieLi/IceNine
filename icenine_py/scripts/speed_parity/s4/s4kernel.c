/* S4 scratch kernels: the table-driven level-0 screen in several loop orders.
 *
 * Compiled by s4kernel.py with the system cc (-O3 -ffp-contract=off) and loaded with ctypes; not
 * part of the library. Semantics are those of screen_kernel.py (Step 2): for peak p and voxel
 * centroid (x, y) of one layer (z folded into c), u = a0 x + a1 y + c per axis in float64, the
 * pixel index is the float -> int truncation of the screen, the detector is lit when the
 * (2r+1)^2 box around that pixel meets the image and bright when the dilated bitmap is set there.
 * Rule 0 (the default qualified-peak rule with 2 detectors): a peak counts iff it is bright on
 * detector 0 (detector 1 cannot change it). Rule 1 (C++ CountQualifiedPeaks, 2 detectors): bright
 * on detector 0 and (bright on detector 1 or not lit on detector 1); detector 1 is tested only
 * when detector 0 is bright.
 *
 * Bitmaps: one byte per pixel (packed = 0) or one bit per pixel, most significant bit first
 * (packed = 1, numpy packbits), shape (n_frames, rows + 2r, stride), padded by r on every side.
 */
#include <stdint.h>
#include <string.h>
#include <math.h>

#pragma clang fp contract(off)

typedef struct {
    const double *a0, *a1, *c;   /* column: a0 x + a1 y + c */
    const double *b0, *b1, *d;   /* row */
    const uint8_t *hit;
    const uint8_t *bm;           /* frame 0 of this detector's bitmaps */
} Det;

typedef struct {
    Det det[2];
    const int32_t *frame;
    int64_t frame_stride;        /* bytes per frame */
    int32_t row_stride;          /* bytes per bitmap row */
    int32_t rows, cols, r, packed, rule;
} Table;

static inline int32_t to_idx(double u) {
    if (u < -1e7) u = -1e7;
    if (u > 1e7) u = 1e7;
    return (int32_t)u; /* truncation toward zero, as the screen's .astype(np.int32) */
}

/* 0: not lit, 1: lit and dark, 3: lit and bright */
static inline int probe(const Table *t, int d, int64_t p, double x, double y) {
    const Det *D = &t->det[d];
    if (!D->hit[p]) return 0;
    int32_t cx = to_idx(D->a0[p] * x + D->a1[p] * y + D->c[p]);
    int32_t cy = to_idx(D->b0[p] * x + D->b1[p] * y + D->d[p]);
    const int32_t r = t->r;
    if (cx < -r || cx >= t->cols + r || cy < -r || cy >= t->rows + r) return 0;
    const uint8_t *fb = D->bm + (int64_t)t->frame[p] * t->frame_stride;
    int32_t X = cx + r, Y = cy + r;
    int b;
    if (t->packed)
        b = (fb[(int64_t)Y * t->row_stride + (X >> 3)] >> (7 - (X & 7))) & 1;
    else
        b = fb[(int64_t)Y * t->row_stride + X] != 0;
    return b ? 3 : 1;
}

static inline int qualified(const Table *t, int64_t p, double x, double y) {
    int s0 = probe(t, 0, p, x, y);
    if (s0 != 3) return 0;
    if (t->rule == 0) return 1;
    int s1 = probe(t, 1, p, x, y);
    return (s1 == 3) || (s1 == 0);
}

/* O2: orientation-major. For each orientation of [o0, o1) and each of its peaks, all voxels of
   [v0, v1). out[(o - o0) * ld + (v - v0)] += qualified. */
void o2_counts(const Table *t, const int64_t *pstart, int64_t o0, int64_t o1, const double *xs,
               const double *ys, int64_t v0, int64_t v1, uint8_t *out, int64_t ld) {
    for (int64_t o = o0; o < o1; o++) {
        uint8_t *row = out + (o - o0) * ld - v0;
        for (int64_t p = pstart[o]; p < pstart[o + 1]; p++) {
            if (!t->det[0].hit[p]) continue;
            const double a0 = t->det[0].a0[p], a1 = t->det[0].a1[p], c = t->det[0].c[p];
            const double b0 = t->det[0].b0[p], b1 = t->det[0].b1[p], d = t->det[0].d[p];
            const uint8_t *fb = t->det[0].bm + (int64_t)t->frame[p] * t->frame_stride;
            const int32_t r = t->r, rs = t->row_stride, cmax = t->cols + r, rmax = t->rows + r;
            if (t->rule == 0 && !t->packed) {
                for (int64_t v = v0; v < v1; v++) {
                    int32_t cx = to_idx(a0 * xs[v] + a1 * ys[v] + c);
                    int32_t cy = to_idx(b0 * xs[v] + b1 * ys[v] + d);
                    if (cx < -r || cx >= cmax || cy < -r || cy >= rmax) continue;
                    row[v] += fb[(int64_t)(cy + r) * rs + (cx + r)] != 0;
                }
            } else if (t->rule == 0) {
                for (int64_t v = v0; v < v1; v++) {
                    int32_t cx = to_idx(a0 * xs[v] + a1 * ys[v] + c);
                    int32_t cy = to_idx(b0 * xs[v] + b1 * ys[v] + d);
                    if (cx < -r || cx >= cmax || cy < -r || cy >= rmax) continue;
                    int32_t X = cx + r;
                    row[v] += (fb[(int64_t)(cy + r) * rs + (X >> 3)] >> (7 - (X & 7))) & 1;
                }
            } else {
                for (int64_t v = v0; v < v1; v++) row[v] += qualified(t, p, xs[v], ys[v]);
            }
        }
    }
}

/* O2, pass only: out[(o - o0) * ld + (v - v0)] = 1 if any peak qualifies; each orientation keeps
   a list of the voxels not yet passed and stops when it is empty. alive: scratch of v1 - v0. */
void o2_pass(const Table *t, const int64_t *pstart, int64_t o0, int64_t o1, const double *xs,
             const double *ys, int64_t v0, int64_t v1, uint8_t *out, int64_t ld, int32_t *alive) {
    for (int64_t o = o0; o < o1; o++) {
        uint8_t *row = out + (o - o0) * ld;
        int64_t n = v1 - v0;
        for (int64_t k = 0; k < n; k++) alive[k] = (int32_t)(v0 + k);
        for (int64_t p = pstart[o]; p < pstart[o + 1] && n > 0; p++) {
            if (!t->det[0].hit[p]) continue;
            int64_t k = 0;
            while (k < n) {
                int32_t v = alive[k];
                if (qualified(t, p, xs[v], ys[v])) {
                    row[v - v0] = 1;
                    alive[k] = alive[--n];
                } else {
                    k++;
                }
            }
        }
    }
}

/* O1: voxel-major (today's order, table-driven). For each voxel, every orientation and peak. */
void o1_counts(const Table *t, const int64_t *pstart, int64_t o0, int64_t o1, const double *xs,
               const double *ys, int64_t v0, int64_t v1, uint8_t *out, int64_t ld) {
    for (int64_t v = v0; v < v1; v++) {
        const double x = xs[v], y = ys[v];
        for (int64_t o = o0; o < o1; o++) {
            int cnt = 0;
            for (int64_t p = pstart[o]; p < pstart[o + 1]; p++) cnt += qualified(t, p, x, y);
            out[(o - o0) * ld + (v - v0)] = (uint8_t)cnt;
        }
    }
}

/* O1, pass only: stop at the first qualified peak of each (orientation, voxel). */
void o1_pass(const Table *t, const int64_t *pstart, int64_t o0, int64_t o1, const double *xs,
             const double *ys, int64_t v0, int64_t v1, uint8_t *out, int64_t ld) {
    for (int64_t v = v0; v < v1; v++) {
        const double x = xs[v], y = ys[v];
        for (int64_t o = o0; o < o1; o++) {
            uint8_t ok = 0;
            for (int64_t p = pstart[o]; p < pstart[o + 1]; p++)
                if (qualified(t, p, x, y)) {
                    ok = 1;
                    break;
                }
            out[(o - o0) * ld + (v - v0)] = ok;
        }
    }
}

/* O3: frame-major. The table t must be sorted by frame (fstart[f] .. fstart[f + 1] the peaks of
   frame f, seg[p] their orientation, out rows indexed by orientation - o0 with all orientations of
   those peaks inside [o0, ...)). Voxel tiles of `tile` inside [v0, v1): for each tile, every frame,
   every peak of the frame, every voxel of the tile. */
void o3_counts(const Table *t, const int64_t *fstart, int64_t f0, int64_t f1, const int32_t *seg,
               int64_t o0, const double *xs, const double *ys, int64_t v0, int64_t v1,
               int64_t tile, uint8_t *out, int64_t ld) {
    for (int64_t vt = v0; vt < v1; vt += tile) {
        int64_t ve = vt + tile < v1 ? vt + tile : v1;
        for (int64_t f = f0; f < f1; f++) {
            for (int64_t p = fstart[f]; p < fstart[f + 1]; p++) {
                uint8_t *row = out + ((int64_t)seg[p] - o0) * ld - v0;
                if (!t->det[0].hit[p]) continue;
                const double a0 = t->det[0].a0[p], a1 = t->det[0].a1[p], c = t->det[0].c[p];
                const double b0 = t->det[0].b0[p], b1 = t->det[0].b1[p], d = t->det[0].d[p];
                const uint8_t *fb = t->det[0].bm + (int64_t)t->frame[p] * t->frame_stride;
                const int32_t r = t->r, rs = t->row_stride, cmax = t->cols + r, rmax = t->rows + r;
                if (t->rule == 0 && !t->packed) {
                    for (int64_t v = vt; v < ve; v++) {
                        int32_t cx = to_idx(a0 * xs[v] + a1 * ys[v] + c);
                        int32_t cy = to_idx(b0 * xs[v] + b1 * ys[v] + d);
                        if (cx < -r || cx >= cmax || cy < -r || cy >= rmax) continue;
                        row[v] += fb[(int64_t)(cy + r) * rs + (cx + r)] != 0;
                    }
                } else if (t->rule == 0) {
                    for (int64_t v = vt; v < ve; v++) {
                        int32_t cx = to_idx(a0 * xs[v] + a1 * ys[v] + c);
                        int32_t cy = to_idx(b0 * xs[v] + b1 * ys[v] + d);
                        if (cx < -r || cx >= cmax || cy < -r || cy >= rmax) continue;
                        int32_t X = cx + r;
                        row[v] += (fb[(int64_t)(cy + r) * rs + (X >> 3)] >> (7 - (X & 7))) & 1;
                    }
                } else {
                    for (int64_t v = vt; v < ve; v++) row[v] += qualified(t, p, xs[v], ys[v]);
                }
            }
        }
    }
}

/* O5: orientation block x voxel tile (tiled O2): for each voxel tile of [v0, v1), the orientations
   of [o0, o1). */
void o5_counts(const Table *t, const int64_t *pstart, int64_t o0, int64_t o1, const double *xs,
               const double *ys, int64_t v0, int64_t v1, int64_t tile, uint8_t *out, int64_t ld) {
    for (int64_t vt = v0; vt < v1; vt += tile) {
        int64_t ve = vt + tile < v1 ? vt + tile : v1;
        o2_counts(t, pstart, o0, o1, xs, ys, vt, ve, out + (vt - v0), ld);
    }
}

/* O4: pixel voting (rule 0, detector 0, packed bitmaps). The voxel centroids lie on two
   triangular lattices (one per triangle type): centroid ~ org[t] + i e1 + j e2, with grid[t] the
   local voxel index at (i - i0[t], j - j0[t]) (-1: none). For each peak, the footprint of the
   patch's bounding box (bx0..bx1, by0..by1 in the sample plane) is scanned row by row in the
   packed dilated bitmap; each bright pixel is back-projected to the lattice points whose pixel
   can be it (a superset: the box of the pixel's preimage, widened by `slack` lattice units), and
   each candidate voxel is tested exactly with the O2 arithmetic. Same counts as O2. */
typedef struct {
    const int32_t *grid[2];
    int32_t i0[2], j0[2], ni[2], nj[2];
    double org[2][2];
    double e1[2], e2[2];
    double bx0, bx1, by0, by1;
    double slack;
} Lattice;

static inline void pre_interval(int32_t q, double *lo, double *hi) {
    /* the u with trunc(u) == q is a subset of [lo, hi] */
    if (q > 0) { *lo = q; *hi = q + 1; }
    else if (q < 0) { *lo = q - 1; *hi = q; }
    else { *lo = -1; *hi = 1; }
}

void o4_counts(const Table *t, const Lattice *L, const int64_t *pstart, int64_t o0, int64_t o1,
               const double *xs, const double *ys, uint8_t *out, int64_t ld, int64_t *stats) {
    const int32_t r = t->r, rs = t->row_stride;
    int64_t n_bright = 0, n_cand = 0;
    for (int64_t o = o0; o < o1; o++) {
        uint8_t *row_out = out + (o - o0) * ld;
        for (int64_t p = pstart[o]; p < pstart[o + 1]; p++) {
            const Det *D = &t->det[0];
            if (!D->hit[p]) continue;
            const double a0 = D->a0[p], a1 = D->a1[p], c = D->c[p];
            const double b0 = D->b0[p], b1 = D->b1[p], d = D->d[p];
            const uint8_t *fb = D->bm + (int64_t)t->frame[p] * t->frame_stride;
            /* footprint quadrilateral of the sample box */
            double qx[4], qy[4];
            const double sx[4] = {L->bx0, L->bx1, L->bx1, L->bx0};
            const double sy[4] = {L->by0, L->by0, L->by1, L->by1};
            double ymin = 1e300, ymax = -1e300;
            for (int k = 0; k < 4; k++) {
                qx[k] = a0 * sx[k] + a1 * sy[k] + c;
                qy[k] = b0 * sx[k] + b1 * sy[k] + d;
                if (qy[k] < ymin) ymin = qy[k];
                if (qy[k] > ymax) ymax = qy[k];
            }
            int32_t yr0 = (int32_t)floor(ymin) - 2, yr1 = (int32_t)floor(ymax) + 2;
            if (yr0 < -r) yr0 = -r;
            if (yr1 > t->rows + r - 1) yr1 = t->rows + r - 1;
            /* lattice maps: u = G (i, j) + u0[t] */
            const double g1x = a0 * L->e1[0] + a1 * L->e1[1], g1y = b0 * L->e1[0] + b1 * L->e1[1];
            const double g2x = a0 * L->e2[0] + a1 * L->e2[1], g2y = b0 * L->e2[0] + b1 * L->e2[1];
            const double det = g1x * g2y - g2x * g1y;
            if (det == 0) continue;
            const double h11 = g2y / det, h12 = -g2x / det, h21 = -g1y / det, h22 = g1x / det;
            double u0x[2], u0y[2];
            for (int tt = 0; tt < 2; tt++) {
                u0x[tt] = a0 * L->org[tt][0] + a1 * L->org[tt][1] + c;
                u0y[tt] = b0 * L->org[tt][0] + b1 * L->org[tt][1] + d;
            }
            for (int32_t y = yr0; y <= yr1; y++) {
                /* x-range of the quadrilateral on the band [y - 1, y + 2] (margin) */
                double xl = 1e300, xh = -1e300;
                const double ylo = y - 1.0, yhi = y + 2.0;
                for (int k = 0; k < 4; k++) {
                    int k2 = (k + 1) & 3;
                    double x1 = qx[k], y1 = qy[k], x2 = qx[k2], y2 = qy[k2];
                    if (y1 >= ylo && y1 <= yhi) { if (x1 < xl) xl = x1; if (x1 > xh) xh = x1; }
                    double ys_[2] = {ylo, yhi};
                    for (int e = 0; e < 2; e++) {
                        double yy = ys_[e];
                        if ((y1 - yy) * (y2 - yy) <= 0 && y1 != y2) {
                            double xx = x1 + (x2 - x1) * (yy - y1) / (y2 - y1);
                            if (xx < xl) xl = xx;
                            if (xx > xh) xh = xx;
                        }
                    }
                }
                if (xl > xh) continue;
                int32_t x0 = (int32_t)floor(xl) - 2, x1i = (int32_t)floor(xh) + 2;
                if (x0 < -r) x0 = -r;
                if (x1i > t->cols + r - 1) x1i = t->cols + r - 1;
                if (x0 > x1i) continue;
                const uint8_t *prow = fb + (int64_t)(y + r) * rs;
                int32_t X0 = x0 + r, X1 = x1i + r;
                for (int32_t B = X0 >> 3; B <= (X1 >> 3); B++) {
                    uint8_t byte = prow[B];
                    if (!byte) continue;
                    for (int bit = 0; bit < 8; bit++) {
                        if (!((byte >> (7 - bit)) & 1)) continue;
                        int32_t X = B * 8 + bit;
                        if (X < X0 || X > X1) continue;
                        int32_t qcol = X - r, qrow = y;
                        n_bright++;
                        double cl, ch, rl, rh;
                        pre_interval(qcol, &cl, &ch);
                        pre_interval(qrow, &rl, &rh);
                        for (int tt = 0; tt < 2; tt++) {
                            double imin = 1e300, imax = -1e300, jmin = 1e300, jmax = -1e300;
                            const double cu[4] = {cl, ch, ch, cl}, cv[4] = {rl, rl, rh, rh};
                            for (int k = 0; k < 4; k++) {
                                double du = cu[k] - u0x[tt], dv = cv[k] - u0y[tt];
                                double ii = h11 * du + h12 * dv, jj = h21 * du + h22 * dv;
                                if (ii < imin) imin = ii;
                                if (ii > imax) imax = ii;
                                if (jj < jmin) jmin = jj;
                                if (jj > jmax) jmax = jj;
                            }
                            int32_t ia = (int32_t)ceil(imin - L->slack) - L->i0[tt];
                            int32_t ib = (int32_t)floor(imax + L->slack) - L->i0[tt];
                            int32_t ja = (int32_t)ceil(jmin - L->slack) - L->j0[tt];
                            int32_t jb = (int32_t)floor(jmax + L->slack) - L->j0[tt];
                            if (ia < 0) ia = 0;
                            if (ja < 0) ja = 0;
                            if (ib >= L->ni[tt]) ib = L->ni[tt] - 1;
                            if (jb >= L->nj[tt]) jb = L->nj[tt] - 1;
                            for (int32_t i = ia; i <= ib; i++)
                                for (int32_t j = ja; j <= jb; j++) {
                                    int32_t v = L->grid[tt][(int64_t)i * L->nj[tt] + j];
                                    if (v < 0) continue;
                                    n_cand++;
                                    int32_t cx = to_idx(a0 * xs[v] + a1 * ys[v] + c);
                                    int32_t cy = to_idx(b0 * xs[v] + b1 * ys[v] + d);
                                    if (cx == qcol && cy == qrow) row_out[v]++;
                                }
                        }
                    }
                }
            }
        }
    }
    stats[0] += n_bright;
    stats[1] += n_cand;
}

/* size of the Table struct, checked from Python */
int64_t table_size(void) { return (int64_t)sizeof(Table); }
