/*
 * Sparse Gaussian elimination modulo a prime p < 2^63, used for the first numerical
 * probe of an IBP system.  It finds which equations the targets actually need, so the
 * later probes can run on a much smaller system.
 *
 * Rows are given in CSR form (rowptr, cols); the value of entry e is the integer
 * combination E[e, :] of the monomial values mv, taken mod p.
 * Column 0 is the most complicated integral; each row is solved for its smallest
 * column, as in the Python code in ff.py.
 *
 * Output: the indices of the rows that produced the pivots the targets depend on
 * (through elimination and back substitution), in increasing order.
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned __int128 u128;
typedef __int128 i128;

static uint64_t inv_mod(uint64_t a, uint64_t p) {
    i128 t = 0, nt = 1, r = p, nr = a;
    while (nr) {
        i128 q = r / nr, tmp;
        tmp = t - q * nt; t = nt; nt = tmp;
        tmp = r - q * nr; r = nr; nr = tmp;
    }
    if (t < 0) t += p;
    return (uint64_t)t;
}

typedef struct { int32_t *a; int n, cap; } ivec;
static void ipush(ivec *v, int32_t x) {
    if (v->n == v->cap) { v->cap = v->cap ? 2 * v->cap : 8; v->a = realloc(v->a, v->cap * sizeof(int32_t)); }
    v->a[v->n++] = x;
}

/* binary min-heap of column indices */
static void hpush(int32_t *h, int *n, int32_t x) {
    int i = (*n)++;
    while (i > 0) { int par = (i - 1) >> 1; if (h[par] <= x) break; h[i] = h[par]; i = par; }
    h[i] = x;
}
static int32_t hpop(int32_t *h, int *n) {
    int32_t top = h[0], x = h[--(*n)];
    int i = 0;
    for (;;) {
        int l = 2 * i + 1, r = l + 1, m = i;
        int32_t best = x;
        if (l < *n && h[l] < best) { m = l; best = h[l]; }
        if (r < *n && h[r] < best) { m = r; best = h[r]; }
        if (m == i) break;
        h[i] = h[m]; i = m;
    }
    if (*n > 0) h[i] = x;
    return top;
}

/* entry e has the value  sum_m E[e*nmon + m] * mv[m]  mod p  (mv = monomial values) */
static uint64_t entry_value(const int64_t *E, int nmon, const uint64_t *mv, uint64_t p, int64_t e) {
    u128 s = 0;
    for (int m = 0; m < nmon; m++) {
        int64_t x = E[e * nmon + m];
        if (!x) continue;
        uint64_t c = x >= 0 ? (uint64_t)x % p : (p - ((uint64_t)(-x) % p)) % p;
        s = (s + (u128)c * mv[m]) % p;
    }
    return (uint64_t)s;
}

int64_t elim_trim(int32_t nrows, int32_t ncols, const int64_t *rowptr, const int32_t *cols,
                  const int64_t *E, int32_t nmon, const uint64_t *mv, uint64_t p,
                  int32_t ntargets, const int32_t *targets, int32_t *out_rows) {
    uint64_t *acc = calloc(ncols, sizeof(uint64_t));
    char *inheap = calloc(ncols, 1);
    int32_t *heap = malloc(sizeof(int32_t) * (ncols + 1));
    int32_t *src = malloc(sizeof(int32_t) * ncols);          /* row that made pivot j, or -1 */
    int32_t **pc = calloc(ncols, sizeof(int32_t *));           /* pivot row columns (after j) */
    uint64_t **pv = calloc(ncols, sizeof(uint64_t *));
    int32_t *plen = calloc(ncols, sizeof(int32_t));
    ivec *hits = calloc(ncols, sizeof(ivec));
    ivec cur = {0, 0, 0};
    int32_t *tmpc = malloc(sizeof(int32_t) * ncols);
    for (int j = 0; j < ncols; j++) src[j] = -1;

    for (int32_t n = 0; n < nrows; n++) {
        int hn = 0;
        cur.n = 0;
        for (int64_t e = rowptr[n]; e < rowptr[n + 1]; e++) {
            int32_t c = cols[e];
            acc[c] = (acc[c] + entry_value(E, nmon, mv, p, e)) % p;
            if (!inheap[c]) { inheap[c] = 1; hpush(heap, &hn, c); }
        }
        while (hn > 0) {
            int32_t j = hpop(heap, &hn);
            inheap[j] = 0;
            uint64_t c = acc[j];
            if (!c) continue;
            if (src[j] >= 0) {                     /* eliminate with pivot j */
                acc[j] = 0;
                ipush(&cur, j);
                uint64_t mc = p - c;
                for (int t = 0; t < plen[j]; t++) {
                    int32_t k = pc[j][t];
                    acc[k] = (uint64_t)(((u128)mc * pv[j][t] + acc[k]) % p);
                    if (!inheap[k]) { inheap[k] = 1; hpush(heap, &hn, k); }
                }
                continue;
            }
            /* new pivot at j: the rest of the heap is the rest of the row */
            uint64_t inv = inv_mod(c, p);
            acc[j] = 0;
            int m = 0;
            while (hn > 0) {
                int32_t k = hpop(heap, &hn);
                inheap[k] = 0;
                if (acc[k]) tmpc[m++] = k;
            }
            pc[j] = malloc(sizeof(int32_t) * (m ? m : 1));
            pv[j] = malloc(sizeof(uint64_t) * (m ? m : 1));
            for (int t = 0; t < m; t++) {
                int32_t k = tmpc[t];
                pc[j][t] = k;
                pv[j][t] = (uint64_t)(((u128)acc[k] * inv) % p);
                acc[k] = 0;
            }
            plen[j] = m;
            src[j] = n;
            hits[j].a = NULL; hits[j].n = 0; hits[j].cap = 0;
            for (int t = 0; t < cur.n; t++) ipush(&hits[j], cur.a[t]);
            break;
        }
        /* a row that reduced to zero leaves acc clean; after a new pivot the heap is empty */
    }

    /* closure from the targets */
    char *need = calloc(ncols, 1);
    int32_t *stack = malloc(sizeof(int32_t) * ncols);
    int sn = 0;
    for (int t = 0; t < ntargets; t++) {
        int32_t j = targets[t];
        if (j >= 0 && j < ncols && src[j] >= 0 && !need[j]) { need[j] = 1; stack[sn++] = j; }
    }
    while (sn > 0) {
        int32_t j = stack[--sn];
        for (int t = 0; t < hits[j].n; t++) {
            int32_t k = hits[j].a[t];
            if (src[k] >= 0 && !need[k]) { need[k] = 1; stack[sn++] = k; }
        }
        for (int t = 0; t < plen[j]; t++) {
            int32_t k = pc[j][t];
            if (src[k] >= 0 && !need[k]) { need[k] = 1; stack[sn++] = k; }
        }
    }
    /* rows in increasing order */
    char *rowneed = calloc(nrows, 1);
    for (int j = 0; j < ncols; j++) if (need[j]) rowneed[src[j]] = 1;
    int64_t cnt = 0;
    for (int32_t n = 0; n < nrows; n++) if (rowneed[n]) out_rows[cnt++] = n;

    for (int j = 0; j < ncols; j++) { free(pc[j]); free(pv[j]); free(hits[j].a); }
    free(acc); free(inheap); free(heap); free(src); free(pc); free(pv); free(plen); free(hits);
    free(cur.a); free(tmpc); free(need); free(stack); free(rowneed);
    return cnt;
}
