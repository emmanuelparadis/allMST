/* Copyright 2024-2025 Emmanuel Paradis */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>

#define sign(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

/* QR decomposition of a rectangular n*p matrix */
int _qr_(double *x, int n, int p)
{
    double *work, *qraux, tol = 1e-8;
    int *jpiv, ldx, i, k;

    ldx = n >= p ? n : p;

    work = (double *)R_alloc(2 * p, sizeof(double));
    qraux = (double *)R_alloc(p, sizeof(double));

    jpiv = (int *)R_alloc(p, sizeof(int));
    for (i = 0; i < p; i++) jpiv[i] = i + 1;

    F77_CALL(dqrdc2)(x, &ldx, &n, &p, &tol, &k, qraux, jpiv, work);

    /* base::determinant() uses F77_CALL(dgetrf) which is much faster.
       see https://github.com/wch/r-source/blob/trunk/src/modules/lapack/Lapack.c
       line 1415.
       See next function. */

    return 0;
}

/* QR decomposition of a square symmetric matrix */
int qr_symm(double *x, int n)
{
    int *jpiv, info;
    jpiv = (int *)R_alloc(n, sizeof(int));
    F77_CALL(dgetrf)(&n, &n, x, &n, jpiv, &info);
    if (info < 0) error("error during QR decomp. of symm. matrix");
    return 0;
}

double det_logscale(double *x, int n, int p, double *sign)
{
    int o, i;
    double res = 0, thesign = 1;

    if (n != p) error("matrix not square");
    if (n == 1) return log(x[0]);

    /* o = _qr_(x, n, p); */
    o = qr_symm(x, n);
    if (o) return -1;

    for (i = 0; i < n*n; i += (n + 1)) {
	thesign *= sign(x[i]);
	if (!thesign) {
	    *sign = 0;
	    return 0;
	}
	res += log(fabs(x[i]));
    }

    if (!(n % 2)) res = -res;
    *sign = thesign;
    return res;
}

double det(double *x, int n, int p)
{
    int o, i;
    double res = 1;

    if (n != p) error("matrix not square");
    if (n == 1) return x[0];

    o = _qr_(x, n, p);
    if (o) return -1;

    for (i = 0; i < n*n; i += (n + 1)) res *= x[i];

    if (!(n % 2)) res = -res;
    return res;
}

void foo(double *x, int *n, int *p, double *res, int *logscale, double *sign)
{
    if (*logscale) {
	res[0] = det_logscale(x, *n, *p, sign);
    } else {
	res[0] = det(x, *n, *p);
    }
}

int any_row_filled_with_zeros(double *x, int n, int p)
{
    int i, k, end = (p - 1)*n;
    for (i = 0; i < n; i++) {
	for (k = 0; k <= end; k += n) {
	    if (x[k]) break;
	}
	if (k > end) return 1;
	x++;
    }
    return 0;
}

int any_column_filled_with_zeros(double *x, int n, int p)
{
    int i, j;
    for (j = 0; j < p; j++) {
	for (i = 0; i < n; i++) {
	    if (x[i]) break;
	}
	if (i == n) return 1;
	x += n;
    }
    return 0;
}

/* this version does det(x) on a log-scale for x square and symmetric */
SEXP bar(SEXP x)/* , SEXP logscale) */
{
    int n;
    double sign, *w, *p;
    SEXP res;
    PROTECT(x = coerceVector(x, REALSXP));
    n = nrows(x);
    if (n != ncols(x)) error("matrix not square");

    PROTECT(res = allocVector(REALSXP, 2));

    /* test if one row or column is filled with 0 */
    p = REAL(x);
    if (any_column_filled_with_zeros(p, n, n) || any_row_filled_with_zeros(p, n, n)) {
	p[0] = 0;
	p[1] = 0;
    } else {
	w = (double *)R_alloc(n * n, sizeof(double));
	memcpy(w, REAL(x), n * n * sizeof(double));
	p[0] = det_logscale(w, n, n, &sign);
	p[1] = sign;
    }

    UNPROTECT(2);
    return res;
}
