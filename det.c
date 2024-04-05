/* Copyright 2024 Emmanuel Paradis */

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

/* this version does det(x) on a log-scale for x square and symmetric */
SEXP bar(SEXP x)/* , SEXP logscale) */
{
    int n;
    double sign[1], *w;
    SEXP res;
    PROTECT(x = coerceVector(x, REALSXP));
    n = nrows(x);
    if (n != ncols(x)) error("matrix not square");
    w = (double *)R_alloc(n * n, sizeof(double));
    /* could use duplicate() */
    memcpy(w, REAL(x), n * n * sizeof(double));
    PROTECT(res = allocVector(REALSXP, 2));
    REAL(res)[0] = det_logscale(w, n, n, sign);
    REAL(res)[1] = sign[0];
    UNPROTECT(2);
    return res;
}
