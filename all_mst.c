/* Copyright 2024 Emmanuel Paradis */

#include <R.h>
#include <Rinternals.h>

double det_logscale(double *x, int n, int p, double *sign);
double det(double *x, int n, int p);

#define STACK_SIZE 1000000
#define EPS 1e-8 /* for the log of determinant */

void _getIandJ_(int k, int n, int *i, int *j)
{
    /* algo taken from pegas 1.2 */
    double b, x, y;
    k++;

    b = n - 0.5;
    x = ceil(b - sqrt(b * b - 2 * k));
    y = n * (1 - x) + (x + 1) * x/2 + k;
    *i = (int)x - 1;
    *j = (int)y - 1;
}

typedef struct DATA {
    double x;
    int o;
} DATA;

static int comp__(const void *a, const void *b)
{
    double x, y;
    x = ((struct DATA *)a)->x;
    y = ((struct DATA *)b)->x;;
    return (x > y) - (x < y);
}

void order_(double *x, int n, int *o)
{
    int i;
    struct DATA *X;

    X = (DATA*)R_alloc(n, sizeof(DATA));

    for (i = 0; i < n; i++) {
	X[i].x = x[i];
	X[i].o = i;
    }

    qsort(X, n, sizeof(struct DATA), comp__);

    for (i = 0; i < n; i++) o[i] = X[i].o;
}

int get_unique_values(int *x, int n, int *u, int N)
{
    int i, j, seen, k = 0;

    if (!N) {
	u[0] = x[0];
	N = 1;
	k = 1;
    }

    for (i = k; i < n; i++) {
	seen = 0;
	for (j = 0; j < N; j++) {
	    if (x[i] == u[j]) {
		seen = 1;
		break;
	    }
	}
	if (seen) continue;
	u[N] = x[i];
	N++;
    }
    return N;
}

int match_(int *x, int n, int v)
{
    int i;

    for (i = 0; i < n; i++)
	if (x[i] == v) return i;

    return -1;
}

int get_clusters(int *MAT, int Npairs, int *clusters, int *nnode, int *u)
{
    int i, j, f, fb, *clus, N, done, *MAT2;

    clus = (int*)R_alloc(STACK_SIZE, sizeof(int));

    MAT2 = MAT + STACK_SIZE;

    N = get_unique_values(MAT, Npairs, u, 0);
    *nnode = get_unique_values(MAT2, Npairs, u, N);

    //    for (i=0;i<*nnode;i++)Rprintf("u(%d) = %d\n", i, u[i]);

    for (i = 0; i < *nnode; i++) clus[i] = i + 1;

    for (i = 0; i < Npairs; i++) {
	/* Rprintf("MAT(%d) : %d %d\n", i, MAT[i], MAT2[i]); */
	f = clus[match_(u, *nnode, MAT[i])];
	fb = clus[match_(u, *nnode, MAT2[i])];
	if (f != fb) {
	    for (j = 0; j < *nnode; j++)
		if (clus[j] == fb) clus[j] = f;
	}
    }
    /* how many clusters? */
    int *u2;
    u2 = (int*)R_alloc(*nnode, sizeof(int));
    N = get_unique_values(clus, *nnode, u2, 0);

    //    for(i=0;i<*nnode;i++) Rprintf("clus[%d]=%d\n", i, clus[i]);

    //    Rprintf("nnode = %d   N = %d\n", *nnode, N);

    for (i = 0; i < *nnode; i++) {
	f = clus[i];
	done = 0;
	for (j = 0; j < N; j++) {
	    if (f == u2[j]) {
		clusters[i] = j;
		done = 1;
		break;
	    }
	    if (done) break;
	}
    }

    return N;
}

double get_det_Laplacian
(int *w_L, int n, int *MAT, int *MAT2, int Npairs, int nnode,
 int *u, int *clusters, int i, int log_scale, double *sign)
{
    double *L, res;
    int j, iw, i2, i3, i4;

    /* 1) build the Laplacian */

    L = (double*)R_Calloc(n * n, double);
    // L = (double*)R_alloc(n * n, sizeof(double));
    memset(L, 0, n * n * sizeof(double));
    memset(w_L, 0, n * sizeof(double));
    iw = 0;
    for (j = 0; j < Npairs; j++) {
	i2 = match_(u, nnode, MAT[j]);
	i3 = match_(u, nnode, MAT2[j]);
	if (clusters[i2] != i) continue;
	if (clusters[i3] != i) continue;
	i4 = i2 - 1;
	while (i4 >= 0) {
	    if (clusters[i4] != i) i2--;
	    i4--;
	}
	i4 = i3 - 1;
	while (i4 >= 0) {
	    if (clusters[i4] != i) i3--;
	    i4--;
	}
	if (i2 >= 0 && i3 >= 0) {
	    L[i2 + i3 * n] = L[i3 + i2 * n] = -1;
	    (L[i2 + i2 * n])++;
	    (L[i3 + i3 * n])++;
	    w_L[iw] = j; /* equivalent to K in the R code */
	    iw++;
	}
    }

    /* for (int iii = 0; iii < cluster_size; iii++) {
       for (int jjj = 0; jjj < cluster_size; jjj++)
       Rprintf("%.0f  ", L[iii + jjj * cluster_size]);
       Rprintf("\n");
       }
       Rprintf("\n"); */

    /* 2) drop the last row and the last column */
    i3 = 1; /* the upward shift to drop the last row */
    i2 = n - 1; /* update the number of cols */
    /* handle the 1st value first because of the modulo
       operation below */
    L[i2] = L[i2 + 1];
    for (j = i2 + 1; j < i2 * i2; j++) {
	L[j] = L[j + i3];
	if (!((j+1) % i2)) i3++;
    }

    /* for (int iii = 0; iii < i2; iii++) {
       for (int jjj = 0; jjj < i2; jjj++)
       Rprintf("%.0f ", L[iii + jjj * i2]);
       Rprintf("\n");
       } */

    /* 3) get the determinant */
    if (log_scale) {
	res = det_logscale(L, i2, i2, sign);
    } else {
	res = det(L, i2, i2);
    }
    R_Free(L);
    return res;
}

SEXP all_mst_C(SEXP D, SEXP size, SEXP LOG)
{
    int e = 0, i, j, k, n, l, z, Nedge, *o, *forest, f1, f2, *v1, *v2, log_scale, Nclusters, nnode;
    double *d, curr_d, *ll, n_MST;
    SEXP res, vtx1, vtx2, link_length;
    unsigned char *alt;

    int MAT[2 * STACK_SIZE], *MAT2;
    MAT2 = MAT + STACK_SIZE;
    int *w, *w_L, *clusters, *u, *rle;
    int Npairs,  pair, imat, discount, i2, i3;
    int cluster_size;
    double  detL, sign;

    PROTECT(D = coerceVector(D, REALSXP));
    PROTECT(size = coerceVector(size, INTSXP));
    PROTECT(LOG = coerceVector(LOG, INTSXP));
    d = REAL(D);
    l = LENGTH(D);
    n = INTEGER(size)[0];
    log_scale = INTEGER(LOG)[0];

    n_MST = !log_scale;

    Nedge = (n - 1) * n / 2; /* worst case */

    v1 = (int*)R_alloc(Nedge, sizeof(int));
    v2 = (int*)R_alloc(Nedge, sizeof(int));
    ll = (double*)R_alloc(Nedge, sizeof(double));

    w = (int*)R_alloc(STACK_SIZE, sizeof(int));
    w_L = (int*)R_alloc(STACK_SIZE, sizeof(int));
    u = (int*)R_alloc(STACK_SIZE, sizeof(int));
    clusters = (int*)R_alloc(STACK_SIZE, sizeof(int));

    alt = (unsigned char*)R_alloc(Nedge, sizeof(unsigned char));
    memset(alt, 0, Nedge * sizeof(unsigned char));

    forest = (int*)R_alloc(n, sizeof(int));
    for (i = 0; i < n; i++) forest[i] = i + 1;

    o = (int*)R_alloc(l, sizeof(int));
    order_(d, l, o);

    /* do RLE of the distances */
    rle = (int*)R_alloc(l, sizeof(int));
    i = k = 0;
    for (;;) {
	curr_d = d[o[i]];
	j = i;
	while (j < l - 1 && d[o[j + 1]] == curr_d) j++;
	rle[k] = j;
	k++;
	i = j + 1;
	if (i >= l) break;
    }

    //    for(i=0; i<k;i++) Rprintf("rle[%d] = %d\n", i, rle[i]);
    //    Rprintf("k=%d  l=%d\n", k, l);

    if (k == l) {
	/* error("no ties in distances: use mst()"); */
	n_MST = 1;
	goto end;
    }

    Npairs = rle[0] + 1;

    //    Rprintf("Npairs = %d\n", Npairs);

    k = 0; /* the next distance to look at */
    e = 0; /* the number of edges done */

    for (;;) {
	curr_d = d[o[k]];
	//	Rprintf("==>> k = %d\tcurr_d = %f\tNpairs = %d\n", k, curr_d, Npairs);
	discount = imat = 0;
	for (pair = 0; pair < Npairs; pair++) {
	    _getIandJ_(o[k], n, &i, &j);
	    if (j > n - 1 || i > n - 1) Rprintf("o[%d]=%d  i=%d  j=%d (n=%d)", k, o[k], i, j, n);
	    f1 = forest[i];
	    f2 = forest[j];
	    if (f1 == f2) {
		discount++;
	    } else {
		v1[e] = i + 1;
		v2[e] = j + 1;
		ll[e] = curr_d;
		e++;
		if (f1 < f2) {
		    MAT[imat] = f1;
		    MAT2[imat] = f2;
		} else {
		    MAT[imat] = f2;
		    MAT2[imat] = f1;
		}
		imat++;
	    }
	    k++;
	}

	Npairs -= discount;

	//	Rprintf("Npairs = %d\tdiscount = %d\n", Npairs, discount);

	if (Npairs < 1) goto here;

	//	for (i2 = 0; i2 < Npairs; i2++)
	//	      Rprintf("MAT[%d] = %d\tMAT2[%d] = %d\n", i2, MAT[i2], i2, MAT2[i2]);

	if (Npairs > 1) {
	    for (i2 = 0; i2 < Npairs; i2++) w[i2] = 1;
	    for (i2 = 0; i2 < Npairs - 1; i2++) {
		for (i3 = i2 + 1; i3 < Npairs; i3++) {
		    if (MAT[i2] == MAT[i3] && MAT2[i2] == MAT2[i3]) {
			/* <FIXME> this line works but not always (not with the woodmouse data) */
			/* It works with the TSPLIB's d15112 data. */
			/*     see also comment in the R code */
			alt[e - Npairs + i3] = 1; /* mark the link as alternative */
			/* </FIXME> */
			w[i2] += w[i3];
			w[i3] = 0;
		    }
		}
	    }

	    //	    for (i2 = 0; i2 < Npairs; i2++)
	    //		Rprintf("w[%d] = %d\n", i2, w[i2]);

	    /* "compress" the matrix of nodes */
	    i3 = 0;
	    for (i2 = 0; i2 < Npairs; i2++) {
		if (w[i2] && !i3) continue;
		if (!w[i2]) {
		    i3++;
		    continue;
		}
		MAT[i2 - i3] = MAT[i2];
		MAT2[i2 - i3] = MAT2[i2];
		w[i2 - i3] = w[i2];
	    }
	    Npairs -= i3;

	    /* get the clusters (independent aggregations */
	    Nclusters = get_clusters(MAT, Npairs, clusters, &nnode, u);

	    /* convert the memberships: */
	    for (i = 0; i < Nclusters; i++) {
		i2 = 0;
		while (clusters[i2] != i) i2++;
		f1 = u[i2];
		for (z = 0; z < n; z++) {
		    f2 = forest[z];
		    for (i3 = 0; i3 < nnode; i3++) {
			if (clusters[i3] != i) continue;
			if (u[i3] == f2) {
			    forest[z] = f1;
			    break;
			}
		    }
		}
	    }
	} else {
	    f1 = MAT[0];
	    f2 = MAT2[0];
	    for (z = 0; z < n; z++)
		if (forest[z] == f2) forest[z] = f1;
	    goto here;
	}

	for (i = 0; i < Nclusters; i++) {
	    cluster_size = 0;
	    for (j = 0; j < nnode; j++) {
		//		Rprintf("clusters[%d] = %d\n", j, clusters[j]);
		if (clusters[j] == i) cluster_size++;
	    }

	    //	    Rprintf("cluster_size = %d\tcurr_d = %f\tNclusters = %d\n", cluster_size, curr_d, Nclusters);
	    //	    for (int iii = 0; iii < Npairs; iii++) Rprintf("w[%d] = %d\n", iii, w[iii]);

	    if (cluster_size < 3) {
		for (j = 0; j < Npairs; j++) {
		    i2 = match_(u, nnode, MAT[j]);
		    i3 = match_(u, nnode, MAT2[j]);
		    if (clusters[i2] != i) continue;
		    if (clusters[i3] != i) continue;
		    if (i2 < 0) continue;
		    if (i3 < 0) continue;
		    /* see below (equivalent to K in the R code) */
		    if (log_scale) {
			n_MST += log(w[j]);
		    } else {
			n_MST *= w[j];
		    }
		}
		continue;
	    }

	    detL = get_det_Laplacian(w_L, cluster_size, MAT, MAT2,
				     Npairs, nnode, u, clusters, i,
				     log_scale, &sign);

	    if (log_scale) {
		//		Rprintf("detL = %.10f\tsign = %.0f\n", detL, sign);
		if (fabs(detL) < EPS) {
		    for (j = 0; j < cluster_size - 1; j++)
			n_MST += log(w[w_L[j]]);
		} else { /* there are loops */
		    /* check if some w > 1 */
		    int sum_w = 0;
		    for (j = 0; j < cluster_size; j++)
			sum_w += w[w_L[j]] - 1;
		    if (sum_w) {
			double ln_b = log(sum_w);
			detL += log1p(exp(ln_b - detL));
		    }
		    n_MST += sign * detL;
		}
	    } else {
		detL = round(detL);
		if (!detL) continue;
		if (detL == 1) {
		    for (j = 0; j < cluster_size - 1; j++)
			n_MST *= w[w_L[j]];
		} else {
		    for (j = 0; j < cluster_size; j++) {
			detL += (w[w_L[j]] - 1);
		    }
		    n_MST *= detL;
		}
	    }
	}
here:
    	/* for (int iii=0;iii<n;iii++)
	     Rprintf("forest[%d] = %d\n", iii, forest[iii]); */

	for (i2 = 1; i2 < n; i2++) {
	    if (forest[i2] != forest[0]) goto there;
	}
	break;
there:	if (e >= Nedge) break;

    	k = rle[0] + 1;
	//	Rprintf("---- k = %d   o[k] = %d\n", k, o[k]);
	if (k >= l) break;
	Npairs = rle[1] - rle[0];
	rle++;
    }
end:
    PROTECT(res = allocVector(VECSXP, 5));

    if (n_MST > 1) {
	PROTECT(vtx1 = allocVector(INTSXP, e));
	memcpy(INTEGER(vtx1), v1, e * sizeof(int));
	SET_VECTOR_ELT(res, 0, vtx1);

	PROTECT(vtx2 = allocVector(INTSXP, e));
	memcpy(INTEGER(vtx2), v2, e * sizeof(int));
	SET_VECTOR_ELT(res, 1, vtx2);

	PROTECT(link_length = allocVector(REALSXP, e));
	memcpy(REAL(link_length), ll, e * sizeof(double));
	SET_VECTOR_ELT(res, 2, link_length);

	SEXP ALT;
	PROTECT(ALT = allocVector(RAWSXP, e));
	memcpy(RAW(ALT), alt, e * sizeof(unsigned char));
	SET_VECTOR_ELT(res, 4, ALT);
    }

    SEXP Number_of_MST;
    PROTECT(Number_of_MST = allocVector(REALSXP, 1));
    REAL(Number_of_MST)[0] = n_MST;
    SET_VECTOR_ELT(res, 3, Number_of_MST);

    n_MST > 1 ? UNPROTECT(9) : UNPROTECT(5);
    return res;
}
