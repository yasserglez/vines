#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include "hinverse.h"

SEXP hinverseNormalCopula(SEXP Rho, SEXP U, SEXP V) {
	double eps = R_pow(DOUBLE_EPS, 0.5);
	double rho = NUMERIC_VALUE(Rho);
	int n = LENGTH(U);
	double *u, *v, *r;
	double vi, ri;
	SEXP R;

	u = NUMERIC_POINTER(U);
	v = NUMERIC_POINTER(V);
	PROTECT(R = NEW_NUMERIC(n));
	r = NUMERIC_POINTER(R);

	for (int i = 0; i < n; i++) {
		if (u[i] <= eps) {
			r[i] = eps;
		} else if (1 - u[i] <= eps) {
			r[i] = 1 - eps;
		} else {
			vi = (v[i] <= eps) ? eps : ((v[i] >= 1 - eps) ? 1 - eps : v[i]);
			ri = pnorm(qnorm(u[i], 0, 1, 1, 0) * sqrt(1 - rho*rho) + rho*qnorm(vi, 0, 1, 1, 0), 0, 1, 1, 0);
			r[i] = (ri <= eps) ? eps : ((ri >= 1 - eps) ? 1 - eps : ri);
		}
	}

	UNPROTECT(1);

	return R;
}
