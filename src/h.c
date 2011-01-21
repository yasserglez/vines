#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include "h.h"


SEXP hNormalCopula(SEXP Rho, SEXP X, SEXP V)
{
	double eps = R_pow(DOUBLE_EPS, 0.5);
	double rho = NUMERIC_VALUE(Rho);
	double *x, *v, *r;
	int i, n;
	double vi, ri;
	SEXP R;

	n = LENGTH(X);
	x = NUMERIC_POINTER(X);
	v = NUMERIC_POINTER(V);
	PROTECT(R = NEW_NUMERIC(n));
	r = NUMERIC_POINTER(R);

	for (i = 0; i < n; i++) {
		if (x[i] <= eps || (rho == 1 && x[i] == v[i] && x[i] != 1)) {
			r[i] = eps;
		} else if (1 - x[i] <= eps || (rho == -1 && 1 - (x[i] + v[i]) <= eps)) {
			r[i] = 1 - eps;
		} else {
			vi = (v[i] <= eps) ? eps : ((v[i] >= 1 - eps) ? 1 - eps : v[i]);
			ri = pnorm((qnorm(x[i], 0, 1, 1, 0) - rho * qnorm(vi, 0, 1, 1, 0)) /
					    sqrt(1 - (rho * rho)), 0, 1, 1, 0);
			r[i] = (ri <= eps) ? eps : ((ri >= 1 - eps) ? 1 - eps : ri);
		}
	}

	UNPROTECT(1);

	return R;
}
