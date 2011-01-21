// vines: R package for multivariate dependence modeling with vines
// Copyright (C) 2010, 2011 Yasser González-Fernández
// Copyright (C) 2010, 2011 Marta Soto
//
// This program is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program. If not, see <http://www.gnu.org/licenses/>.

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include "h.h"

SEXP hNormalCopula(SEXP Rho, SEXP X, SEXP V)
{
	double eps = R_pow(DOUBLE_EPS, 0.5);
	double rho = NUMERIC_VALUE(Rho);
	int n = LENGTH(X);
	double *x, *v, *r;
	double vi, ri;
	SEXP R;

	x = NUMERIC_POINTER(X);
	v = NUMERIC_POINTER(V);
	PROTECT(R = NEW_NUMERIC(n));
	r = NUMERIC_POINTER(R);

	for (int i = 0; i < n; i++) {
		if (x[i] <= eps || (rho == 1 && x[i] == v[i] && x[i] != 1)) {
			r[i] = eps;
		} else if (1 - x[i] <= eps || (rho == -1 && 1 - (x[i] + v[i]) <= eps)) {
			r[i] = 1 - eps;
		} else {
			vi = (v[i] <= eps) ? eps : ((v[i] >= 1 - eps) ? 1 - eps : v[i]);
			ri = pnorm((qnorm(x[i], 0, 1, TRUE, FALSE) - rho*qnorm(vi, 0, 1, TRUE, FALSE)) /
					   sqrt(1 - rho*rho), 0, 1, TRUE, FALSE);
			r[i] = (ri <= eps) ? eps : ((ri >= 1 - eps) ? 1 - eps : ri);
		}
	}

	UNPROTECT(1);

	return R;
}

SEXP hIndepCopula(SEXP X, SEXP V)
{
	return X;
}

SEXP hTCopula(SEXP Rho, SEXP Df, SEXP X, SEXP V)
{
	double eps = R_pow(DOUBLE_EPS, 0.5);
	double rho = NUMERIC_VALUE(Rho);
	double df = NUMERIC_VALUE(Df);
	int n = LENGTH(X);
	double *x, *v, *r;
	double b2, vi, ri;
	SEXP R;

	x = NUMERIC_POINTER(X);
	v = NUMERIC_POINTER(V);
	PROTECT(R = NEW_NUMERIC(n));
	r = NUMERIC_POINTER(R);

	for (int i = 0; i < n; i++) {
		if (x[i] <= eps || (rho == 1 && x[i] == v[i] && x[i] != 1)) {
			r[i] = eps;
		} else if (1 - x[i] <= eps || (rho == -1 && 1 - (x[i] + v[i]) <= eps)) {
			r[i] = 1 - eps;
		} else {
			vi = (v[i] <= eps) ? eps : ((v[i] >= 1 - eps) ? 1 - eps : v[i]);
			b2 = qt(vi, df, TRUE, FALSE);
			ri = pt((qt(x[i], df, TRUE, FALSE) - rho*b2) /
					sqrt(((df + b2*b2) * (1 - rho*rho)) / (df+1)), df+1, TRUE, FALSE);
			r[i] = (ri <= eps) ? eps : ((ri >= 1 - eps) ? 1 - eps : ri);
		}
	}

	UNPROTECT(1);

	return R;
}

SEXP hClaytonCopula(SEXP Theta, SEXP X, SEXP V)
{
	double eps = R_pow(DOUBLE_EPS, 0.15);
	double theta = NUMERIC_VALUE(Theta);
	int n = LENGTH(X);
	double *x, *v, *r;
	double xi, vi, ri;
	SEXP R;

	if (theta <= eps) {
		return X;
	} else {
		x = NUMERIC_POINTER(X);
		v = NUMERIC_POINTER(V);
		PROTECT(R = NEW_NUMERIC(n));
		r = NUMERIC_POINTER(R);

		for (int i = 0; i < n; i++) {
			if (x[i] <= eps) {
				r[i] = eps;
			} else if (1 - x[i] <= eps) {
				r[i] = 1 - eps;
			} else {
				xi = (x[i] <= eps) ? eps : ((x[i] >= 1 - eps) ? 1 - eps : x[i]);
				vi = (v[i] <= eps) ? eps : ((v[i] >= 1 - eps) ? 1 - eps : v[i]);
				ri = R_pow(vi, -theta-1) * R_pow(R_pow(xi, -theta) + R_pow(vi, -theta) - 1, -1-1/theta);
				r[i] = (ri <= eps) ? eps : ((ri >= 1 - eps) ? 1 - eps : ri);
			}
		}

		UNPROTECT(1);

		return R;
	}
}
