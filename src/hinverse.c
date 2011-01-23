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

#include "hinverse.h"

SEXP hinverseNormalCopula(SEXP RHO, SEXP U, SEXP V) {
    double eps = R_pow(DOUBLE_EPS, 0.5);
    double rho = NUMERIC_VALUE(RHO);
    int n = LENGTH(U);
    double *u, *v, *x;
    double vi, xi;
    SEXP X;

    u = NUMERIC_POINTER(U);
    v = NUMERIC_POINTER(V);
    PROTECT(X = NEW_NUMERIC(n));
    x = NUMERIC_POINTER(X);

    for (int i = 0; i < n; i++) {
        if (u[i] <= eps) {
            x[i] = eps;
        } else if (1 - u[i] <= eps) {
            x[i] = 1 - eps;
        } else {
            vi = (v[i] <= eps) ? eps : ((v[i] >= 1 - eps) ? 1 - eps : v[i]);
            xi = pnorm(qnorm(u[i], 0, 1, TRUE, FALSE) * sqrt(1 - rho * rho) +
                    rho * qnorm(vi, 0, 1, 1, 0), 0, 1, TRUE, FALSE);
            x[i] = (xi <= eps) ? eps : ((xi >= 1 - eps) ? 1 - eps : xi);
        }
    }

    UNPROTECT(1);

    return X;
}

SEXP hinverseIndepCopula(SEXP U, SEXP V) {
    return U;
}

SEXP hinverseTCopula(SEXP RHO, SEXP DF, SEXP U, SEXP V) {
    double eps = R_pow(DOUBLE_EPS, 0.5);
    double rho = NUMERIC_VALUE(RHO);
    double df = NUMERIC_VALUE(DF);
    int n = LENGTH(U);
    double *u, *v, *x;
    double b2, vi, xi;
    SEXP X;

    u = NUMERIC_POINTER(U);
    v = NUMERIC_POINTER(V);
    PROTECT(X = NEW_NUMERIC(n));
    x = NUMERIC_POINTER(X);

    for (int i = 0; i < n; i++) {
        if (u[i] <= eps) {
            x[i] = eps;
        } else if (1 - u[i] <= eps) {
            x[i] = 1 - eps;
        } else {
            vi = (v[i] <= eps) ? eps : ((v[i] >= 1 - eps) ? 1 - eps : v[i]);
            b2 = qt(vi, df, TRUE, FALSE);
            xi = pt(qt(u[i], df + 1, TRUE, FALSE) *
                    sqrt(((df + b2 * b2) * (1 - rho * rho)) / (df + 1)) +
                    rho * b2, df, TRUE, FALSE);
            x[i] = (xi <= eps) ? eps : ((xi >= 1 - eps) ? 1 - eps : xi);
        }
    }

    UNPROTECT(1);

    return X;
}

SEXP hinverseClaytonCopula(SEXP THETA, SEXP U, SEXP V) {
    double eps = R_pow(DOUBLE_EPS, 0.15);
    double theta = NUMERIC_VALUE(THETA);
    int n = LENGTH(U);
    double *u, *v, *x;
    double vi, xi;
    SEXP X;

    if (theta <= eps) {
        X = U;
    } else {
        u = NUMERIC_POINTER(U);
        v = NUMERIC_POINTER(V);
        PROTECT(X = NEW_NUMERIC(n));
        x = NUMERIC_POINTER(X);

        for (int i = 0; i < n; i++) {
            if (u[i] <= eps) {
                x[i] = eps;
            } else if (1 - u[i] <= eps) {
                x[i] = 1 - eps;
            } else {
                vi = (v[i] <= eps) ? eps : ((v[i] >= 1 - eps) ? 1 - eps : v[i]);
                xi = R_pow(R_pow(u[i] * R_pow(vi, theta + 1), -theta/(theta+1)) +
                        1 - R_pow(vi, -theta), -1 / theta);
                x[i] = (xi <= eps) ? eps : ((xi >= 1 - eps) ? 1 - eps : xi);
            }
        }

        UNPROTECT(1);
    }

    return X;
}
