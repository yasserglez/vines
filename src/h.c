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
#include <Rinternals.h>
#include <Rmath.h>

SEXP hIndepCopula(SEXP X, SEXP V) {
    return X;
}

// See Aas, K., Czado, C., Frigessi, A. & Bakken, H. Pair-copula constructions
// of multiple dependence. Insurance Mathematics and Economics, 2009, Vol. 44,
// pp. 182-198 for the expression for the Gaussian, Student's t, Clayton and
// Gumbel copulas.

SEXP hNormalCopula(SEXP RHO, SEXP X, SEXP V) {
    double eps = R_pow(DOUBLE_EPS, 0.5);
    double rho = asReal(RHO);
    int n = LENGTH(X);
    double *x, *v, *u;
    double vi, ui;
    SEXP U;

    x = REAL(X);
    v = REAL(V);
    PROTECT(U = allocVector(REALSXP, n));
    u = REAL(U);

    for (int i = 0; i < n; i++) {
        if (x[i] <= eps || (rho == 1 && x[i] == v[i] && x[i] != 1)) {
            u[i] = eps;
        } else if (1 - x[i] <= eps || (rho == -1 && 1 - (x[i] + v[i]) <= eps)) {
            u[i] = 1 - eps;
        } else {
            vi = (v[i] <= eps) ? eps : ((v[i] >= 1 - eps) ? 1 - eps : v[i]);
            ui = pnorm((qnorm(x[i], 0, 1, TRUE, FALSE) -
                        rho * qnorm(vi, 0, 1, TRUE, FALSE)) /
                    sqrt(1 - rho * rho), 0, 1, TRUE, FALSE);
            u[i] = (ui <= eps) ? eps : ((ui >= 1 - eps) ? 1 - eps : ui);
        }
    }

    UNPROTECT(1);

    return U;
}

SEXP hTCopula(SEXP RHO, SEXP DF, SEXP X, SEXP V) {
    double eps = R_pow(DOUBLE_EPS, 0.5);
    double rho = asReal(RHO);
    double df = asReal(DF);
    int n = LENGTH(X);
    double *x, *v, *u;
    double vi, ui;
    double tmp;
    SEXP U;

    x = REAL(X);
    v = REAL(V);
    PROTECT(U = allocVector(REALSXP, n));
    u = REAL(U);

    for (int i = 0; i < n; i++) {
        if (x[i] <= eps || (rho == 1 && x[i] == v[i] && x[i] != 1)) {
            u[i] = eps;
        } else if (1 - x[i] <= eps || (rho == -1 && 1 - (x[i] + v[i]) <= eps)) {
            u[i] = 1 - eps;
        } else {
            vi = (v[i] <= eps) ? eps : ((v[i] >= 1 - eps) ? 1 - eps : v[i]);
            tmp = qt(vi, df, TRUE, FALSE);
            ui = pt((qt(x[i], df, TRUE, FALSE) - rho * tmp) /
                    sqrt(((df + tmp * tmp) * (1 - rho * rho)) / (df + 1)),
                    df + 1, TRUE, FALSE);
            u[i] = (ui <= eps) ? eps : ((ui >= 1 - eps) ? 1 - eps : ui);
        }
    }

    UNPROTECT(1);

    return U;
}

SEXP hClaytonCopula(SEXP THETA, SEXP X, SEXP V) {
    double eps = R_pow(DOUBLE_EPS, 0.15);
    double theta = asReal(THETA);
    int n = LENGTH(X);
    double *x, *v, *u;
    double vi, ui;
    SEXP U;

    if (theta <= eps) {
        U = X;
    } else {
        x = REAL(X);
        v = REAL(V);
        PROTECT(U = allocVector(REALSXP, n));
        u = REAL(U);

        for (int i = 0; i < n; i++) {
            if (x[i] <= eps) {
                u[i] = eps;
            } else if (1 - x[i] <= eps) {
                u[i] = 1 - eps;
            } else {
                vi = (v[i] <= eps) ? eps : ((v[i] >= 1 - eps) ? 1 - eps : v[i]);
                ui = R_pow(vi, -theta - 1) *
                        R_pow(R_pow(x[i], -theta) + R_pow(vi, -theta) - 1,
                                -1 - 1 / theta);
                u[i] = (ui <= eps) ? eps : ((ui >= 1 - eps) ? 1 - eps : ui);
            }
        }

        UNPROTECT(1);
    }

    return U;
}

SEXP hGumbelCopula(SEXP THETA, SEXP X, SEXP V)
{
    double eps = R_pow(DOUBLE_EPS, 0.5);
    double theta = asReal(THETA);
    int n = LENGTH(X);
    double *x, *v, *u;
    double vi, ui;
    double mlogxi, mlogvi, tmp;
    SEXP U;

    if (theta <= eps) {
        U = X;
    } else {
        x = REAL(X);
        v = REAL(V);
        PROTECT(U = allocVector(REALSXP, n));
        u = REAL(U);

        for (int i = 0; i < n; i++) {
            if (x[i] <= eps) {
                u[i] = eps;
            } else if (1 - x[i] <= eps) {
                u[i] = 1 - eps;
            } else {
                vi = (v[i] <= eps) ? eps : ((v[i] >= 1 - eps) ? 1 - eps : v[i]);
                mlogxi = -log(x[i]);
                mlogvi = -log(vi);
                tmp = R_pow(mlogxi, theta) + R_pow(mlogvi, theta);
                ui = exp(-R_pow(tmp, 1/theta)) * R_pow(mlogvi, theta-1) *
                        R_pow(tmp, 1/theta-1) / vi;
                u[i] = (ui <= eps) ? eps : ((ui >= 1 - eps) ? 1 - eps : ui);
            }
        }

        UNPROTECT(1);
    }

    return U;
}
