/* vines: R package for multivariate dependence modeling with vines
 * Copyright (C) 2010-2011 Yasser González-Fernández
 * Copyright (C) 2010-2011 Marta Soto
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "itervine.h"
#include "common.h"


/* Based on the Algorithm 3 from Aas, K., Czado, C., Frigessi, A. & Bakken, H.
 * Pair-copula constructions of multiple dependence. Insurance Mathematics
 * and Economics, 2009, Vol. 44, pp. 182-198.
 */

SEXP iterCVine(SEXP Vine, SEXP Data, SEXP Fit, SEXP Eval) {
    int trees, d, n, nEvals, currentEval;
    SEXP Copulas, Evals, V, X, Y;
    SEXP R, Rnames;

    PROTECT(Vine = duplicate(Vine));
    trees = asInteger(R_do_slot(Vine, mkString("trees")));
    d = asInteger(R_do_slot(Vine, mkString("dimension")));
    Copulas = R_do_slot(Vine, mkString("copulas"));
    n = nrows(Data);
    nEvals = isNull(Eval) ? 0 : (d-1)*d/2 - (d-trees-1)*(d-trees)/2;
    currentEval = 0;

    PROTECT(R = allocVector(VECSXP, 2));
    PROTECT(Rnames = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(R, 0, Vine);
    SET_VECTOR_ELT(Rnames, 0, mkString("vine"));
    PROTECT(Evals = allocVector(VECSXP, nEvals));
    SET_VECTOR_ELT(R, 1, Evals);
    SET_VECTOR_ELT(Rnames, 1, mkString("evals"));
    namesgets(R, Rnames);

    if (trees == 0) { /* Independence Canonical vine. */
        UNPROTECT(4);
        return R;
    }

    /* The indexes of the second dimension of the v array differs with
     * the indexes of the first dimension of the v array in Algorithm 3
     * because of the 1-based indexing. */

    PROTECT(V = alloc3DArray(REALSXP, n, d-1, d));
    PROTECT(X = allocVector(REALSXP, n));
    PROTECT(Y = allocVector(REALSXP, n));
    for (int i = 1; i <= d; i++) {
        for (int k = 1; k <= n; k++) {
            SET_REAL_3D(V, k, 1, i, n, d-1, GET_REAL_2D(Data, k, i, n));
        }
    }
    for (int j = 1; j <= trees; j++) {
        R_CheckUserInterrupt(); /* Allowing interrupts. */

        for (int i = 1; i <= d-j; i++) {
            for (int k = 1; k <= n; k++) {
                SET_REAL_1D(X, k, GET_REAL_3D(V, k, j, 1, n, d-1));
                SET_REAL_1D(Y, k, GET_REAL_3D(V, k, j, i+1, n, d-1));
            }
            if (!isNull(Fit)) {
                SET_VECTOR_2D(Copulas, j, i, d-1,
                        evalFunction5(Fit, Vine,
                                ScalarInteger(j), ScalarInteger(i), X, Y));
            }
            if (!isNull(Eval)) {
                SET_VECTOR_ELT(Evals, currentEval++,
                        evalFunction5(Eval, Vine,
                                ScalarInteger(j), ScalarInteger(i), X, Y));
            }
        }
        if (j == trees) {
            break;
        }
        for (int i = 1; i <= d-j; i++) { /* Observations for the next tree. */
            for (int k = 1; k <= n; k++) {
                SET_REAL_3D(V, k, j+1, i, n, d-1,
                        asReal(h(GET_VECTOR_2D(Copulas, j, i, d-1),
                               ScalarReal(GET_REAL_3D(V, k , j, i+1, n, d-1)),
                               ScalarReal(GET_REAL_3D(V, k , j, 1, n, d-1)))));
            }
        }
    }

    UNPROTECT(7);

    return R;
}


/* Based on the Algorithm 4 from Aas, K., Czado, C., Frigessi, A. & Bakken, H.
 * Pair-copula constructions of multiple dependence. Insurance Mathematics
 * and Economics, 2009, Vol. 44, pp. 182-198.
 */

SEXP iterDVine(SEXP Vine, SEXP Data, SEXP Fit, SEXP Eval) {
    int trees, d, n, nEvals, currentEval;
    SEXP Copulas, Evals, V, X, Y;
    SEXP R, Rnames;

    PROTECT(Vine = duplicate(Vine));
    trees = asInteger(R_do_slot(Vine, mkString("trees")));
    d = asInteger(R_do_slot(Vine, mkString("dimension")));
    Copulas = R_do_slot(Vine, mkString("copulas"));
    n = nrows(Data);
    nEvals = isNull(Eval) ? 0 : (d-1)*d/2 - (d-trees-1)*(d-trees)/2;
    currentEval = 0;

    PROTECT(R = allocVector(VECSXP, 2));
    PROTECT(Rnames = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(R, 0, Vine);
    SET_VECTOR_ELT(Rnames, 0, mkString("vine"));
    PROTECT(Evals = allocVector(VECSXP, nEvals));
    SET_VECTOR_ELT(R, 1, Evals);
    SET_VECTOR_ELT(Rnames, 1, mkString("evals"));
    namesgets(R, Rnames);

    if (trees == 0) { /* Independence D-vine. */
        UNPROTECT(4);
        return R;
    }

    /* The indexes of the second dimension of the v array differs with
     * the indexes of the first dimension of the v array in Algorithm 4
     * because of the 1-based indexing. */

    PROTECT(V = alloc3DArray(REALSXP, n, d, imax2(2*d-4, d)));
    PROTECT(X = allocVector(REALSXP, n));
    PROTECT(Y = allocVector(REALSXP, n));
    for (int i = 1; i <= d; i++) {
        for (int l = 1; l <= n; l++) {
            SET_REAL_3D(V, l, 1, i, n, d-1, GET_REAL_2D(Data, l, i, n));
        }
    }
    for (int i = 1; i <= d-1; i++) {
        for (int l = 1; l <= n; l++) {
            SET_REAL_1D(X, l, GET_REAL_3D(V, l, 1, i, n, d));
            SET_REAL_1D(Y, l, GET_REAL_3D(V, l, 1, i+1, n, d));
        }
        if (!isNull(Fit)) {
            SET_VECTOR_2D(Copulas, 1, i, d-1,
                    evalFunction5(Fit, Vine,
                            ScalarInteger(1), ScalarInteger(i), X, Y));
        }
        if (!isNull(Eval)) {
            SET_VECTOR_ELT(Evals, currentEval++,
                    evalFunction5(Eval, Vine,
                            ScalarInteger(1), ScalarInteger(i), X, Y));
        }
    }
    for (int l = 1; l <= n; l++) {
        SET_REAL_3D(V, l, 2, 1, n, d,
                asReal(h(GET_VECTOR_2D(Copulas, 1, 1, d-1),
                       ScalarReal(GET_REAL_3D(V, l, 1, 1, n, d)),
                       ScalarReal(GET_REAL_3D(V, l, 1, 2, n, d)))));
    }
    for (int k = 1; k <= imax2(d-3, 0); k++) {
        for (int l = 1; l <= n; l++) {
            SET_REAL_3D(V, l, 2, 2*k, n, d,
                    asReal(h(GET_VECTOR_2D(Copulas, 1, k+1, d-1),
                           ScalarReal(GET_REAL_3D(V, l, 1, k+2, n, d)),
                           ScalarReal(GET_REAL_3D(V, l, 1, k+1, n, d)))));
        }
        for (int l = 1; l <= n; l++) {
            SET_REAL_3D(V, l, 2, 2*k+1, n, d,
                    asReal(h(GET_VECTOR_2D(Copulas, 1, k+1, d-1),
                           ScalarReal(GET_REAL_3D(V, l, 1, k+1, n, d)),
                           ScalarReal(GET_REAL_3D(V, l, 1, k+2, n, d)))));
        }
    }
    for (int l = 1; l <= n; l++) {
        SET_REAL_3D(V, l, 2, 2*d-4, n, d,
                asReal(h(GET_VECTOR_2D(Copulas, 1, d-1, d-1),
                       ScalarReal(GET_REAL_3D(V, l, 1, d, n, d)),
                       ScalarReal(GET_REAL_3D(V, l, 1, d-1, n, d)))));
    }
    for (int j = 2; j <= trees; j++) {
        R_CheckUserInterrupt(); /* Allowing interrupts. */

        for (int i = 1; i <= d-j; i++) {
            for (int l = 1; l <= n; l++) {
                SET_REAL_1D(X, l, GET_REAL_3D(V, l, j, 2*i-1, n, d));
                SET_REAL_1D(Y, l, GET_REAL_3D(V, l, j, 2*i, n, d));
            }
            if (!isNull(Fit)) {
                SET_VECTOR_2D(Copulas, j, i, d-1,
                        evalFunction5(Fit, Vine,
                                ScalarInteger(j), ScalarInteger(i), X, Y));
            }
            if (!isNull(Eval)) {
                SET_VECTOR_ELT(Evals, currentEval++,
                        evalFunction5(Eval, Vine,
                                ScalarInteger(j), ScalarInteger(i), X, Y));
            }
        }
        if (j == trees) {
            break;
        }
        /* Observations for the next tree. */
        for (int l = 1; l <= n; l++) {
            SET_REAL_3D(V, l, j+1, 1, n, d,
                    asReal(h(GET_VECTOR_2D(Copulas, j, 1, d-1),
                           ScalarReal(GET_REAL_3D(V, l, j, 1, n, d)),
                           ScalarReal(GET_REAL_3D(V, l, j, 2, n, d)))));
        }
        if (d > 4) {
            for (int i = 1; i <= d-j-2; i++) {
                for (int l = 1; l <= n; l++) {
                    SET_REAL_3D(V, l, j+1, 2*i, n, d,
                            asReal(h(GET_VECTOR_2D(Copulas, j, i+1, d-1),
                                   ScalarReal(GET_REAL_3D(V, l, j, 2*i+2, n, d)),
                                   ScalarReal(GET_REAL_3D(V, l, j, 2*i+1, n, d)))));
                }
                for (int l = 1; l <= n; l++) {
                    SET_REAL_3D(V, l, j+1, 2*i+1, n, d,
                            asReal(h(GET_VECTOR_2D(Copulas, j, i+1, d-1),
                                   ScalarReal(GET_REAL_3D(V, l, j, 2*i+1, n, d)),
                                   ScalarReal(GET_REAL_3D(V, l, j, 2*i+2, n, d)))));
                }
            }
        }
        for (int l = 1; l <= n; l++) {
            SET_REAL_3D(V, l, j+1, 2*d-2*j-2, n, d,
                    asReal(h(GET_VECTOR_2D(Copulas, j, d-j, d-1),
                           ScalarReal(GET_REAL_3D(V, l, j, 2*d-2*j, n, d)),
                           ScalarReal(GET_REAL_3D(V, l, j, 2*d-2*j-1, n, d)))));
        }
    }

    UNPROTECT(7);

    return R;
}
