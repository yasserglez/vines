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

#include "common.h"
#include "itervine.h"


/* Algorithm 3 from Aas, K., Czado, C., Frigessi, A. & Bakken, H.
 * Pair-copula constructions of multiple dependence. Insurance
 * Mathematics and Economics, 2009, Vol. 44, pp. 182-198.
 */

SEXP iterCVine(SEXP vine, SEXP data, SEXP fFit, SEXP fEval) {
    int trees, d, n;
    int totalEvals, currentEval;
    SEXP copulas, evals, v, x, y;
    SEXP ans, ansNames;

    PROTECT(vine = duplicate(vine));
    d = asInteger(R_do_slot(vine, mkString("dimension")));
    trees = asInteger(R_do_slot(vine, mkString("trees")));
    copulas = R_do_slot(vine, mkString("copulas"));
    n = nrows(data);
    totalEvals = isNull(fEval) ? 0 : (d-1)*d/2 - (d-trees-1)*(d-trees)/2;
    currentEval = 0;

    PROTECT(ans = allocVector(VECSXP, 2));
    PROTECT(ansNames = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, vine);
    SET_VECTOR_ELT(ansNames, 0, mkString("vine"));
    PROTECT(evals = allocVector(VECSXP, totalEvals));
    SET_VECTOR_ELT(ans, 1, evals);
    SET_VECTOR_ELT(ansNames, 1, mkString("evals"));
    namesgets(ans, ansNames);

    if (trees == 0) { /* Independence Canonical vine. */
        UNPROTECT(4);
        return ans;
    }

    /* The indexes of the second dimension of the v array differs with
     * the ones from the first dimension of the v array in Algorithm 3
     * because we are numbering the indexes starting at 1.
     */

    PROTECT(v = alloc3DArray(REALSXP, n, d-1, d));
    PROTECT(x = allocVector(REALSXP, n));
    PROTECT(y = allocVector(REALSXP, n));

    for (int i = 1; i <= d; i++) {
        for (int k = 1; k <= n; k++) {
            SET_REAL_3D(v, k, 1, i, n, d-1, GET_REAL_2D(data, k, i, n));
        }
    }
    for (int j = 1; j <= trees; j++) {
        R_CheckUserInterrupt(); /* Allowing to user to stop the execution. */

        for (int i = 1; i <= d-j; i++) {
            for (int t = 1; t <= n; t++) {
                SET_REAL_1D(x, t, GET_REAL_3D(v, t, j, 1, n, d-1));
                SET_REAL_1D(y, t, GET_REAL_3D(v, t, j, i+1, n, d-1));
            }
            if (!isNull(fFit)) {
                SET_VECTOR_2D(copulas, j, i, d-1,
                        evalFunction5(fFit, vine,
                                      ScalarInteger(j), ScalarInteger(i),
                                      x, y));
            }
            if (!isNull(fEval)) {
                SET_VECTOR_ELT(evals, currentEval++,
                        evalFunction5(fEval, vine,
                                      ScalarInteger(j), ScalarInteger(i),
                                      x, y));
            }
        }
        if (j == trees) {
            break;
        }
        /* Observations for the next tree. */
        for (int i = 1; i <= d-j; i++) {
            for (int t = 1; t <= n; t++) {
                SET_REAL_3D(v, t, j+1, i, n, d-1,
                        asReal(h(GET_VECTOR_2D(copulas, j, i, d-1),
                                 ScalarReal(GET_REAL_3D(v, t , j, i+1, n, d-1)),
                                 ScalarReal(GET_REAL_3D(v, t , j, 1, n, d-1)))));
            }
        }
    }

    UNPROTECT(7);

    return ans;
}


/* Algorithm 4 from Aas, K., Czado, C., Frigessi, A. & Bakken, H.
 * Pair-copula constructions of multiple dependence. Insurance
 * Mathematics and Economics, 2009, Vol. 44, pp. 182-198.
 */

SEXP iterDVine(SEXP vine, SEXP data, SEXP fFit, SEXP fEval) {
    int d, n, trees;
    int totalEvals, currentEval;
    SEXP copulas, evals, v, x, y;
    SEXP ans, ansNames;

    PROTECT(vine = duplicate(vine));
    d = asInteger(R_do_slot(vine, mkString("dimension")));
    trees = asInteger(R_do_slot(vine, mkString("trees")));
    copulas = R_do_slot(vine, mkString("copulas"));
    n = nrows(data);
    totalEvals = isNull(fEval) ? 0 : (d-1)*d/2 - (d-trees-1)*(d-trees)/2;
    currentEval = 0;

    PROTECT(ans = allocVector(VECSXP, 2));
    PROTECT(ansNames = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, vine);
    SET_VECTOR_ELT(ansNames, 0, mkString("vine"));
    PROTECT(evals = allocVector(VECSXP, totalEvals));
    SET_VECTOR_ELT(ans, 1, evals);
    SET_VECTOR_ELT(ansNames, 1, mkString("evals"));
    namesgets(ans, ansNames);

    if (trees == 0) { /* Independence D-vine. */
        UNPROTECT(4);
        return ans;
    }

    /* The indexes of the second dimension of the v array differs with
     * the ones from the first dimension of the v array in Algorithm 4
     * because we are numbering the indexes starting at 1.
     */

    PROTECT(v = alloc3DArray(REALSXP, n, d, imax2(2*d-4, d)));
    PROTECT(x = allocVector(REALSXP, n));
    PROTECT(y = allocVector(REALSXP, n));

    for (int i = 1; i <= d; i++) {
        for (int t = 1; t <= n; t++) {
            SET_REAL_3D(v, t, 1, i, n, d, GET_REAL_2D(data, t, i, n));
        }
    }
    for (int i = 1; i <= d-1; i++) {
        for (int t = 1; t <= n; t++) {
            SET_REAL_1D(x, t, GET_REAL_3D(v, t, 1, i, n, d));
            SET_REAL_1D(y, t, GET_REAL_3D(v, t, 1, i+1, n, d));
        }
        if (!isNull(fFit)) {
            SET_VECTOR_2D(copulas, 1, i, d-1,
                    evalFunction5(fFit, vine,
                                  ScalarInteger(1), ScalarInteger(i),
                                  x, y));
        }
        if (!isNull(fEval)) {
            SET_VECTOR_ELT(evals, currentEval++,
                    evalFunction5(fEval, vine,
                                  ScalarInteger(1), ScalarInteger(i),
                                  x, y));
        }
    }
    for (int t = 1; t <= n; t++) {
        SET_REAL_3D(v, t, 2, 1, n, d,
                asReal(h(GET_VECTOR_2D(copulas, 1, 1, d-1),
                         ScalarReal(GET_REAL_3D(v, t, 1, 1, n, d)),
                         ScalarReal(GET_REAL_3D(v, t, 1, 2, n, d)))));
    }
    for (int k = 1; k <= d-3; k++) {
        for (int t = 1; t <= n; t++) {
            SET_REAL_3D(v, t, 2, 2*k, n, d,
                    asReal(h(GET_VECTOR_2D(copulas, 1, k+1, d-1),
                             ScalarReal(GET_REAL_3D(v, t, 1, k+2, n, d)),
                             ScalarReal(GET_REAL_3D(v, t, 1, k+1, n, d)))));
        }
        for (int t = 1; t <= n; t++) {
            SET_REAL_3D(v, t, 2, 2*k+1, n, d,
                    asReal(h(GET_VECTOR_2D(copulas, 1, k+1, d-1),
                             ScalarReal(GET_REAL_3D(v, t, 1, k+1, n, d)),
                             ScalarReal(GET_REAL_3D(v, t, 1, k+2, n, d)))));
        }
    }
    for (int t = 1; t <= n; t++) {
        SET_REAL_3D(v, t, 2, 2*d-4, n, d,
                asReal(h(GET_VECTOR_2D(copulas, 1, d-1, d-1),
                         ScalarReal(GET_REAL_3D(v, t, 1, d, n, d)),
                         ScalarReal(GET_REAL_3D(v, t, 1, d-1, n, d)))));
    }
    for (int j = 2; j <= trees; j++) {
        R_CheckUserInterrupt(); /* Allowing to user to stop the execution. */

        for (int i = 1; i <= d-j; i++) {
            for (int t = 1; t <= n; t++) {
                SET_REAL_1D(x, t, GET_REAL_3D(v, t, j, 2*i-1, n, d));
                SET_REAL_1D(y, t, GET_REAL_3D(v, t, j, 2*i, n, d));
            }
            if (!isNull(fFit)) {
                SET_VECTOR_2D(copulas, j, i, d-1,
                        evalFunction5(fFit, vine,
                                      ScalarInteger(j), ScalarInteger(i),
                                      x, y));
            }
            if (!isNull(fEval)) {
                SET_VECTOR_ELT(evals, currentEval++,
                        evalFunction5(fEval, vine,
                                      ScalarInteger(j), ScalarInteger(i),
                                      x, y));
            }
        }
        if (j == trees) {
            break;
        }
        /* Observations for the next tree. */
        for (int t = 1; t <= n; t++) {
            SET_REAL_3D(v, t, j+1, 1, n, d,
                    asReal(h(GET_VECTOR_2D(copulas, j, 1, d-1),
                             ScalarReal(GET_REAL_3D(v, t, j, 1, n, d)),
                             ScalarReal(GET_REAL_3D(v, t, j, 2, n, d)))));
        }
        if (d > 4) {
            for (int i = 1; i <= d-j-2; i++) {
                for (int t = 1; t <= n; t++) {
                    SET_REAL_3D(v, t, j+1, 2*i, n, d,
                            asReal(h(GET_VECTOR_2D(copulas, j, i+1, d-1),
                                     ScalarReal(GET_REAL_3D(v, t, j, 2*i+2, n, d)),
                                     ScalarReal(GET_REAL_3D(v, t, j, 2*i+1, n, d)))));
                }
                for (int t = 1; t <= n; t++) {
                    SET_REAL_3D(v, t, j+1, 2*i+1, n, d,
                            asReal(h(GET_VECTOR_2D(copulas, j, i+1, d-1),
                                     ScalarReal(GET_REAL_3D(v, t, j, 2*i+1, n, d)),
                                     ScalarReal(GET_REAL_3D(v, t, j, 2*i+2, n, d)))));
                }
            }
        }
        for (int t = 1; t <= n; t++) {
            SET_REAL_3D(v, t, j+1, 2*d-2*j-2, n, d,
                    asReal(h(GET_VECTOR_2D(copulas, j, d-j, d-1),
                             ScalarReal(GET_REAL_3D(v, t, j, 2*d-2*j, n, d)),
                             ScalarReal(GET_REAL_3D(v, t, j, 2*d-2*j-1, n, d)))));
        }
    }

    UNPROTECT(7);

    return ans;
}
