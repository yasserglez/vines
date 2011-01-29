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
#include <R_ext/Utils.h>

#include "common.h"
#include "rvine.h"


/* Algorithm 1 from Aas, K., Czado, C., Frigessi, A. & Bakken, H.
 * Pair-copula constructions of multiple dependence. Insurance
 * Mathematics and Economics, 2009, Vol. 44, pp. 182-198.
 */

SEXP rCVine(SEXP vine, SEXP samples) {
    int n, d, trees;
    SEXP copulas, w, v;
    SEXP ans;

    d = asInteger(R_do_slot(vine, mkString("dimension")));
    trees = asInteger(R_do_slot(vine, mkString("trees")));
    copulas = R_do_slot(vine, mkString("copulas"));
    n = asInteger(samples);

    PROTECT(w = allocMatrix(REALSXP, n, d));
    GetRNGstate();
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= d; j++) {
            SET_REAL_2D(w, i, j, n, runif(0, 1));
        }
    }
    PutRNGstate();

    if (trees == 0) { /* Independence Canonical vine. */
        UNPROTECT(1);
        return w;
    }

    PROTECT(v = allocMatrix(REALSXP, d, d));
    PROTECT(ans = allocMatrix(REALSXP, n, d));

    for (int i = 1; i <= n; i++) { /* Values of the first variable. */
        SET_REAL_2D(ans, i, 1, n, GET_REAL_2D(w, i, 1, n));
    }

    for (int s = 1; s <= n; s++) { /* Loop over samples. */
        R_CheckUserInterrupt(); /* Allowing to user to stop the execution. */

        SET_REAL_2D(v, 1, 1, d, GET_REAL_2D(ans, s, 1, n));

        for (int i = 2; i <= d; i++) { /* Loop over the other variables. */
            SET_REAL_2D(v, i, 1, d, GET_REAL_2D(w, s, i, n));
            for (int k = imin2(trees, i-1); k >= 1; k--) {
                SET_REAL_2D(v, i, 1, d,
                        asReal(hinverse(GET_VECTOR_2D(copulas, k, i-k, d-1),
                                        ScalarReal(GET_REAL_2D(v, i, 1, d)),
                                        ScalarReal(GET_REAL_2D(v, k, k, d)))));
            }
            SET_REAL_2D(ans, s, i, n, GET_REAL_2D(v, i, 1, d));

            if (i == d) {
                break;
            }

            for (int j = 1; j <= imin2(trees, i-1); j++) {
                SET_REAL_2D(v, i, j+1, d,
                        asReal(h(GET_VECTOR_2D(copulas, j, i-j, d-1),
                                 ScalarReal(GET_REAL_2D(v, i, j, d)),
                                 ScalarReal(GET_REAL_2D(v, j, j, d)))));
            }
        }
    }

    UNPROTECT(3);

    return ans;
}


/* Algorithm 2 from Aas, K., Czado, C., Frigessi, A. & Bakken, H.
 * Pair-copula constructions of multiple dependence. Insurance
 * Mathematics and Economics, 2009, Vol. 44, pp. 182-198.
 */

SEXP rDVine(SEXP vine, SEXP samples) {
    int n, d, trees;
    SEXP copulas, w, v;
    SEXP ans;

    d = asInteger(R_do_slot(vine, mkString("dimension")));
    trees = asInteger(R_do_slot(vine, mkString("trees")));
    copulas = R_do_slot(vine, mkString("copulas"));
    n = asInteger(samples);

    PROTECT(w = allocMatrix(REALSXP, n, d));
    GetRNGstate();
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= d; j++) {
            SET_REAL_2D(w, i, j, n, runif(0, 1));
        }
    }
    PutRNGstate();

    PROTECT(v = allocMatrix(REALSXP, d, imax2(2*d-4, d)));
    PROTECT(ans = allocMatrix(REALSXP, n, d));

    for (int i = 1; i <= n; i++) { /* Values of the first variable. */
        SET_REAL_2D(ans, i, 1, n, GET_REAL_2D(w, i, 1, n));
    }
    for (int i = 1; i <= n; i++) { /* Values of the second variable. */
        SET_REAL_2D(ans, i, 2, n,
                asReal(hinverse(GET_VECTOR_2D(copulas, 1, 1, d-1),
                                ScalarReal(GET_REAL_2D(w, i, 2, n)),
                                ScalarReal(GET_REAL_2D(w, i, 1, n)))));
    }

    if (d == 2) { /* Stop if there are only two variables. */
        UNPROTECT(3);
        return ans;
    }

    for (int s = 1; s <= n; s++) { /* Loop over samples. */
        R_CheckUserInterrupt(); /* Allowing to user to stop the execution. */

        SET_REAL_2D(v, 1, 1, d, GET_REAL_2D(ans, s, 1, n));
        SET_REAL_2D(v, 2, 1, d, GET_REAL_2D(ans, s, 2, n));
        SET_REAL_2D(v, 2, 2, d,
                asReal(h(GET_VECTOR_2D(copulas, 1, 1, d-1),
                         ScalarReal(GET_REAL_2D(v, 1, 1, d)),
                         ScalarReal(GET_REAL_2D(v, 2, 1, d)))));

        for (int i = 3; i <= d; i++) { /* Loop over the other variables. */
            SET_REAL_2D(v, i, 1, d, GET_REAL_2D(w, s, i, n));

            if (trees >= 2) {
                for (int k = imin2(trees, i-1); k >= 2; k--) {
                    SET_REAL_2D(v, i, 1, d,
                            asReal(hinverse(GET_VECTOR_2D(copulas, k, i-k, d-1),
                                            ScalarReal(GET_REAL_2D(v, i, 1, d)),
                                            ScalarReal(GET_REAL_2D(v, i-1, 2*k-2, d)))));
                }
            }
            SET_REAL_2D(v, i, 1, d,
                    asReal(hinverse(GET_VECTOR_2D(copulas, 1, i-1, d-1),
                                    ScalarReal(GET_REAL_2D(v, i, 1, d)),
                                    ScalarReal(GET_REAL_2D(v, i-1, 1, d)))));
            SET_REAL_2D(ans, s, i, n, GET_REAL_2D(v, i, 1, d));

            if (i == d) {
                break;
            }

            if (trees >= 2) {
                SET_REAL_2D(v, i, 2, d,
                        asReal(h(GET_VECTOR_2D(copulas, 1, i-1, d-1),
                                 ScalarReal(GET_REAL_2D(v, i-1, 1, d)),
                                 ScalarReal(GET_REAL_2D(v, i, 1, d)))));
            }
            if (trees >= 3) {
                SET_REAL_2D(v, i, 3, d,
                        asReal(h(GET_VECTOR_2D(copulas, 1, i-1, d-1),
                                 ScalarReal(GET_REAL_2D(v, i, 1, d)),
                                 ScalarReal(GET_REAL_2D(v, i-1, 1, d)))));
            }
            if (trees >= 3 && i > 3) {
                for (int j = 2; j <= imin2(trees-1, i-2); j++) {
                    SET_REAL_2D(v, i, 2*j, d,
                            asReal(h(GET_VECTOR_2D(copulas, j, i-j, d-1),
                                     ScalarReal(GET_REAL_2D(v, i-1, 2*j-2, d)),
                                     ScalarReal(GET_REAL_2D(v, i, 2*j-1, d)))));
                    SET_REAL_2D(v, i, 2*j+1, d,
                            asReal(h(GET_VECTOR_2D(copulas, j, i-j, d-1),
                                     ScalarReal(GET_REAL_2D(v, i, 2*j-1, d)),
                                     ScalarReal(GET_REAL_2D(v, i-1, 2*j-2, d)))));
                }
            }
            if (trees >= i) {
                SET_REAL_2D(v, i, 2*i-2, d,
                        asReal(h(GET_VECTOR_2D(copulas, i-1, 1, d-1),
                                 ScalarReal(GET_REAL_2D(v, i-1, 2*i-4, d)),
                                 ScalarReal(GET_REAL_2D(v, i, 2*i-3, d)))));
            }
        }
    }

    UNPROTECT(3);

    return ans;
}
