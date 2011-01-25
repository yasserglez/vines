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


/* Based on the Algorithm 1 from Aas, K., Czado, C., Frigessi, A. & Bakken, H.
 * Pair-copula constructions of multiple dependence. Insurance Mathematics
 * and Economics, 2009, Vol. 44, pp. 182-198.
 */

SEXP rCVine(SEXP Vine, SEXP N) {
    int n = asInteger(N);
    int d = asInteger(R_do_slot(Vine, mkString("dimension")));
    int trees = asInteger(R_do_slot(Vine, mkString("trees")));
    SEXP Copulas = R_do_slot(Vine, mkString("copulas"));
    SEXP W, V, R;

    PROTECT(W = allocMatrix(REALSXP, n, d));
    GetRNGstate();
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= d; j++) {
            SET_REAL_2D(W, i, j, n, runif(0, 1));
        }
    }
    PutRNGstate();

    if (trees == 0) {
        UNPROTECT(1);
        return W;
    }

    PROTECT(V = allocMatrix(REALSXP, d, d));
    PROTECT(R = allocMatrix(REALSXP, n, d));

    for (int i = 1; i <= n; i++) { /* Values of the first variable. */
        SET_REAL_2D(R, i, 1, n, GET_REAL_2D(W, i, 1, n));
    }

    for (int s = 1; s <= n; s++) { /* Loop over samples. */
        SET_REAL_2D(V, 1, 1, d, GET_REAL_2D(R, s, 1, n));

        for (int i = 2; i <= d; i++) { /* Loop over the other variables. */
            SET_REAL_2D(V, i, 1, d, GET_REAL_2D(W, s, i, n));
            for (int k = imin2(trees, i - 1); k >= 1; k--) {
                SET_REAL_2D(V, i, 1, d,
                        asReal(hinverse(GET_VECTOR_2D(Copulas, k, i - k, d - 1),
                               ScalarReal(GET_REAL_2D(V, i, 1, d)),
                               ScalarReal(GET_REAL_2D(V, k, k, d)))));
            }
            SET_REAL_2D(R, s, i, n, GET_REAL_2D(V, i, 1, d));

            if (i == d) {
                break;
            }

            for (int j = 1; j <= imin2(trees, i - 1); j++) {
                SET_REAL_2D(V, i, j + 1, d,
                        asReal(h(GET_VECTOR_2D(Copulas, j, i - j, d - 1),
                               ScalarReal(GET_REAL_2D(V, i, j, d)),
                               ScalarReal(GET_REAL_2D(V, j, j, d)))));
            }
        }

        R_CheckUserInterrupt(); /* Allowing interrupts. */
    }

    UNPROTECT(3);

    return R;
}
