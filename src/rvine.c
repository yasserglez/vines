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

#include "common.h"
#include "rvine.h"


// Based on the Algorithm 1 from Aas, K., Czado, C., Frigessi, A. & Bakken, H.
// Pair-copula constructions of multiple dependence. Insurance Mathematics
// and Economics, 2009, Vol. 44, pp. 182-198.

SEXP rCVine(SEXP VINE, SEXP N) {
    int n = asInteger(N);
    int trees = asInteger(R_do_slot(VINE, mkString("trees")));
    int d = asInteger(R_do_slot(VINE, mkString("dimension")));
    SEXP COPULAS = R_do_slot(VINE, mkString("copulas"));
    SEXP W, V, R;

    PROTECT(W = allocMatrix(REALSXP, n, d));
    GetRNGstate();
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= d; j++) {
            SET_REAL2(W, i, j, n, runif(0, 1));
        }
    }
    PutRNGstate();

    if (trees == 0) {
        UNPROTECT(1);
        return W;
    }

    PROTECT(V = allocMatrix(REALSXP, d, d));
    PROTECT(R = allocMatrix(REALSXP, n, d));

    for (int i = 1; i <= n; i++) { // Values of the first variable.
        SET_REAL2(R, i, 1, n, REAL2(W, i, 1, n));
    }

    for (int s = 1; s <= n; s++) { // Loop over samples.
        SET_REAL2(V, 1, 1, d, REAL2(R, s, 1, n));

        for (int i = 2; i <= d; i++) { // Loop over the other variables.
            SET_REAL2(V, i, 1, d, REAL2(W, s, i, n));
            for (int k = imin2(trees, i - 1); k >= 1; k--) {
                SET_REAL2(V, i, 1, d,
                        asReal(hinverse(VECTOR2_ELT(COPULAS, k, i - k, d - 1),
                                ScalarReal(REAL2(V, i, 1, d)),
                                ScalarReal(REAL2(V, k, k, d)))));
            }
            SET_REAL2(R, s, i, n, REAL2(V, i, 1, d));

            if (i == d) {
                break;
            }

            for (int j = 1; j <= imin2(trees, i - 1); j++) {
                SET_REAL2(V, i, j + 1, d,
                        asReal(h(VECTOR2_ELT(COPULAS, j, i - j, d - 1),
                                ScalarReal(REAL2(V, i, j, d)),
                                ScalarReal(REAL2(V, j, j, d)))));
            }
        }
    }

    UNPROTECT(3);

    return R;
}
