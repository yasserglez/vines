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


SEXP h(SEXP Copula, SEXP X, SEXP V) {
    SEXP F, Call, U;

    PROTECT(F = findFun(install("h"), R_GlobalEnv));
    PROTECT(Call = lang4(F, Copula, X, V));
    PROTECT(U = coerceVector(eval(Call, R_GlobalEnv), REALSXP));
    UNPROTECT(3);

    return U;
}

SEXP hinverse(SEXP Copula, SEXP U, SEXP V) {
    SEXP F, Call, X;

    PROTECT(F = findFun(install("hinverse"), R_GlobalEnv));
    PROTECT(Call = lang4(F, Copula, U, V));
    PROTECT(X = coerceVector(eval(Call, R_GlobalEnv), REALSXP));
    UNPROTECT(3);

    return X;
}

SEXP evalFunction5(SEXP F, SEXP X1, SEXP X2, SEXP X3, SEXP X4, SEXP X5) {
    SEXP Call, T, R;

    PROTECT(Call = T = allocList(6));
    SET_TYPEOF(Call, LANGSXP);
    SETCAR(T, F); T = CDR(T);
    SETCAR(T, X1); T = CDR(T);
    SETCAR(T, X2); T = CDR(T);
    SETCAR(T, X3); T = CDR(T);
    SETCAR(T, X4); T = CDR(T);
    SETCAR(T, X5);
    R = eval(Call, R_GlobalEnv);
    UNPROTECT(1);

    return R;
}
