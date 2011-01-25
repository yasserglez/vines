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


SEXP h(SEXP COPULA, SEXP X, SEXP V) {
    SEXP H, H_CALL, R;

    PROTECT(H = findFun(install("h"), R_GlobalEnv));
    PROTECT(H_CALL = lang4(H, COPULA, X, V));
    PROTECT(R = coerceVector(eval(H_CALL, R_GlobalEnv), REALSXP));

    UNPROTECT(3);

    return R;
}

SEXP hinverse(SEXP COPULA, SEXP U, SEXP V) {
    SEXP HINVERSE, HINVERSE_CALL, R;

    PROTECT(HINVERSE = findFun(install("hinverse"), R_GlobalEnv));
    PROTECT(HINVERSE_CALL = lang4(HINVERSE, COPULA, U, V));
    PROTECT(R = coerceVector(eval(HINVERSE_CALL, R_GlobalEnv), REALSXP));

    UNPROTECT(3);

    return R;
}
