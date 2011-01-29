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


SEXP h(SEXP copula, SEXP x, SEXP v) {
    SEXP f, fCall, ans;
    PROTECT_INDEX anspi;

    PROTECT(f = findFun(install("h"), R_GlobalEnv));
    PROTECT(fCall = lang4(f, copula, x, v));
    PROTECT_WITH_INDEX(ans = eval(fCall, R_GlobalEnv), &anspi);
    REPROTECT(ans = coerceVector(ans, REALSXP), anspi);

    UNPROTECT(3);

    return ans;
}

SEXP hinverse(SEXP copula, SEXP u, SEXP v) {
    SEXP f, fCall, ans;
    PROTECT_INDEX anspi;

    PROTECT(f = findFun(install("hinverse"), R_GlobalEnv));
    PROTECT(fCall = lang4(f, copula, u, v));
    PROTECT_WITH_INDEX(ans = eval(fCall, R_GlobalEnv), &anspi);
    REPROTECT(ans = coerceVector(ans, REALSXP), anspi);

    UNPROTECT(3);

    return ans;
}

SEXP evalFunction5(SEXP f, SEXP x1, SEXP x2, SEXP x3, SEXP x4, SEXP x5) {
    SEXP fCall, t, ans;

    PROTECT(fCall = t = allocList(6));
    SET_TYPEOF(fCall, LANGSXP);
    SETCAR(t, f); t = CDR(t);
    SETCAR(t, x1); t = CDR(t);
    SETCAR(t, x2); t = CDR(t);
    SETCAR(t, x3); t = CDR(t);
    SETCAR(t, x4); t = CDR(t);
    SETCAR(t, x5);
    PROTECT(ans = eval(fCall, R_GlobalEnv));

    UNPROTECT(2);

    return ans;
}
