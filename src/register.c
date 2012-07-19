/* vines: Multivariate Dependence Modeling with Vines
 * Copyright (C) 2010-2012 Yasser González Fernández <ygonzalezfernandez@gmail.com>
 * Copyright (C) 2010-2012 Marta Rosa Soto Ortiz <mrosa@icimaf.cu>
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
#include <R_ext/Rdynload.h>

#include "h.h"
#include "hinverse.h"


R_CallMethodDef callMethods[] = {
    { "hIndepCopula", (DL_FUNC) &hIndepCopula, 2 },
    { "hNormalCopula", (DL_FUNC) &hNormalCopula, 3 },
    { "hinverseNormalCopula", (DL_FUNC) &hinverseNormalCopula, 3 },
    { "hinverseIndepCopula", (DL_FUNC) &hinverseIndepCopula, 2 },
    { "hTCopula", (DL_FUNC) &hTCopula, 4 },
    { "hinverseTCopula", (DL_FUNC) &hinverseTCopula, 4 },
    { "hClaytonCopula", (DL_FUNC) &hClaytonCopula, 3 },
    { "hinverseClaytonCopula", (DL_FUNC) &hinverseClaytonCopula, 3 },
    { "hGumbelCopula", (DL_FUNC) &hGumbelCopula, 3 },
    { "hFGMCopula", (DL_FUNC) &hFGMCopula, 3 },
    { "hGalambosCopula", (DL_FUNC) &hGalambosCopula, 3 },
    { "hFrankCopula", (DL_FUNC) &hFrankCopula, 3 },
    { "hinverseFrankCopula", (DL_FUNC) &hinverseFrankCopula, 3 },
    { NULL, NULL, 0 }
};

void R_init_vines(DllInfo *info) {
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
