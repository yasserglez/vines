/* vines: Multivariate Dependence Modeling with Vines
 * Copyright (C) 2010, 2011 Yasser González-Fernández <ygf@icimaf.cu>
 * Copyright (C) 2010, 2011 Marta Soto <mrosa@icimaf.cu>
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

SEXP hinverseCopula(SEXP Copula, SEXP U, SEXP V);
SEXP hinverseNormalCopula(SEXP Rho, SEXP U, SEXP V);
SEXP hinverseIndepCopula(SEXP U, SEXP V);
SEXP hinverseTCopula(SEXP Rho, SEXP Df, SEXP U, SEXP V);
SEXP hinverseClaytonCopula(SEXP THETA, SEXP U, SEXP V);
SEXP hinverseFrankCopula(SEXP Theta, SEXP U, SEXP V);
