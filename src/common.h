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

#define GET_REAL_1D(x,i) REAL(x)[(i)-1)]
#define GET_REAL_2D(x,i,j,ni) REAL(x)[((i)-1)+(ni)*((j)-1)]
#define GET_REAL_3D(x,i,j,k,ni,nj) REAL(x)[((i)-1)+(ni)*((j)-1)+(nj)*((k)-1)]

#define SET_REAL_1D(x,i,v) GET_REAL_1D(x,i) = (v)
#define SET_REAL_2D(x,i,j,ni,v) GET_REAL_2D(x,i,j,ni) = (v)
#define SET_REAL_3D(x,i,j,k,ni,nj,v) GET_REAL_3D(x,i,j,k,ni,nj) = (v)

#define GET_VECTOR_1D(x,i) VECTOR_ELT(x,(i)-1)
#define GET_VECTOR_2D(x,i,j,ni) VECTOR_ELT(x,((i)-1)+(ni)*((j)-1))
#define GET_VECTOR_3D(x,i,j,k,ni,nj) VECTOR_ELT(x,((i)-1)+(ni)*((j)-1)+(nj)*((k)-1))

#define SET_VECTOR_1D(x,i,v) SET_VECTOR_ELT(x,(i)-1,v)
#define SET_VECTOR_2D(x,i,j,ni,v) SET_VECTOR_ELT(x,((i)-1)+(ni)*((j)-1),v)
#define SET_VECTOR_3D(x,i,j,k,ni,nj,v) SET_VECTOR_ELT(x,((i)-1)+(ni)*((j)-1)+(nj)*((k)-1),v)


SEXP h(SEXP Copula, SEXP X, SEXP V);
SEXP hinverse(SEXP Copula, SEXP U, SEXP V);
