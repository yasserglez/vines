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

#define REAL1(x,i) REAL(x)[(i)-1)]
#define REAL2(x,i,j,ni) REAL(x)[((i)-1)+(ni)*((j)-1)]
#define REAL3(x,i,j,k,ni,nj) REAL(x)[((i)-1)+(ni)*((j)-1)+(nj)*((k)-1)]

#define SET_REAL1(x,i,v) REAL1(x,i) = (v)
#define SET_REAL2(x,i,j,ni,v) REAL2(x,i,j,ni) = (v)
#define SET_REAL3(x,i,j,k,ni,nj,v) REAL3(x,i,j,k,ni,nj) = (v)

#define VECTOR1_ELT(x,i) VECTOR_ELT(x,(i)-1)
#define VECTOR2_ELT(x,i,j,ni) VECTOR_ELT(x,((i)-1)+(ni)*((j)-1))
#define VECTOR3_ELT(x,i,j,k,ni,nj) VECTOR_ELT(x,((i)-1)+(ni)*((j)-1)+(nj)*((k)-1))

#define SET_VECTOR1_ELT(x,i,v) SET_VECTOR_ELT(x,(i)-1,v)
#define SET_VECTOR2_ELT(x,i,j,ni,v) SET_VECTOR_ELT(x,((i)-1)+(ni)*((j)-1),v)
#define SET_VECTOR3_ELT(x,i,j,k,ni,nj,v) SET_VECTOR_ELT(x,((i)-1)+(ni)*((j)-1)+(nj)*((k)-1),v)

SEXP h(SEXP COPULA, SEXP X, SEXP V);
SEXP hinverse(SEXP COPULA, SEXP U, SEXP V);
