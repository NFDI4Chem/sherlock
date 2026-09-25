/*
    file: OutlsdSrc/matutil.c
    Copyright(C)2000 CNRS - UMR 7312 - Jean-Marc Nuzillard

    This file is part of LSD.

    LSD is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    LSD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with LSD; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 * These routines are simplified versions of those
 * provided in "Numerical Recipes in C" for dynamic
 * objects allocation. Indexes start at 0.
 */

#include "matutil.h"

#include <stdlib.h>

double *dvector(unsigned int sz)
{
	return (double *) malloc((size_t) (sz * sizeof(double))) ;
}

void free_dvector(double *v)
{
	free((void *) (v)) ;
}

double **dmatrix(unsigned int sz1, unsigned int sz2)
{
	unsigned int i;
	double **m;

	m = (double **) malloc((size_t)(sz1 * sizeof(double*)));
	if(!m) {
		return m;
	}

	m[0] = (double *) malloc((size_t)(sz1 * sz2 * sizeof(double)));
	if(!m[0]) {
		free((void *) (m));
		return (double **)m[0];
	}
	for(i=1 ; i<sz1 ; i++) m[i] = m[i-1] + sz2;

	return m;
}

void free_dmatrix(double **m)
{
	free((void *) (m[0]));
	free((void *) (m));
}

