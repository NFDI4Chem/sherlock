/*
    file: LsdSrc/distance.c
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

#include "defs.h"
#include <stdio.h>
#define X 1000

extern struct 	atom at[MAX_ATOMES] ; /* lsd.c */
extern struct	pile heap[MAX_HEAP] ; /* lsd.c */
extern struct 	lcorr propcorr[MAX_CORR] ; /* lsd.c */
extern int	iheap ; /* heap.c */
extern int	max_atomes ; /* lsd.c */
extern int	maxj ; /* lsd.c */
extern int	ln[MAX_LN][MAX_LNBUF] ;
extern int	lnsize[MAX_LN] ;
extern int	flverb ;

extern void	wrln(int i) ; /* start.c */

int d1[MAX_ATOMES][MAX_ATOMES] ;
int d2[MAX_ATOMES][MAX_ATOMES] ;

int maxjok(void) ;
void distances(void) ;
int checkdist(int x, int y, int c) ;

int maxjok(void)
{
	int ok ;
	struct pilehist *a ;
	int i ;

	distances() ;
/* d1 is the distance matrix of the molecule */
	ok = TRUE ;
	for(i = 1 ; (i<=iheap) && ok ; i++) {
/* searches the heap for eliminated correlations */
		a = &heap[i].pinfo.phist ;
		if(heap[i].ptyp == 6) {
			ok = checkdist(a->h1, a->h2, a->corr) ;
/* checking for distance constraint */
		}
	}
	return ok ;
}

void distances(void)
{
	int i, j, k, l1, l2 ;

	for(i=1 ; i<=max_atomes ; i++)
		for(j=1 ; j<=max_atomes ; j++)
			d1[i][j] = X ;
/* all atoms are far away */
	for(i=1 ; i<=max_atomes ; i++)
		for(j=0 ; j<at[i].nlia ; j++)
			d1[i][at[i].lia[j]] = 1 ;
/* 1 is the distance between bonded atom */
	for(k=1 ; k<=max_atomes ; k++) {
		for(i=1 ; i<=max_atomes ; i++)
			for(j=1 ; j<=max_atomes ; j++) {
				l1 = d1[i][j] ;
				l2 = d1[i][k] + d1[k][j] ;
				d2[i][j] = (l1 < l2) ? l1 : l2 ;
			}
		for(i=1 ; i<=max_atomes ; i++)
			for(j=1 ; j<=max_atomes ; j++)
				d1[i][j] = d2[i][j] ;
	}
	for(i=1 ; i<=max_atomes ; i++)
		d1[i][i] = 0 ;
/* d1 is the distance matrix */
	if (flverb > 1) {
/* new with 3.1.2 */
		printf("matrix of atom distances\n") ;
		printf("    ") ;
		for (j = 1 ; j <= max_atomes ; j++) {
			printf("%3d", j) ;
		}
		printf("\n") ;
		for(i=1 ; i<=max_atomes ; i++) {
			printf("%3d:", i) ;
			for(j=1 ; j<=max_atomes ; j++)
				printf("%3d", d1[i][j]) ;
			printf("\n") ;
		}
	}
}

int checkdist(int x, int y, int c)
{
	int i ;
	int y0 ;
	int d, dmin, dmax ;

	dmin = propcorr[c].min ;
	dmax = propcorr[c].max ;
/* printf("x=%d y=%d c=%d dmin=%d dmax=%d\n", x, y, c, dmin, dmax) ; */
	if(y>0) {
/* correlation x-y is not ambiguous */
		d = d1[x][y] ;
		if(flverb > 1) {
			if(!dmax) {
				printf("distance(%d,%d) is %d and must be greater or equal to %d.\n", x, y, d, dmin) ;
			} else {
				printf("distance(%d,%d) is %d and must be in the [%d-%d] interval.\n", x, y, d, dmin, dmax) ;
			}
		}
		return d >= dmin && ((!dmax) ? TRUE : d <= dmax) ;
/* checks distance constraint */
	} else {
		y = -y ;
/* was buggy in 3.1.1 */
		if(flverb > 1) {
			printf ("distance between %d and one among ", x) ;
			wrln(y) ;
			if(!dmax) {
				printf("must be greater or equal to %d.\n", dmin) ;
			} else {
				printf("must be in the [%d-%d] interval.\n", dmin, dmax) ;
			}
		}
		for(i=0 ; i<lnsize[y] ; i++) {
/* correlation x-y is ambiguous */
			y0 = ln[y][i] ;
			d = d1[x][y0] ;
			if(d >= dmin && ((!dmax) ? TRUE : d <= dmax)) {
/* was buggy in 3.1.1 */
				if(flverb > 1) {
					printf("distance(%d,%d) is %d.\n", x, y0, d) ;
				}
				return TRUE ;
/* there is at least one alternative that validates the constraint */
			}
		}
		return FALSE ;
/* all distances are too high */
	}
}

