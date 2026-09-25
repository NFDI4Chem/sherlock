/*
    file: LsdSrc/remcorr.c
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

extern struct atom at[MAX_ATOMES] ; /* lsd.c */
extern struct 	lcorr propcorr[MAX_CORR] ; /* lsd.c */
extern int	valcorr[MAX_CORR] ; /* lsd.c */
extern int	base[MAX_ATOMES] ; /* lsd.c */
extern int	nbase ; /* lsd.c */
extern int	flverb ; /* lsd.c */
extern int	flcorrvar ; /* lsd.c */
extern int	lnsize[MAX_LN] ; /* lsd.c */
extern int	ln[MAX_LN][MAX_LNBUF] ; /* lsd.c */
extern int	nivheap ; /* lsd.c */

extern void	empile0(int *a, int v) ; /* heap.c */
extern void	empile(int t, int a, int v) ; /* heap.c */
extern void	wrln2(int y) ; /* start1.c */

void	remcorr3(int fl, int n, int x, int y) ;
int		remcorr2(int x, int z, int d) ;
void	remcorrvar(int x, int z, int d) ;
int		inutil1(int x, int y) ;
int		inutil2(int x, int z, int y) ;
int		inutil4(int z, int x, int y) ;

void	remcorr3(int fl, int n, int x, int y)	
/*
 * removes the nth correlation, which is between x and y
 * if fl is true, it means that the correlation has to be removed because
 * it becomes obsolete. This happens upon bond formation when
 * x and y become visible to each other (either directly bound, or both 
 * bound to a common atom. The correlation is then considered as a 
 * confirmation information. fl is false when the correlation has been
 * really used in order to establish bonds.
 */
{
	int	*a ;

	empile0(&(valcorr[n]), FALSE) ;	
/* the correlation is marked as unusable */
	a = &at[x].nc ;
/* address where x's number of valid correlations is stored */
	empile0(a, (*a) - 1) ;
/* one less valid correlation for x */
	if (base[x] && !(*a)) {
/* if x belongs to the base and has no more correlation */
		empile0(&(base[x]), FALSE) ;
		empile0(&nbase , nbase - 1) ;
/* it is removed from the base */
	}
	if (y > 0) {
/* the correlation has no variant : y is not a numeric list reference.
 * y'th number of correlations and the base content is modified as for x */
		a = &at[y].nc ;
		empile0(a, (*a) - 1) ;
		if (base[y] && !(*a)) {
			empile0(&(base[y]), FALSE) ;
			empile0(&nbase , nbase - 1) ;
		}
	}
	if (fl) {
/* for confirmation correlations */
		empile(3, x, y) ;
		if (flverb > 1) {
/* message for confirmation correlations */
			printf("%d : removes correlation %d - ", nivheap, x) ;
			wrln2(y) ;		
		}
	}
}


int		remcorr2(int x, int z, int d) 
/*
 * removes the x-z correlation, if any (not known yet).
 * x and z are atom numbers.
 * when d equals 1, returns FALSE if there is a correlation
 * between x and z which is not j1ok.
 * when d equals 2, returns FALSE if there is a correlation
 * between x and z which is neither j1ok nor j2ok (jnok only).
 */
{
	int	i, c ;
	struct lcorr *pc ;
	
	if ((!at[x].ingrp) && (!at[x].nc)) {
		return TRUE ;
	}
	if ((!at[z].ingrp) && (!at[z].nc)) {
		return TRUE ;
	}
/* if x is not a member of a group and has no correlation, there is no
 * correlation to remove. As well for z. */
	for (i = 0 ; i < at[x].nctot ; i++) {
		c = at[x].ex[i] ;
		pc = &propcorr[c] ;
		if (valcorr[c] && (z == at[x].other[i])) {
			if(((d == 1) && pc->j1ok) || ((d == 2) && pc->j2ok)) {
/* the correlation is really a confirmation one because its min/max coupling path lengths do not conflict... */			
				remcorr3(TRUE, c, x, z) ;
/* ...so it is removed */
			} else if( (d == 2) && pc->j1ok ) {
/* the correlation does not invalidate the last hypothesis because cyclopropane is allowed... */
/* ...and it must not be removed */
			} else {
/* the correlation invalidates the last hypothesis because of its min/max coupling path lengths... */
				return FALSE ;
/* ...so backtrack */
			}
/* correlation c between x and z must be removed */
		}
	}
	if (flcorrvar) {
/* either x or z may be member of a group */
		remcorrvar(x, z, d) ; 
		remcorrvar(z, x, d) ; 
	}
	return TRUE ;
}


void	remcorrvar(int x, int z, int d)
/*
 * removes correlations between x and a group that contains z.
 * if d equals 1 only the j1ok correlations are removed.
 * if d equals 2 only the j2ok correlations are removed.
 */
{
	int	c, z0, s, i, j ;
	struct lcorr *pc ;
	
	if (at[x].corrgrp && at[z].ingrp) {
/* x correlates with a group and z belongs to a group */
		for (i = 0 ; i < at[x].nctot ; i++) {
			c = at[x].ex[i] ;
			pc = &propcorr[c] ;
			z0 = at[x].other[i] ;
			if (valcorr[c] && (z0 < 0)) {
/* correlation c is valid and correlates x with the group z0 */
				s = lnsize[-z0] ;
				for (j = 0 ; j < s && ln[-z0][j] != z ; j++) 
					;
				if (j != s) {
/* if z belongs to this group, the x-z0 correlation must be removed */
					if( ((d == 1) && pc->j1ok) || ((d == 2) && pc->j2ok) ) {
/* the correlation is a confirmation one because its min/max coupling path lengths do not conflict... */	
						remcorr3(TRUE, c, x, z0) ;
/* ...so it is removed */
					}
				}
			}
		}
	}
}


int		inutil1(int x, int y)
/*
 * if there is already a bond between y and an atom z, an eventual
 * x-z correlation only brings a confirmation.
 * returns FALSE if there is a correlation of x with a neighbour of y
 * which is neither j1ok nor j2ok (jnok only)
 */
{
	int	i, z ;

	if(!remcorr2(x, y, 1)) return FALSE ;
/* if x and y are correlating (corrected in v3.4.1) */
	for (i = 0 ; i < at[y].nlia ; i++) {
		z = at[y].lia[i] ;
		if (z != x) {
			if(!remcorr2(x, z, 2)) return FALSE ;
		}
	}
	return TRUE ;
}


int		inutil2(int x, int z, int y)
/*
 * a new bond x-z is created. the bond z-y may or may not be already created 
 * x --------- z - - - - - y
 */
{
	int fl ;
	
	fl = inutil1(z, x) ;
/* looks for t bound to x and for a z-t correlation */
	fl = fl && remcorr2(x, y, 2) ;
/* if x and y are correlating (corrected in v3.4.1) */
	fl = fl && remcorr2(z, x, 1) ;
/* if x and z are correlating */
	fl = fl && inutil4(z, x, y) ;
/* looks for t bound to z, and to t-x and t-z correlations */
	return fl  ;
}


int		inutil4(int z, int x, int y)
/*
 * looks for t bound to z, and to t-x correlations because x-z is a new bond
 * returns FALSE if there is a correlation of x with a neighbor of z
 * that is neither j1ok nor j2ok (jnok only)
 */
{
	int	i, n, t ;

	n = at[z].nlia ;
	if (n <= 2) { 
		return TRUE ;
	}
/* if z is already bound to t (different from x and y), this point is
 * necesserely reached because z has then at least 3 bonds. It would save
 * time to ensure that for a given (x,y,z) triplet this point is only 
 * reached one time */
	for (i = 0 ; i < n ; i++) {
		t = at[z].lia[i] ;
		if (t != x && t != y) {
			if(!remcorr2(x, t, 2)) return FALSE ;
/* removes eventual x-t correlations. must be done because x-z is a new bond.
 * if the y-z bond already exists, the y-t correlations have already been removed.
 * if the y-z bond does not already exist, it will created and y-t correlations
 * will subsequently be removed */
		}
	}
	return TRUE ;
}
