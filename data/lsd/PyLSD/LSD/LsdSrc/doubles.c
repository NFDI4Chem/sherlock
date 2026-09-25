/*
    file: LsdSrc/doubles.c
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
#include "stdio.h"
extern int	nivheap ; /* lsd.c */

extern void	marque(void) ;
extern void	empile0(int *a, int v) ;
extern void	relache(void) ;

extern struct atom at[MAX_ATOMES] ;
extern int	max_atomes ;
extern int	flkeep ;
extern int	flverb ; /* lsd.c */

int	dou[MAX_ATOMES] ; /* true when the sp2 atom has got its double bond or when an sp atom has got only one double bond */
int tri[MAX_ATOMES] ; /* true when the sp atom has got its triple bond */
int cum[MAX_ATOMES] ; /* true when the sp atom has got two double bond (cumulene) */
int	fldou ; /* true when the distribution of multiple bonds succeeds */
int	lan[MAX_ATOMES] ; /* logical "ancient" atom set - see abredt() */
int	lne[MAX_ATOMES] ; /* logical "new" atom set - see abredt()  */
int	lfu[MAX_ATOMES] ; /* logical "future" atom set - see abredt()  */
int numsp2 ; /* number of sp2 atoms */
int	numeqdbbonds = 0 ; /* the double of the number of double bond equivalents (a triple bond is equivalent to two double bonds) */

int	multiplelia(void) ;
void	makemb(void) ;
int		mb(int yo[]) ;
int	antibredt(void) ;
int	abredt(int a, int b, int c, int d, int e) ;

int	multiplelia(void)
/*
 * initialises atoms for multiple bonds search
 * calls makemb, that sets fldou
 */
{
	int	i, x ;

	fldou = FALSE ;
	for (x = 1 ; x <= max_atomes ; x++) {
		dou[x] = FALSE ;
		tri[x] = FALSE ;
		cum[x] = FALSE ;
		for (i = 0 ; i < at[x].nlia ; i++) {
			at[x].ordre[i] = 1 ;
		}
	}
	makemb() ;
	flkeep = FALSE ; /* heap manager back to normal state */
	return fldou ;
}


void	makemb(void)
/*
 * builds the double and triple bonds set recursively
 */
{
	int	i, j, y, z, o ;
	int yo[2] ;

	if (mb(yo)) { /* mb returns 0 when all double and triple bonds have been placed */
		y = yo[0] ; /* y is sp or sp2 and needs a multiple bond */
		o = yo[1] ; /* 2 for a double (y is sp2), 3 for a triple (y is sp) and 22 for a second double bond (y is sp) */
		for (i = 0 ; i < at[y].nlia ; i++) {
			if (fldou) {
/* a solution has been found, nothing to do anymore */
				break ;
			}
			z = at[y].lia[i] ; 
/* z is y's ith neighbour */
			if ( (at[z].hybrid == 0) || ((at[z].hybrid == 1) && dou[z]) || ((at[z].hybrid == 2) && tri[z]) || ((at[z].hybrid == 2) && cum[z]) ) {
/* z has to be sp2 and not involved in a double bond yet
 * or sp and not involved in a triple bond yet
 * or sp and not involved in two double bonds yet */
				continue ;
			}
			for (j = 0 ; (j < at[z].nlia) && (at[z].lia[j] != y) ; j++)
/* y is z's jth neighbour */
				;
			marque() ;
/* a new group in the heap : y-z and z-y bonds are double or triple
 * y and z are involved in a multiple bond */
			if(o == 2) {
/* double bond between y and z. y is sp2, z may be sp2 or sp */
				if(flverb > 1) {
					printf("%d : double bond between %d and %d\n", nivheap, y, z) ;
				}
				empile0(&(at[y].ordre[i]), 2) ;
				empile0(&(at[z].ordre[j]), 2) ;
				empile0(&(dou[y]), TRUE) ;
				empile0(&(dou[z]), TRUE) ;
			} else if(o == 3) {
/* triple bond between y and z, two sp atoms */
				if(flverb > 1) {
					printf("%d : triple bond between %d and %d\n", nivheap, y, z) ;
				}
				empile0(&(at[y].ordre[i]), 3) ;
				empile0(&(at[z].ordre[j]), 3) ;
				empile0(&(tri[y]), TRUE) ;
				empile0(&(tri[z]), TRUE) ;
			} else {
/* double bond (cumulene) between y and z. y is sp, z may be sp or sp2 */
/* y is now involved in two double bonds */
				if(flverb > 1) {
					printf("%d : double bond between %d and %d\n", nivheap, y, z) ;
				}
				empile0(&(at[y].ordre[i]), 2) ;
				empile0(&(at[z].ordre[j]), 2) ;
				empile0(&(cum[y]), TRUE) ;
				empile0(&(dou[z]), TRUE) ;
			}
			makemb() ;
/* find the next multiple bond */
			relache() ;
/* releases the data from the heap unless flkeep has become true, as it
 * happens when all double and triple bonds have been found */
		}
	} else {
/* preserves information upon backtracking
 * multiple bonds have all been found */
		flkeep = fldou = TRUE ;
	}
}


int		mb(int yo[])
/*
 * 1) Tries first to look for an sp atom y, involved in only one double bond. 
 *    This is intended to discover cumulene-type systems.
 * 2) Tries to look for an atom y bound to x, involved in a double bond.
 *    y must be sp2 and not involved in a double bond.
 *    This is intended to discover the systems of conjugated double bonds.
 * 3) If no such x can be found, y is chosen as an sp2 atom not involved
 *    in a double bond.
 * 4) and 5) idem as 2) and 3) with triple bonds and sp atoms.
 *    If no y can be found the job is done.
 */
{
	int	go, x, y, i ;
	
	for (y = 1 ; y <= max_atomes ; y++) {
		if (at[y].hybrid != 2 || !dou[y] || cum[y]) {
			continue ;
		}
		yo[0] = y ;
		yo[1] = 22 ;
		return TRUE ;
/* y is sp and is involved in only one double bond */
	}
	for (x = 1 ; x <= max_atomes ; x++) {
		if (!dou[x]) {
			continue ;
		}
/* look for an atom y bound to x, involved in a double bond */
		for (i = 0 ; i < at[x].nlia ; i++) {
			y = at[x].lia[i] ;
			go = dou[y] || (at[y].hybrid != 1) ;
			if (!go) {
				yo[0] = y ;
				yo[1] = 2 ;
				return TRUE ;
			}
		}
/* y is sp2 and is not involved in a double bond */
	}
	for (y = 1 ; y <= max_atomes ; y++) {
		go = dou[y] || (at[y].hybrid != 1) ;
		if (!go) {
			yo[0] = y ;
			yo[1] = 2 ;
			return TRUE ;
		}
/* y is sp2 and is not involved in a double bond */
	}
	for (x = 1 ; x <= max_atomes ; x++) {
		if (!tri[x]) {
			continue ;
		}
/* look for an atom y bound to x, involved in a triple bond */
		for (i = 0 ; i < at[x].nlia ; i++) {
			y = at[x].lia[i] ;
			go = tri[y] || (at[y].hybrid != 2) ;
			if (!go) {
				yo[0] = y ;
				yo[1] = 3 ;
				return TRUE ;
			}
		}
/* y is sp and is not involved in a triple bond */
	}
	for (y = 1 ; y <= max_atomes ; y++) {
		if(!tri[y] && (at[y].hybrid == 2) && !dou[y] && !cum[y]) {
			yo[0] = y ;
			yo[1] = 3 ;
			return TRUE ;
		}
/* y is sp and is not involved in a triple bond */
	}
	return FALSE ;
/* no y can be found */
}


int	antibredt(void)
/* 
 * looking for a-b=c(-d1)-d2) arrangements, and test them against Bredt'rule
 * see abredt
 */

{
	int	atab[3], a, ia, na ;
	int	b, nb ;
	int	ctab[2], c, ic, nc ;
	int	dtab[3], d, d2, id, id2, nd ;
	int	i ;

	for (b = 1 ; b <= max_atomes ; b++) {
		if (!at[b].utile || at[b].hybrid == 0 || (nb = at[b].nlia) < 2) {
			continue ;
		}
/* b has nb neighbours with nb >=2 */
		for (i = 0, nc = 0 ; i < nb ; i++) {
			if (at[b].ordre[i] == 2) {
/* looking for the double bonds of atom b */
				c = at[b].lia[i];
				if (at[c].nlia >= 3) {
					ctab[nc++] = c ;
/* atoms in ctab share a double bond with b and have at least three neighbors */
				}
			}
		}
		if (nc == 0) {
			continue ;
		}
/* at least one double bond of b with elements in ctab, nc such elements */
		for (ic = 0 ; ic < nc ; ic++) {
			c = ctab[ic] ;
/* considering now the b=c double bond */
			if (flverb > 1) {
				printf("Checking Bredt's rule for %d=%d\n", b, c) ;
			}
			for (i = 0, na = 0 ; i < nb ; i++) {
				a = at[b].lia[i] ;
				if (a != c && at[a].nlia > 1) {
					atab[na++] = a ;
				}
			}
			if (na < 1) {
				continue ;
			}
/* atab contains the neighbors of b that are different from c, na such atoms
with more than one neighbor, na at least equal to 1 */
			for (i = 0, nd = 0 ; i < at[c].nlia ; i++) {
				d = at[c].lia[i] ;
				if (d != b && at[d].nlia > 1) {
					dtab[nd++] = d ;
				}
			}
			if (nd < 2) {
				continue ;
			}
/* dtab contains the neighbors of c that are different from b, nd such atoms
with more than one neighbor, nd at least equal to 2 */
			for (ia = 0 ; ia < na ; ia++) {
				a = atab[ia] ;
/* a is a neighbor of b */
				for (id = 1 ; id < nd ; id++) {
					d = dtab[id] ;
/* d is a neighbor of c */
					for (id2 = 0 ; id2 < id ; id2++) {
						d2 = dtab[id2] ;
/* d2 is another neighbor of c */
						if (abredt(a, b, c, d, d2)) {
							return TRUE ;
						}
					}
				}
			}
		}
	}
/* the structure is not antibredt. It is a good one */
	return FALSE ;

}


int	abredt(int a, int b, int c, int d, int e)
/*
 * looks for all the set (lan) of atoms separated from a by 4 or less bonds.
 * if d and e are simultaneously in this set, it means that the double bond
 * b=c is embedded in 2 rings of size less or equal to 7. The structure
 * is then declared as antibredt.
 */
{
	int	i, j, x, y, cf ;

	if (flverb > 1) {
		printf ("Checking Bredt's rule for %d-%d=%d(-%d)-%d\n", a, b, c, d, e) ;
	}
	for (x = 1 ; x <= max_atomes ; x++) {
		lan[x] = lne[x] = FALSE ;
	}
	lne[a] = lan[a] = TRUE ;
/* lan and lne only contain a */
	for (i = 0 ; i < 4 ; i++) { /* 4 times the step of set enlargement */
		for (x = 1 ; x <= max_atomes ; x++) {
			lfu[x] = FALSE ;
/* reset the list of future set members */
		}
		cf = 0 ;
/* counts the number of elements in lfu */
		for (x = 1 ; x <= max_atomes ; x++) {
			if (lne[x]) {
/* x is an atom that has entered in lan during the last iteration
 * or it is a (for the first iteration */
				for (j = 0 ; j < at[x].nlia ; j++) {
					y = at[x].lia[j] ;
/* y is its jth neighbour */
					if ((y != b) && (y != c) && (!lan[y])) {
						lfu[y] = lan[y] = TRUE ;
						cf++;
					}
/* y is neither in lan nor it is b or c. it can enter into lan. */
				}
			}
		}
		if (cf == 0) {
			break ;
/* no new atom during the iteration. No need to go further */
		}
		for (x = 1 ; x <= max_atomes ; x++) {
			lne[x] = lfu[x] ;
/* the list of the atoms discovered during the iteration will be the
 * list of new atoms for the next iteration */
		}
	}
	return lan[d] && lan[e] ;
/* d and e must not be in lan, else abredt returns true */
}


