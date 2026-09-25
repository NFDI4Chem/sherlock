/*
    file: LsdSrc/equiv.c
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
extern int	max_atomes ; /* lsd.c */
extern int	flCclasses ; /* lsd.c */
extern int	fldupl ; /* lsd.c */

extern void	myexit(int i) ; /* lsd.c */
extern void	empile0(int *a, int v) ; /* heap.c */

struct co codes[MAX_EQU] ;
int	max_codes = 0 ;

void	initequ(void) ;
int	calcstat(struct atom *a) ;
int	testequ(int x) ;
void	modequ(int x) ;

void	initequ(void)
/*
 * Two or more unbound uncorrelating atoms of identical status can be
 * interchanged without modifing the currently growing structure. This
 * will result in chemically identical structures only differing
 * by the atoms numbering. In order to save time, these atoms have to be
 * grouped and only one at a time can be involved in the bond formation
 * process.
 * Initequ builts equivalence classes. Unbound uncorrelating atoms of identical
 * status are grouped toghether. Classes are numbered from 0. at[x].equstat
 * is this class number, unless x is bound or correlating or alone in its
 * class. The class number gives access in codes[] to the class's atoms status
 * (valco) and to the number of atoms already belonging to the class (ico).
 * at[x].ieq is the index of x within its class, starting to 1.
 * Once all atoms have been inspected with respect to the status equivalence,
 * all the codes[i].ico are reset to 1. codes[i].ico will be used afterwards
 * to tell which is the atom within the ith class that can be used in order
 * to establish bonds.
 */
{
	int	i, x, s, *e ;
	struct atom *a ;

/*
	if((fldupl > 0) && (flCclasses == FALSE)) {
		fldupl = 1 ;
	}
*/
/* avoids the elimination of duplicated solution when carbon classes are not allowed */
	for (x = 1 ; x <= max_atomes ; x++) {
		if (!at[x].utile) {
			continue ;
		}
		a = &at[x] ;
		if (((a->element == 0) && (flCclasses == FALSE)) || (a->nc != 0) || (a->nlia != 0) || a->ingrp) {
			a->equstat = -1 ;
/* x is carbon and no carbon classes is allowed or x is bound or x correlates */
		} else {
/* x neither bound nor correlating */
			s = calcstat(a) ;
/* coded status */
			for (i = 0 ; i < max_codes && codes[i].valco != s ; i++)
				;
			if (i == max_codes) {
/* unknown coded status */
				codes[i].valco = s ;
/* new code */
				codes[i].ico = 1 ;
/* first occurence of this code */
				a->equstat = i ;
/* storing x's class number */
				a->iequ = 1 ;
/* x's index within the class */
				if (max_codes == MAX_EQU) {
					printf("increase MAX_EQU \n") ;
					myexit(-1) ;
/* too many classes */
				}
				max_codes++;
/* one more class */
			} else {
/* x belongs to an already defined class */
				a->equstat = i ;
/* x belongs to the ith class */
				a->iequ = ++codes[i].ico ;
/* x's index within this class */
			}
		}
	}

	for (x = 1 ; x <= max_atomes ; x++) {
		if (at[x].utile) {
			e = &at[x].equstat ;
			if ((*e >= 0) && (codes[*e].ico == 1)) {
/* at[x].equstat turns to -1 if it is alone in its class */
				*e = -1 ;
			}
		}
	}

	for (i = 0 ; i < max_codes ; i++) {
		codes[i].ico = 1 ;
/* only the fist element of each class is accessible for bonding processes */
	}
}


int	calcstat(struct atom *a) 	
/*
 * translates x's status in a numerical code 
 */
{
	return a->element + 20 * (a->hybrid + 20 * (a->mult + 20 * a->charge)) ;
}


int	testequ(int x)
/*
 * tells if x can be involved in a bond due to equivalence constraints.
 * does not modify equivalence classes
 */
{
	struct atom *a ;

	a = &at[x] ;
	if (a->equstat < 0) {
		return TRUE ;
/* no equivalence constraint for x */
	}
	return (a->iequ == codes[a->equstat].ico) ;
/* returns true if x is presently the atom of lowest index within its 
 * equivalence class */
}


void	modequ(int x)
/*
 * x is now bound. modequ modifies its eventual equivalence class 
 */
{
	struct atom *a ;
	int	*e ;

	a = &at[x] ;
	if (a->equstat < 0) {
		return ;
/* x is not equivalent to an other atom : nothing to do */
	}
	e = &(codes[a->equstat].ico) ;
	empile0(e, (*e) + 1) ;
/* allows the next member of the class to be involved in bonds */
	empile0(&(a->equstat), -1) ;
/* x is not anymore equivalent to an other atom */
}
