/*
    file: LsdSrc/connexe.c
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
#define DEBUG
*/

#include "defs.h"

#ifdef DEBUG
#include <stdio.h>
#endif

int	setcnx(int x, int y, int fl) ;

extern void	non_nul(int errnum, int v, ...) ; /* lecture.c */
extern void	empile0(int *a, int v) ; /* heap.c */
extern int	insatur(int x) ; /* lsd.c */

extern struct	atom at[MAX_ATOMES] ; /* lsd.c */
extern int	max_atomes ; /* lsd.c */
extern int	flcnx ; /* lsd.c */

int	ncnx ; /* next available connected unit index */
int	natomes ; /* the number of atoms in the solution(s) */

/*
 * Udates atom connectivity data.
 * The "connectivity" of a solution is true or false
 * depending on whether it is in one piece or not.
 * The French word for that is "connexity".
 * Do not mix "connected unit" (here) and "connected set" (in lsd.c).
 * setcnx() always returns TRUE if flcnx is FALSE (connectivity checking off).
 * If flcnx is TRUE, returns TRUE if the binding of 
 * atoms x and y does not create a fragment in which 
 * all atoms are complete.
 * If fl is TRUE then all variable modifications are backtraking.
 */
int	setcnx(int x, int y, int fl)
{
	int cx, cy, cz ;
	int i ;
	int z ;
	int good ;
	int nat ;
	int result ;
	
	if(!flcnx) {
		return TRUE ;
	}
/* no connectivity constraint if flcnx is false */
	cx = at[x].icnx ;
	cy = at[y].icnx ;
#ifdef DEBUG
fprintf(stderr, "*** connecting %d(%d) with %d(%d)\n", x, cx, y, cy) ;
#endif
	good = FALSE ;
	nat = 0 ;
	for(z = 1 ; z <= max_atomes ; z++) {
		cz = at[z].icnx ;
#ifdef DEBUG
fprintf(stderr, "what about %d(%d)?\n", z, cz) ;
#endif
		if(cz == cx || cz == cy) {
/* atom z is in the connected unit of either x or y. */
			nat++ ;
#ifdef DEBUG
fprintf(stderr, "adding %d to connected unit %d, that has now %d atoms\n", z, ncnx, nat) ;
#endif
/* one more atom in the new connected unit */
			good = good || insatur(z) ;
#ifdef DEBUG
fprintf(stderr, "%d: nlmax: %d nlia: %d\n", z, at[z].nlmax, at[z].nlia) ;
fprintf(stderr, "good is now %s\n", good ? "TRUE" : "FALSE") ;
#endif
/* good becomes true if at least one atom that will get ncnx as icnx is
   not complete */
#ifdef DEBUG
fprintf(stderr, "backtracking is %s\n", fl ? "TRUE" : "FALSE") ;
#endif
			if(fl) {
				empile0(&at[z].icnx, ncnx) ;
			} else {
				at[z].icnx = ncnx ;
			}
		}
	}
	if(fl) {
		empile0(&ncnx, ncnx+1) ;
	} else {
		ncnx++ ;
	}
/* sets the next available connected unit index */
	result = good || (nat == natomes) ;
#ifdef DEBUG
fprintf(stderr, "result is %s\n", result ? "TRUE" : "FALSE") ;
#endif
	return result ;
/* if the number of atoms in the connected unit of x and y is equal
   to the total number of atoms, there cannot be two connected units
   and setcnx() must return TRUE */
}
