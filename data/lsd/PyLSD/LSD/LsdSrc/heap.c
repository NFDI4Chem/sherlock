/*
    file: LsdSrc/heap.c
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
extern struct pile heap[MAX_HEAP] ; /* lsd.c */
extern int	nivheap ; /* lsd.c */
extern int	flhistoire ; /* lsd.c */
extern int	supniv ; /* lsd.c */
extern int	nivdisp ; /* lsd.c */
extern int	fldisp ; /* lsd.c */

extern void	display(int p) ; /* display.c */
extern void	wrln2(int y) ; /* start1.c */
extern void	myexit(int i) ; /* lsd.c */

int	flkeep = FALSE ; /* controls the behaviour of the heap manager */
int	iheap = 0 ; /* heap index */

void	marque(void) ;
void	empile0(int *a, int v) ;
void	empile(int t, int a, int v) ;
void	empile6(int a, int v, int c) ;
void	relache(void) ;
void	histoire(void) ;


void	marque(void)
/*
 * the lsd algorithm is recursive in its very nature. In order to prevent
 * the duplication of the structure's state at each recursive call,
 * the state is kept in global variables. The values that must be
 * recovered when lsd is backtracking are pushed into the heap. The heap
 * is also used to store choices that have been made since the beginning of
 * the resolution. These informations are listed upon call to histoire().
 * All the calls to empile (means "push to the top of the heap") taking
 * place between marque() ("mark") and relache ("release") are prone to
 * be released simultaneously. When double bonds are placed there is no point
 * to release the modified values. Only when the global flag flkeep is false,
 * the modified values are put back to their initial location. The heap is an
 * array of records of type pile (means "stack"). The filling of the heap
 * is controlled by the iheap index. The ptyp field is 0 if the heap is used
 * for recursivity simulation, and different from 0 if only history
 * information are pushed onto. The niv field is the same for all items
 * pushed between marque() and relache(). the pinfo field is a union
 * of records pmem (ptyp is 0) and phist (else). In pmem adr strores the
 * address of the value to be modified and val is the value to be preserved.
 * in phist h1 and h2 contains the values the user wants to store.
 */
{
	nivheap++;
/* new bloc of information in the heap */
	if (nivheap > supniv) {
		supniv = nivheap ;
/* the highest bloc index */
	}
	if (nivdisp && (supniv == nivdisp)) {
/* stops according to the MLEV command */
		display(fldisp) ;
/* displays the current structure */
		if (flhistoire) {
/* and history upon request */
			histoire() ;
		}
		myexit(1) ;
/* and stops here */
	}
}

/*
 * empile0 is new with version 3.2.1. 
 * For compatibility with 64 bit computers in which
 * the size of an (int) is not equal to the size of an (int *).
 * Calls of empile with 0 as first argument are replaced
 * by calls to empile0 throughout the code.
 */ 
void	empile0(int *a, int v)
/*
 * pushes an address/value pair on the stack
 */
{
	iheap++;
/* increments heap index */
	if (iheap == MAX_HEAP) {
/* the heap is too small */
		printf("stack is full - increase MAX_HEAP\n") ;
		myexit(-1) ;
	}
	heap[iheap].pinfo.pmem.adr = a ;
	heap[iheap].pinfo.pmem.val = *a ;
/* substitutes the old value at address a by the new one (v) */
	*a = v ;
	heap[iheap].niv = nivheap ;
/* stores bloc number */
	heap[iheap].ptyp = 0 ;
/* stores the information type */
}


void	empile(int t, int a, int v)
/*
 * pushes a pair of integer on the stack
 */
{
	iheap++;
/* increments heap index */
	if (iheap == MAX_HEAP) {
/* the heap is too small */
		printf("stack is full - increase MAX_HEAP\n") ;
		myexit(-1) ;
	}
/* a and v are stored into the heap */
	heap[iheap].pinfo.phist.h1 = a ;
	heap[iheap].pinfo.phist.h2 = v ;
	heap[iheap].niv = nivheap ;
/* stores bloc number */
	heap[iheap].ptyp = t ;
/* stores the information type */
}


void	empile6(int a, int v, int c)
/*
 * pushes a pair of integer and 
 * the corresponding correlation index on the stack.
 * Calls of empile with 6 as first argument are replaced
 * by calls to empile6 throughout the code.
 */
{
	iheap++;
/* increments heap index */
	if (iheap == MAX_HEAP) {
/* the heap is too small */
		printf("stack is full - increase MAX_HEAP\n") ;
		myexit(-1) ;
	}
/* a, v and c are stored into the heap */
	heap[iheap].pinfo.phist.h1 = a ;
	heap[iheap].pinfo.phist.h2 = v ;
	heap[iheap].pinfo.phist.corr = c ;
	heap[iheap].niv = nivheap ;
/* stores bloc number */
	heap[iheap].ptyp = 6 ;
/* stores the information type */
}


void	relache(void)
/*
 * pops
 */
{
	int	*a, v ;

	while (heap[iheap].niv == nivheap) {
/* deals only with the most recent bloc of informations */
		if ((!heap[iheap].ptyp) && (!flkeep)) {
/* goes back to previous state only for ptyp 0 information and when
 * flkeep is false */
			a = heap[iheap].pinfo.pmem.adr ;
			v = heap[iheap].pinfo.pmem.val ;
			*a = v ;
/* done */
		}
		iheap--;
/* the heap index is decremented */
	}
	nivheap--; 
/* the block has been completely removed */
}


void	histoire(void)
/*
 * goes through the heap from bottom to top in order to list chronologically
 * the events of the resolution 
 */
{
	struct pilehist *a ;
	int	i ;

	for (i = 1 ; i <= iheap ; i++) {
/* from bottom to top */
		a = &heap[i].pinfo.phist ;
		switch (heap[i].ptyp) {
/* history items have ptype in the range 1-6 */
		case 1 :
			printf("\nuses correlation %d - ", a->h1) ;
			wrln2(a->h2) ;
			break ;
		case 2 :
			printf("binds atoms %d and %d\n", a->h1, a->h2) ;
			break ;
		case 3 :
			printf("confirms by correlation %d - ", a->h1) ;
			wrln2(a->h2) ;
			break ;
		case 4 :
			printf("\nbinds atoms %d and %d\n", a->h1, a->h2) ;
			break ;
		case 5 :
			printf("hypothetic correlation %d - %d\n", a->h1, a->h2) ;
			break ;
		case 6 :
			printf("\neliminates correlation %d - ", a->h1) ;
			wrln2(a->h2) ;
			break ;
		}
/* self-explanatory */
	}
}
