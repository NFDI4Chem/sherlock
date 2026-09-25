/*
    file: LsdSrc/chklsd.c
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
extern struct noy nucl[MAX_NUCL] ; /* lecture.c */
extern int 	numeqdbbonds ; /* doubles.c */
extern int	max_atomes ; /* lsd.c */
extern int 	numskelatoms ; /* lsd.c */
extern int 	numskelbonds ; /* lsd.c */

int	numrings ; 

int 	chklsd(void) ;
void	setiidx(void) ;

/*
 * prints basic molecule information. returns the sum of atom
 * valencies, that should be even.
 * Defines the iidx field of atoms, used for inchi calculation.
 * iidx index start at 0 and only concern useful (utile) atoms
 * so that 0 <= iidx < numskelatoms.
 */
int chklsd(void)
{
	struct atom *a ;
	int numatomes = 0 ;
	int numbonds = 0 ;
	float curmass ;
	float totalmass = 0 ;
	int numinsat ;
	int nH ;
	int badmass = FALSE ;
	int i ;
	int n ;
	int x ;
	int totcharge = 0 ;
	int poscharge = 0 ;
	int negcharge = 0 ;
	
	setiidx() ;
/* sets the iidx field of useful (utile) atoms and calculates numskelatoms */
	numeqdbbonds /= 2 ;
/* numeqdbbonds was previously the number of sp2 equivalent atoms (an sp atom is equivalent to two sp2 atoms),
 * it must be an even number */
	for(i = 1 ; i <= max_atomes ; i++) {
		a = &at[i] ;
		if (a->utile == TRUE) {
			numbonds += a->valence ;
			if(a->charge) {
				totcharge++ ;
				(a->charge > 0) ? poscharge++ : negcharge++ ;
			}
		}
	}
	numbonds += nucl[MAX_NUCL-1].count ;
	if ((numbonds % 2) == 1) return numbonds ;
	numbonds /= 2 ;
	for (i = 0 ; i < MAX_NUCL ; i++) {
		n = nucl[i].count ;
		if (n > 0) {
			numatomes += n ;
			curmass = nucl[i].mass ;
			if (curmass < 0.0) {
				badmass = TRUE ;
				printf("Mass of %s is unknown\n", nucl[i].sym) ;
			} else {
				totalmass += n * curmass ;
			}
		}
	}
	numinsat = numbonds - numatomes + 1 ;
	numrings = numinsat - numeqdbbonds ;
	nH = nucl[MAX_NUCL-1].count ;
	numskelbonds = numbonds - numeqdbbonds - nH ;
	for (i = 0 ; i < (MAX_NUCL-1) ; i++) {
		n = nucl[i].count ;
		if (n > 0) {
			printf("%s %d ", nucl[i].sym, n) ;
		}
	}
	printf("H %d ", nH) ;
	printf("\n") ;
	printf("%d atom%s\n", numatomes, (numatomes > 1) ? "s" : "") ;
	printf("%d bond%s\n", numbonds, (numbonds > 1) ? "s" : "") ;
	printf("%d skeletal atom%s\n", numskelatoms, (numskelatoms > 1) ? "s" : "") ;
	printf("%d skeletal bond%s\n", numskelbonds, (numskelbonds > 1) ? "s" : "") ;
	printf("Degree of unsaturation: %d\n", numinsat) ;
	printf("%d ring%s\n", numrings, (numrings > 1) ? "s" : "") ;
	printf("Total number of charges: %d\n", totcharge) ;
	if(totcharge) {
		printf("%d positive charge%s\n", poscharge, (poscharge > 1) ? "s" : "") ;
		printf("%d negative charge%s\n", negcharge, (negcharge > 1) ? "s" : "") ;
	}
	if (!badmass) {
		printf("M = %.5f\n", totalmass) ;
	}
	return numbonds * 2 ;
}

void	setiidx(void)
/* sets the iidx field of useful (utile) atoms and calculates numskelatoms */
{
	int x ;
	struct atom *al ;
	
	for (x = 1 ; x <= max_atomes ; x++) {
		al = &at[x] ;
		if (!al->utile) {
			continue ;
		}
		al->iidx = numskelatoms ;
		numskelatoms++ ;
	}

}
