/*
    file: LsdSrc/display.c
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
#include <signal.h>

extern struct atom at[MAX_ATOMES] ;
extern struct noy nucl[MAX_NUCL] ;
extern int	isol ;
extern int	iph2 ;
extern int	flverb ;
extern int	max_atomes ;
extern FILE *fic2 ;

extern void controlC(int) ;

void	display(int p) ;
void	outlist(FILE *f) ;
void	outdump(void) ;


void	display(int p)

/* 
 * p = 0 : displays as bonds list
 * p = 1 : displays complete information for OUTLSD
 */

{
	signal(SIGINT, SIG_IGN) ;
/* avoids interruption during solution printout */
	switch (p) {
	case 0 :
		outlist(fic2) ;
		if (flverb && (fic2 != stdout)) {
			outlist(stdout) ;
		}
		break ;
	case 1 : 
		outdump() ; 
		break ;
	}
	signal(SIGINT, controlC) ;
/* interruption allowed now */
}


void	outlist(FILE *f)

/*
 * displays as lists
 */

{
	int	i, x, y, n ;
	struct atom *a ;

	fprintf(f, "\n*** solution %d *** %d\n", isol, iph2) ;
	for (x = 1 ; x <= max_atomes ; x++) {
		if (at[x].utile) {
			a = &at[x] ;
			n = a->nlia ;
			for (i = 0 ; i < n ; i++) {
				y = a->lia[i] ;
				if (x < y) {
					fprintf(f, "%4d%4d%2s\n", 
					x, 
					y, 
					( (a->ordre[i] == 1) ? " -" : ((a->ordre[i] == 2) ? " =" : " #") )) ;
				}
			}
		}
	}
}


void	outdump(void)

/*
 * displays for OUTLSD
 */

{
	int	i, x ;
	struct atom *a ;

	if(isol == 1) fprintf(fic2, "OUTLSD\n") ;
	fprintf(fic2, "%d %d\n", max_atomes, isol) ;
	for (x = 1 ; x <= max_atomes ; x++) {
		a = &at[x] ;
		fprintf(fic2, " %1d", a->utile) ;
		if (a->utile == TRUE) {
			fprintf(fic2, " %2s %1d %1d %1d %1d %2d", 
				nucl[a->element].symout, 
				a->valence, 
				a->mult, 
				a->nlmax, 
				a->nlia,
				a->charge ) ;
			for (i = 0 ; i < a->nlia ; i++) {
				fprintf(fic2, " %3d %1d", 
					a->lia[i] , 
					a->ordre[i] ) ;
			}
			for (i = a->nlia ; i < 4 ; i++) {
				fprintf(fic2, " %3d %1d", 0, 0) ;
			}
		}
		fprintf(fic2, "\n") ;
	}
}


