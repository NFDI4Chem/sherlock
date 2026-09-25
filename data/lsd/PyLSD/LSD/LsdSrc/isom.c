/*
    file: LsdSrc/isom.c
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

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "Uniciter/uniciter.h"
#include "Inchizer/inchizer.h"

static void	*pu1 = NULL ; /* pointer to uniciter when DUPL is 1 */ 
static void	*pu2 = NULL ; /* pointer to uniciter when DUPL is 2 (for inchi) */ 
static char	*solclas = NULL ; /* compressed solution */

int	isom(void) ;
int	identical(void) ;
void	isom_free(void) ;
static void	err(int e) ;
static void	feedinchiAtoms(void) ;
static void	feedinchiBonds(void) ;
static void	feedsolclas(void) ;

extern void	myexit(int i) ; /* lsd.c */
extern void	non_nul(int errnum, int v, ...) ; /* lecture.c */

extern struct atom at[MAX_ATOMES] ; /* lsd.c */
extern struct noy nucl[MAX_NUCL] ; /* lecture.c */
extern int	max_atomes ; /* lsd.c */
extern int	numskelatoms ; /* lsd.c */
extern int	numskelbonds ; /* lsd.c */
extern void	classe(int n, int t[]) ; /* lsd.c */


static void err(int e)
/* print an error message and exit */
{
	char *msgs[] = {"Memory allocation failed in DUPL.",
			"AVL internal error."} ; 
	char *msg ;
	switch (e) {
	case 1:
		msg = msgs[0] ;
		break ;
	case 2:
	case 3:
		msg = msgs[1] ;
		break ;
	}
	fprintf(stderr, "%s (status is %d)\n", msg, e) ;
	myexit(-1) ;
}

int	isom(void)
/* returns TRUE only if the current structure is an isomorph of a previous one */
{
	int status ;
	char *istring, *iaux ;
	int iresult ;
	int uniq ;
	
	if (pu2 == NULL) {
/* first use of isom() */
		pu2 = new_uniciter() ;
		if(status = uniciter_status()) err(status) ;
/* new AVL tree for old structure InChI search */
		/* newInChI(numskelatoms) ; */
/* allocates atoms as input to InChI */
		/* feedinchiAtoms() ; */
/* fills the fixed part of the InChI input */
	}
	/* modif */
	newInChI(numskelatoms) ;
	feedinchiAtoms() ;
	/* modif */
	feedinchiBonds() ;
/* fills the bond part of the InChI input */
	iresult = stringify(&istring, &iaux) ;
/* get InChI string in istring and error code in iresult */
	non_nul(506, (iresult == 0) || (iresult == 1), istring) ;
/* error: no InChI has been created */
#ifdef DEBUG
	printf("retval: %d\n", iresult) ;
	printf("result: %s\n", istring) ;
#endif
	uniq = is_unique(pu2, istring) ;
	if (status = uniciter_status()) err(status) ;
/* tries to insert current InChI string in the AVL tree and check status */
#ifdef DEBUG
	printf("is %s\n", uniq ? "new" : "old") ;
#endif
	cleanOutInChI() ;
/* free all the strings allocated by InChI for result output */	
	/* modif */
	cleanInpInChI() ;
	/* modif */
	return uniq ? FALSE : TRUE ;
/* isom() returns true if the solution has already been inserted in its AVL tree */
}

int	identical(void)
/* returns true is the current solution has already been found */
{
	int status ;
	int uniq ;
	
	if (pu1 == NULL) {
/* first use of identical */
		pu1 = new_uniciter() ;
		if (status = uniciter_status()) err(status) ;
/* create a new AVL tree */
		solclas = (char *)malloc(2 * numskelbonds + 1) ;
		if (solclas == NULL) err(1) ;
/* create room for a string version of a solution */
	}
	feedsolclas() ;
/* convert the current solution to a string */
	uniq = is_unique(pu1, solclas) ;
	if (status = uniciter_status()) err(status) ;
/* try to insert the current solution string in its AVL tree */
	return uniq ? FALSE : TRUE ;
/* identical() returns true if the solution has already been inserted in its AVL tree */
}

void	isom_free(void)
/* frees localy allocated memory */
{
	free_uniciter(pu1) ;
	if (solclas != NULL) free(solclas) ;
	free_uniciter(pu2) ;
	cleanInpInChI() ;
}

static void	feedinchiAtoms(void)
/* fill InChI input with atom data */
{
	int x ;
	struct noy *pnuc ;
	struct atom *ax ;
	
	for (x = 1 ; x <= max_atomes ; x++) {
		ax = &at[x] ;
		if (!ax->utile) {
			continue ;
		}
		pnuc = nucl + ax->element ;
		addAtom(ax->iidx, pnuc->symout, ax->charge) ;
	}	
}

static void	feedinchiBonds(void)
/* fill InChI input with bond data */
{
	int x, y, aliidx ;
	int i ;
	struct atom *ax ;
	
	remBonds() ;
/* remove all already entered bonds, if any */
	for (x = 1 ; x <= max_atomes ; x++) {
		ax = &at[x] ;
		if (!ax->utile) {
			continue ;
		}
		aliidx = ax->iidx ;
		for (i = 0 ; i < ax->nlia ; i++) {
			y = ax->lia[i] ;
			if (y < x) {
				addBond(at[y].iidx, aliidx, ax->ordre[i]) ;
			}
		}
	}	
}

static void	feedsolclas(void)
/* converts current solution to a character string */
{
	int x, y, clas[4] ;
	struct atom *ax ;
	int i, j ;
	int icl ;


	icl = 0 ;
	for (x = 1 ; x <= max_atomes ; x++) {
		ax = &at[x] ;
		if (ax->utile) {
			for (i = 0, j = 0 ; i < at[x].nlia ; i++) {
				if (x < (y = ax->lia[i])) {
					clas[j++] = y ;
				}
			}
			classe(j, clas) ;
/* x's neighbours have been ranked in ascending order */
			for (i = 0 ; i < j ; i++) {
				solclas[icl++] = (ucar)x ;
				solclas[icl++] = (ucar)clas[i] ; 
			}
		} /* ordered solution */
	}
	solclas[icl] = 0 ;
/* end of string solclas */
}
