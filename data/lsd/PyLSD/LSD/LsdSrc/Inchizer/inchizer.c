/*
    file: LsdSrc/Inchizer/inchizer.c
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "inchi_api.h"
#include "e_inchi_atom.h"

/*
 * See inchizer.h for the interface definition of inchizer
 */
static inchi_Input Inp ;
static inchi_Output Out ;

int 	newInChI(int n) ;
void 	addAtom(int ia, char *elt, int charge) ;
void 	addBond(int x, int y, int bondtype) ;
void	remBonds() ;
int	stringify(char **ps, char **paux) ;
void	cleanInpInChI(void) ;
void	cleanOutInChI(void) ;

static void	bond_from_x_to_y(inchi_Atom *lx, int y, int bondtype) ;

int 	newInChI(int n)
/* returns 0 if memory allocation for n atoms failed */
{
	inchi_Atom *lat, *la ;
	int i, j ;
	
	Inp.atom = lat = e_CreateInchi_Atom(n) ;
	if (lat == NULL) {
		return 0 ;
	}
/* memory allocation from atoms */
	for(i=0, la=lat ; i<n ; i++,la++) {
		la->x = la->y = la->z = 0.0 ;
		la->num_bonds = 0 ;
		la->radical = INCHI_RADICAL_NONE ;
		la->isotopic_mass = 0 ;
		for(j = 1 ; j < NUM_H_ISOTOPES+1 ; j++) {
			la->num_iso_H[j] = 0 ;
		}
		la->num_iso_H[0] = -1 ;
	}
/* defaut values of atoms */
	Inp.num_atoms = n ;
	Inp.num_stereo0D = 0 ;
	Inp.stereo0D = NULL ;
	Inp.szOptions = "-SNon" ;
/* InChI computation parameters */
	return 1 ;
}

void	addAtom(int x, char *elt, int charge)
/* add new atom number x of given element elt and charge charge */
{
	inchi_Atom *la ;
	
	la = Inp.atom + x ;
	strcpy(la->elname, elt) ;
	la->charge = charge ;
}

void	addBond(int x, int y, int bondtype)
/* add a new bond between atoms x and y of type bondtype (1: single, 2:double, 3:triple) */
{
	bond_from_x_to_y(Inp.atom + x, y, bondtype) ;
	bond_from_x_to_y(Inp.atom + y, x, bondtype) ;
}

void	remBonds()
/* reset bond count to 0 for each atom */
{
	inchi_Atom *la ;
	int n ;
	int i ;
	
	n = Inp.num_atoms ;
	la = Inp.atom ;
	for(i=0 ; i<n ; i++,la++) {
		la->num_bonds = 0 ;
	}
}

static void	bond_from_x_to_y(inchi_Atom *lx, int y, int bondtype)
/* make of y a new neigbor of atom x with bondtype multiplicity */
{
	int nb ;
	
	nb = lx->num_bonds ;
	lx->neighbor[nb] = y ;
	lx->bond_type[nb] = bondtype ;
	lx->bond_stereo[nb] = INCHI_BOND_STEREO_NONE ;
	lx->num_bonds = nb + 1 ;
}

int	stringify(char **ps, char **paux)
/* returns 0 if it was possible to put the InChI of the current molecule at *ps.
Otherwise returns the InChI error code and set *ps to an error message.
*paux contains auxiliary information (numbering, equivalence classes) */
{
	int retval ;
	char *what ;
	
	retval = GetINCHI (&Inp, &Out) ;
	what = ((retval==0) || (retval==1)) ? Out.szInChI : Out.szMessage;
/*
fprintf(stderr, "InChI:   %s\n", Out.szInChI) ;
fprintf(stderr, "Message: %s\n", Out.szMessage) ;
fprintf(stderr, "AuxInfo: %s\n\n", Out.szAuxInfo) ;
*/
	*ps = what ;
	*paux = Out.szAuxInfo ;

	return retval ;
}

void	cleanOutInChI(void)
/* frees the InchI output data */
{
	FreeINCHI(&Out) ;
}

void	cleanInpInChI(void)
/* frees the InchI input data */
{
	e_FreeInchi_Input(&Inp) ;
}
