/*
    file: LsdSrc/Inchizer/inchizer.h
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
 * The inchizer library is a ultra-thin software layer above the InChI library
 * (see inchi_api.h)
 */
 
extern int 	newInChI(int n) ;
/* returns 0 if memory allocation for n atoms failed */

extern void 	addAtom(int x, char *elt, int charge) ;
/* add new atom number x of given element elt and charge charge */

extern void 	addBond(int x, int y, int bondtype) ;
/* add a new bond between atoms x and y of type bondtype (1: single, 2:double, 3:triple) */

extern void	remBonds() ;
/* reset bond count to 0 for each atom */

extern int	stringify(char **ps, char **paux) ;
/* returns 0 if it was possible to put the InChI of the current molecule at *ps.
Otherwise returns the InChI error code and set *ps to an error message.
*paux contains auxiliary information (numbering, equivalence classes) */

extern void	cleanInpInChI(void) ;
/* frees the InchI input data */

extern void	cleanOutInChI(void) ;
/* frees the InchI output data */
