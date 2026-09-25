/*
    file: LsdSrc/Inchizer/demo.c
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
#include "inchizer.h"

int main()
{
	char *result = NULL ;
	char *aux = NULL ;
	int retval ;
	
	newInChI(3) ;
	
	addAtom(0, "C", 0) ;
	addAtom(1, "C", 0) ;
	addAtom(2, "O", 0) ;
	addBond(0, 1, 1) ;
	addBond(1, 2, 1) ;
	addBond(2, 0, 1) ;
	
	retval = stringify(&result, &aux) ;
	printf("\nOxirane\n") ;
	printf("retval: %d\n", retval) ;
	printf("result: %s\n", result) ;
	cleanOutInChI() ;
	
	remBonds() ;
	addBond(0, 1, 1) ;
	addBond(1, 2, 1) ;

	retval = stringify(&result, &aux) ;
	printf("\nEthanol\n") ;
	printf("retval: %d\n", retval) ;
	printf("result: %s\n", result) ;
	cleanOutInChI() ;

	cleanInpInChI() ;

	return 0 ;
}
