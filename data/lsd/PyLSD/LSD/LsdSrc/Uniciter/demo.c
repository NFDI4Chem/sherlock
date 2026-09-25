/*
    file: LsdSrc/Uniciter/demo.c
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
#include <stdlib.h>
#include "uniciter.h"

/*
 * demo arg1 arg2 arg3 ... argn
 * tells for each argi if there is an argj (j<i) identical to argi (old arg)
 * or not (new arg)
 */

char *msgs[] = {"memory allocation failed.", "uniciter internal error"} ; 
void err(int e)
/* print an error message and exit */
{
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
	exit(1) ;
}

int main(int argc, char *argv[])
{
	void *pu = NULL ;
	int i ;
	char *s ;
	int uniq ;
	int status ;
	
	pu = new_uniciter() ;
	status = uniciter_status() ;
	if(status) err(status) ;
/* get a new uniciter address and status checking */
	for(i=1 ; i<argc ; i++) {
		s = argv[i] ;
		uniq = is_unique(pu, s) ;
/* uniq is true if s has still not been seen on the command line */
		status = uniciter_status() ;
		if(status) err(status) ;
/* check status */
		printf("%s is %s\n", s, uniq ? "new" : "old") ;
	}
	
	free_uniciter(pu) ;
	status = uniciter_status() ;
	if(status) err(status) ;
/* free uniciter data and check status */
	
	return 0 ;
}
