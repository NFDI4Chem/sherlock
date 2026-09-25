/*
    file: Genpos/genpos.c
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
 * This programs builts a postscript file from 2D coordinates generated
 * by the "outlsd 6 " command.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "defs.h"

#define SC 300 
#define P2 0.01
#define CORPS 20
#define TRUE 1
#define FALSE 0

FILE *fp ;
int	isol ;
int	max_atomes ;
struct at3D at[MAX_AT3D] ;
int	ides = 0 ;
int ncharge ;

int	lire(void) ;
int	nombres(void) ;
void	wr2d(int n) ;
void	segment(int i, int j, int order) ;
void	numero(int i) ;
void	charge(int i) ;
void	entete(void) ;


int main(int argc, char **argv)
{
	char	magic[5] ;

	if (argc != 1) {
		fprintf(stderr, "Usage : genpos\n") ;
		return 1  ;
	}
/* genpos read and writes from and to  standard devices */
	if ( !fscanf(stdin, "%4s", magic) || strcmp(magic, "DRAW") ) {
		fprintf(stderr, "genpos: input data is not valid.\n") ;
		exit(1) ;
	}
/* input comes from "outlsd 6" output */
	while (lire()) {
		wr2d(max_atomes) ;
/* does the job */
	}
	printf("%%%%Trailer\n") ;
	printf("%%%%Pages: %d\n", ides - 1) ;
	printf("%%%%EOF\n") ;
/* final page number */
	exit(0) ;
}


int	lire(void)
{
	int	fl ;

	fl = nombres() ;
/* reads one molecule */
	if ((!fl) && (ides == 0)) {
		fprintf(stderr, "genpos: input data is not valid.\n") ;
		exit(1) ; 
	}
/* an input file should at least contain one molecule before the EOF */
	if (ides == 0) {
		entete() ;
/* Postscript header */
	}
	ides++;
	return fl ;
}


int	nombres(void)
/* reads a molecule */
{
	int	i, x, v, dum ;
	struct at3D *a ;
	
	ncharge = 0 ;
/* number of atom with charge initialisation */
	v = fscanf(stdin, "%d %d", &max_atomes, &isol) ;
/* number of atoms and solution's index (used as title) */
	if (!v || !max_atomes) {
		return FALSE ;
	}
	for (x = 0 ; x < max_atomes ; x++) {
		a = &at[x] ;
		v = fscanf(stdin, "%d %d %2s %d %d", &dum, &a->loc, a->nom, &a->charge, &a->nlien) ;
/* atom order number, lsd atom reference, symbol, charge, number of bonds */
		if (v != 5) {
			return FALSE ;
		}
		if (a->charge) {
			ncharge++ ;
/* one more atom with charge */
		}
		for (i = 0 ; i < 6 ; i++) {
			v = fscanf(stdin, " %d %d", &a->lien[i] , &a->ordl[i]) ;
/* neighbour's index and bond order */
			if (v != 2) {
				return FALSE ;
			}
		}
		for (i = 0 ; i < 2 ; i++) {
			v = fscanf(stdin, "%lf", &a->coor[i] ) ;
			if (v != 1) {
				return FALSE ;
			}
/* 2D coordinates */
		}
	}
	return TRUE ;
/* so far so good */
}


void	wr2d(int n)
{
	int	i, j, k ;
	double	maxcoor, fact, x, y ;
	char	titre[5] ;

	sprintf(titre, "%d", isol) ;
	printf("%%%%Page: ? %d\n", isol) ;
	printf("bigfont setfont\n") ;
	printf("%7.1f%7.1f m (%s) show\n", 1.5 * SC, 1.85 * SC, titre) ;
/* x and y coordinates are in the range 0 - 2*SC */

	maxcoor = 0.0 ;
	for (i = 0 ; i < n ; i++) {
		x = fabs(at[i].coor[0]) ;
		y = fabs(at[i].coor[1]) ;
		if (x > maxcoor) {
			maxcoor = x ;
		}
		if (y > maxcoor) {
			maxcoor = y ;
		}
	}
	maxcoor *= 1.2 ;
	fact = SC / maxcoor ;
/* factor from atom coordinates to PS coordinates */
	for (i = 0 ; i < n ; i++) {
		at[i].coor[0] *= fact ;
		at[i].coor[1] *= fact ;
		at[i].coor[0] += SC ;
		at[i].coor[1] += SC ;
	}
/* computes PS coordinates */
	printf("smallfont setfont\n") ;
	for (i = 0 ; i < n ; i++) {
		for (k = 0 ; k < at[i].nlien ; k++) {
			j = at[i].lien[k] ;
			if (i < j) {
				segment(i, j, at[i].ordl[k]) ;
/* draws bonds */
			}
		}
	}
	for (i = 0 ; i < n ; i++) {
		numero(i) ;
/* draws atoms numbers */
	}
	if (ncharge) {
		printf("minifont setfont\n") ;
	}
	for (i = 0 ; i < n ; i++) {
		if(at[i].charge) {
			charge(i) ;
		}
	}
/* draws atoms charges */
	printf("showpage\n") ;
}


void	segment(int i, int j, int order)
/*
 * draws bonds 
 */
{
	struct at3D *a1, *a2 ;
	int	k ;
	double	x1, y1, x2, y2, dx, dy, dr, nx, ny ;
	double	p1[2], p2[2], p3[2] ;

	a1 = &at[i] ; 
	a2 = &at[j] ;
	x1 = a1->coor[0] ; 
	y1 = a1->coor[1] ;
/* x is the starting point */
	x2 = a2->coor[0] ; 
	y2 = a2->coor[1] ;
/* y is the end point */
	dx = (x2 - x1) ; 
	dy = (y2 - y1) ; 
/* the x-y segment */
	dr = hypot(dx, dy) ;
/* the x-y distance */
	nx = P2 * dy / dr * SC ; 
	ny = -(P2 * dx / dr) * SC ;
/* the normal to the x-y segment whose length is the distance between 
 * the segments forming a multiple bond */
	x1 -= (order - 1.0) / 2.0 * nx ; 
	y1 -= (order - 1.0) / 2.0 * ny ;
	x2 -= (order - 1.0) / 2.0 * nx ; 
	y2 -= (order - 1.0) / 2.0 * ny ;
/* x and y are moved in order to draw the first segment of a multiple bond */
	p1[0] = x1 ; 
	p1[1] = y1 ; 
	p2[0] = x2 ; 
	p2[1] = y2 ;

	for (k = 0 ; k < order ; k++) {
/* as many segments as the order of the bond */
		printf("n%7.1f%7.1f m%7.1f%7.1f l s\n",  p1[0], p1[1], p2[0], p2[1]) ;
/* draws it */
		p1[0] += nx ; 
		p1[1] += ny ;
		p2[0] += nx ; 
		p2[1] += ny ;
	}
}


void	numero(int i)
/*
 * draws the lsd reference number at the atom location
 */
{
	char	s[5] ;
	double	x, y ;

	sprintf(s, "%s", at[i].nom) ;
	sprintf(s + strlen(s), "%d", at[i].loc) ;
	x = at[i].coor[0] ;
	y = at[i].coor[1] - CORPS / 2;
	printf("%7.1f%7.1f m (%s) c\n", x, y, s) ;
}


void	charge(int i)
/*
 * draws the atom charge
 */
{
	char	s[3] ;
	double	x, y ;
	int c ;
	
	c = at[i].charge ;
	if(c > 0) {
/* positive charge */
		if(c == 1) {
			sprintf(s, "+") ;
		} else {
			sprintf(s, "%d+", c) ;
		}
	} else {
/* negative charge */
		if(c == -1) {
			sprintf(s, "-") ;
		} else {
			sprintf(s, "%d-", abs(c)) ;
		}
	}
	x = at[i].coor[0] + CORPS ;
	y = at[i].coor[1] + CORPS / 4 ;
/* sets coordinates of charge */
	printf("%7.1f%7.1f m (%s) c\n", x, y, s) ;
}


void	entete(void)
/*
 * Postscript header
 */
{
	printf("%%!PS-Adobe-3.0\n") ;
/* gsview in now able to forwards *and* backwards navigate in the file */
	printf("%%%%Creator: GENPOS 2.0\n") ;
	printf("%%%%Pages: (atend)\n") ;
	printf("%%%%EndComments\n") ;
	printf("%%%%BeginProlog\n") ;
	printf(" /l /lineto load def\n") ;
	printf(" /m /moveto load def\n") ;
	printf(" /s /stroke load def\n") ;
	printf(" /n /newpath load def\n") ;
	printf(" /c { dup stringwidth pop 2 div neg 0 rmoveto show } def\n") ;
	printf(" /bigfont /Times-Bold findfont 30 scalefont def\n") ;
	printf(" /smallfont /Times-Bold findfont %d scalefont def\n", CORPS) ;
	printf(" /minifont /Times-Bold findfont %d scalefont def\n", CORPS - 5) ;
	printf("%%%%EndProlog\n") ;

}


