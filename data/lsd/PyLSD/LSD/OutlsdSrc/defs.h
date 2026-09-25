/*
    file: OutlsdSrc/defs.h
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

#define MAX_AT 200 
#define MAX_ITER 80
#define MAX_RMS 0.001
#define MAX_TT 400
#define MAX_NN 4

#define X 1000
#define MAX_LEN 300
#define STRL 16
#define MAX_NB 300
#define MAX_CU MAX_AT

#define FALSE 0
#define TRUE 1

struct atom {
	char elt[3] ;
	int valence ;
	int mult ;
	int hybrid ;
	int charge ;
	int nbmax ;
	int nb ;
	int b[MAX_NB] ;
	int bo[MAX_NB] ;
	int n2o ;
	int o2n ;
	int savnb ;
	int typx ;
	double coor[3] ;
	double dsp[3] ;	
} ;

