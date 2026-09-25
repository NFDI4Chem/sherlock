/*
    file: OutlsdSrc/outlsd.c
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
 * outlsd p < file.sol
 *
 * p = 1: bond lists
 * p = 2: ARGOS format (obsolete)
 * p = 3: JMNMOL format (obsolete)
 * p = 4: Ye Olde MacroModel format (obsolete)
 * p = 5: SMILES chains
 * p = 6: LSD 2D DRAW format (.coo)
 * p = 7: SDF 2D format (.mol)
 * p = 8: SDF 3D format with H atoms (.sdf)
 * p = 9: SDF 3D format without H atoms (.sdf)
 * p = 10: SDF 0D format without H atoms (.sdf)
 */

#include "defs.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "matutil.h"
#include "jacobi.h"

#define SQR(X) ((X)*(X))

#define A -0.3
#define B 0.000001
#define C 0.2
#define E 0.5
#define THRS 1.0

struct atom at[MAX_AT] ; /* atoms */
int nl; /* number of atoms, invalid ones included */
FILE *fin ; /* input file */
int p ; /* command line parameter */
int na ; /* number of OULSD atoms (valid LSD atoms) */
int nawH ; /* number of OULSD atoms, including H */
int isol ; /* solution index */
char bsym[3] = {'-', '=', '#'} ; /* bond symbols */
int tt[MAX_TT][3] ; /* Triangle and Tetrahedron supplementary bonds */
int ntt ; /* number of tt bonds */
// int bonded [MAX_AT][MAX_AT] ; /* connectivity matrix */


void start(int argc, char *argv[]) ;
void pr_use() ;
void start_reading() ;
void start_writing() ;
int readsol() ;
void too_many_atoms() ;
void procsol() ;
void gen3D(int hastt) ;
void addtt() ;
double myrand (double range) ;
void minimiz(int nat, int fltt, int dim) ;
int bonded(int i, int j) ;
void gen3DwH() ;
void gen2D() ;
void smash() ;
void writesol() ;
void writelist() ;
extern void writesmiles() ;
void writedraw() ;
void writesdf(int wH) ;
void checkcoo(int nat) ;

int main(int argc, char *argv[])
{
	start(argc, argv) ;
	while(readsol()) {
		procsol() ;
		writesol() ;
	}
}

/*
 * void start(int argc, char *argv[])
 * read command line argument and set input stream ready for reading solutions
 */
void start(int argc, char *argv[])
{
	if(argc != 2) {
/* only one command line parameter */
		pr_use() ;
	} else {
		p = atoi(argv[1]) ;
		if(p < 1 || p > 10) {
			pr_use() ;
		}
		start_reading() ;
/* prepare input file for solution reading */
		start_writing() ;
/* write first line for some formats */
		srand(1) ;
/* prepares the random number generator for 2D/3D coordinate generation */
	}
}

/*
 * void pr_use() print usage message
 */
void pr_use()
{
	printf("outlsd: usage: outlsd p\n") ;
	printf("p = 1: bond lists\n") ;
	printf("p = 2: ARGOS format (obsolete format)\n") ;
	printf("p = 3: JMNMOL format (obsolete format)\n") ;
	printf("p = 4: Ye Olde MACROMODEL format (obsolete format)\n") ;
	printf("p = 5: SMILES chains\n") ;
	printf("p = 6: LSD 2D DRAW format (.coo)\n") ;
	printf("p = 7: SDF 2D format (.mol)\n") ;
	printf("p = 8: SDF 3D format without H atoms (.sdf)\n") ;
	printf("p = 9: SDF 3D format with H atoms (.sdf)\n") ;
	printf("p = 10: SDF 0D format without H atoms (.sdf)\n") ;
	exit(1) ;
}

/*
 * void start_reading() reads the comment part of the input file, if any,
 * and checks the presence of the OUTLSD header.
 */
void start_reading()
{
	int c ; /* current character */
	char header[7] ; /* storage for OUTLSD header string */
	
	c = fgetc(stdin) ;
	if(c == EOF) {
/* empty file */
		printf("outlsd: input file is empty.\n") ;
		pr_use() ;
	}
	if(c == '#') {
/* there is a comment, reads until end of comment */
		while( ((c = fgetc(stdin)) != EOF) && (c != '#') ) ;
		if(c == EOF) {
/* a second line that starts with a # is expected */
			printf("outlsd: comment has no end.\n") ;
			pr_use() ;
		}
		c = fgetc(stdin) ; 
/* skip to next line */
		if(c == EOF) {
			printf("outlsd: unexpected end-of-file.\n") ;
			pr_use() ;
		}
	} else {
/* no comment in input file */
		ungetc(c, stdin) ;
/* put back the character that was read to detect a comment */
	}
	scanf("%6s", header) ;
	if(strcmp(header, "OUTLSD")) {
		printf("outlsd: This is not a file for OUTLSD.\n") ;
		exit(1) ;
	}
}

/*
 * void start_writing() writes global header of files for JMNMOL
 * and for m_edit and genpos (.coo) input files
 */
void start_writing()
{
	switch(p) {
		case 3 : printf("JMNMOL\n") ; break ;
		case 6 : printf("DRAW\n") ; break ;
		default : break ;
	}
}

/*
 * void readsol()
 */
int readsol()
{
	int ia, ib, in, x, v, ok ;
	struct atom *a ;

	ia = 0 ;
/* na is the index of the current OUTLSD atom */
	v = scanf("%d %d", &nl, &isol) ;
/* first line, 2 ints: number of atom description lines, solution index */
	if(nl >= MAX_AT) {
/* the number of LSD atoms is limited to MAX_AT-1 */
		too_many_atoms() ;
	}
	if(!v || !nl) {
		return FALSE ;
	}
/* v == 0 means bad reading (final 0 missing), nl == 0 means correct end of solution set */
	for(x=1 ; x<=nl ; x++) {
/* x is the index of the current LSD atom */
		if(ia == MAX_AT) {
/* the number of OUTLSD atoms is limited to MAX_AT */
			too_many_atoms() ;
		}
		a = &at[ia] ;
		scanf(" %1d", &ok) ;
/* does LSD atom number nl exist? */
		if(ok == TRUE) {
/* a line that starts with 1 corresponds to an existing atom */
			v = scanf(" %2s %1d %1d %1d %1d %2d",
				a->elt,
				&a->valence,
				&a->mult,
				&a->nbmax,
				&a->nb,
				&a->charge ) ;
			if(v != 6) {
				return FALSE ;
			}
/* element name, valence, number of attached H, number of expected neighbors and of real ones, charge */
			a->savnb = a->nb ;
/* nb could be modified when implicit H are added, initial value in savnb */
			a->hybrid = a->valence - a->mult - a->nbmax ;
/* hybridization is 0 for sp3, 1 for sp2 and 2 for sp atoms */
			for(ib=0 ; ib<MAX_NN ; ib++) {
/* read neighbor index (LSD) and bond order */
				v = scanf(" %3d %1d",
					&a->b[ib] ,
					&a->bo[ib] ) ;
				if(v != 2) {
					return FALSE ;
				}
			}
			a->n2o = x ;
/* OUTLSD ("new") to LSD ("old") atom index conversion */
			at[x].o2n = ia ;
/* LSD ("old") to OUTLSD ("new") atom index conversion */
			ia++ ;
		}
	}
	na = ia ;
/* na is the true number of atoms in the molecule */
	for(ia = 0 ; ia < na ; ia++) {
/* scan each atom ... */
		a = &at[ia] ;
		for(ib = 0 ; ib < a->nb ; ib++) {
/* ... and scan each of its bonds ... */
			in = a->b[ib] = at[a->b[ib]].o2n ;
/* ... to replace LSD atoms indexes by OUTLSD ones. */
		}
	}
	return TRUE ;
}

/*
 * void too_many_atoms()
 */
void too_many_atoms()
{
	printf("outlsd: Too many atoms. Recompile with higher MAX_AT in defs.h.\n") ;
	exit(1) ;
}

/*
 * void procsol() is the general solution processor
 */
void procsol()
{
	switch(p) {
		case 8 :
		case 9 : gen3DwH() ; break ;
		case 6 :
		case 7 : gen2D() ; break ;
		default : break ;
	}
}

/*
 * void gen3D(int hastt) generates 3D coordinates without H
 */
void gen3D(int hastt)
{
	double edge ;
	int k, ia ;
	struct atom *a ;
	double *coo ;
	
	edge = exp(log((double)nl)/3.0) ;
/* the volume of a cube of edge lenth edge is na, the number of atomes */
	for(ia = 0 ; ia < na ; ia++) {
		coo = at[ia].coor ;
		for(k = 0 ; k < 3 ; k++) {
			coo[k] = myrand(edge) ;
/* initializes atom coordinates with random numbers in the -edge/2 .. edge/2 range */
		}
	}
	ntt = 0 ;
	if(hastt) {
		addtt() ;
	}
	minimiz(na, hastt, 3) ; 
}

/*
 * void addtt() creates tetrahedral, triangular (and linear) fictitious bons
 */
void addtt()
{
	int i, j, k, itt, aj, ak, ia, h, nn, sw ;
	struct atom *a ;
	
	itt = ntt ;
/* ready to add a new tt fictitious bond. ntt has been initialized at 0 */
	for(ia = 0 ; ia < na ; ia++) {
		a = &at[ia] ;
		h = a->hybrid ;
		nn = a->nb ;
		for(j = 0 ; j <nn - 1 ; j++) {
			for(k = j + 1 ; k < nn ; k++) {
/* consider all pairs (aj, ak) of neighbors of atom ia */
				aj = a->b[j] ;
				ak = a->b[k] ;
				if(aj > ak) {
/* set aj < ak */
					sw = aj ;
					aj = ak ;
					ak = sw ;
				}
				for(i = 0 ; i < itt && !(tt[i][0]==aj && tt[i][1]==ak) ; i++) ;
/* check that this bond is original */
				if(i == itt) {
/* new tt bond */
					if(itt==MAX_TT) {
						printf("outlsd: recompile with a higer value of MAX_TT\n") ;
						exit(1) ;
					}
/* ok to store the new tt bond */
					tt[itt][0] = aj ;
					tt[itt][1] = ak ;
/* tt bond extremities */
					tt[itt][2] = h ;
/* geometry  around atom ia, 0 for tetrahedric, 1 for triangular, 2 for linear */
					itt++ ;
/* ready for a new tt bond */
				}
			}
		}
	}
	ntt = itt ;
}

/*
 * double myrand (double range) returns a random number 
 * in the [-range/2, range/2] interval ;
 */
double myrand (double range)
{
	return ((double)rand()/(double)RAND_MAX - 0.5)*range ;
}

/*
 * void minimiz(int nat, int fltt, int dim) operates a geometry optimization
 * over atoms 0..(nat-1) in at. If fltt is true, the fictitious bonds
 * in tt are taken into account. if dim is 2, only the first two coordinates
 * in the coor field af atoms are considered.
 */
void minimiz(int nat, int fltt, int dim)
{
	struct atom *ai, *aj ;
	int i, j, k, n, fli, flj, a1, a2, itt, h ;
	double rms, v[3], D, f, r2, di ;
	static double Dtab[3] = {8.0/3.0, 3.0, 4.0} ;
	
// printf("Initial state\n") ;
// checkcoo(nat) ;
	rms = 1.0 ;
/* so that the loop below can be entered... */
	for(n = 0 ; n < MAX_ITER && rms > MAX_RMS ; n++) {
/* iterates until atoms do not move or allowed iteration number not reached */ 
		for(i = 0 ; i < nat ; i++) {
/* loop over atoms */
			ai = &at[i] ;
			for(k = 0 ; k < dim ; k++) {
				ai->dsp[k] = 0.0 ;
			}
/* reset atom displacement */
			for(j = 0 ; j < nat ; j++) {
				if(i == j) {
					continue ;
				}
/* all (i, j), i different of j, atom pairs are now considered */
				aj = &at[j] ;
				r2 = 0.0 ;
				for(k = 0 ; k < dim ; k++) {
					v[k] = aj->coor[k] - ai->coor[k] ;
					r2 += SQR(v[k]) ;
				}
//				r2 = SQR(v[0]) + SQR(v[1]) + SQR(v[2]) ;
/* r2 = SQR(d(i,j)) with Euclidian distance, v is the interatomic vector */
				if(bonded(i, j)) {
					fli = ai->elt[0] == 'H' ;
					flj = aj->elt[0] == 'H' ;
					D = (fli || flj) ? 0.25 : 1.0 ;
/* D is the square of the equilibium distance of the i-j bond */ 
					f = C*(1-1/(r2/D + B)) ;
/* f is 0 when r2 is equal to D */
				} else {
					f = A/SQR(r2+B) ;
/* repulsive term */
				}
				for(k = 0 ; k < dim ; k++) {
					ai->dsp[k] += f * v[k] ;
				}
/* accumulate displacement for i according to bonded/non-bonded interactions */
			}
		}
		if(fltt) {
/* incorporate pseudo-bonded interactions that define the angles around atoms */
			for(itt = 0 ; itt < ntt ; itt++) {
				a1 = tt[itt][0] ;
				a2 = tt[itt][1] ;
/* consider the a1-a2 pseudo-bond */
				r2 = 0.0 ;
				for(k = 0 ; k < dim ; k++) {
					v[k] = at[a2].coor[k] - at[a1].coor[k] ;
					r2 += SQR(v[k]) ;
				}
/* r2 = SQR(d(a1,a2)) */
				h = tt[itt][2] ;
				D = Dtab[h] ;
/* D is the square of the a1-a2 equilibrium distance */
				f = E*(1-1/(B + r2/D)) ;
				for(k = 0 ; k < dim ; k++) {
					di = v[k] * f ;
					at[a2].dsp[k] -= di ;
					at[a1].dsp[k] += di ;
/* displacements due to pseudo-bonds are taken into account */
				}
			}
		}
		rms = 0.0 ;
/* reset sum of squares of displacements */
		for(i = 0 ; i < nat ; i++) {
			for(k = 0 ; k < dim ; k++) {
				rms += SQR(at[i].dsp[k]) ;
			}
		}
		rms = sqrt(rms/nat) ;
/* rms displacement */
		f = (rms > THRS) ? THRS/rms : 1.0 ;
/* displacement scaling factor */
		for(i = 0 ; i < nat ; i++) {
			for(k = 0 ; k < dim ; k++) {
				at[i].coor[k] += f * (at[i].dsp[k]) ;
			}
		}
// printf("after iteration %d  rms: %25.20f\n", n, rms) ;
// checkcoo(nat) ;
/* coordinates update */
	}
/*	for(i = 0 ; i < itt ; i++) {
		a1 = tt[i][0] ;
		a2 = tt[i][1] ;
		for(k = 0 ; k < dim ; k++)
			v[k] = at[a2].coor[k] - at[a1].coor[k] ;
		r2 = sqrt(car(v[0]) + car(v[1]) + car(v[2])) ;
		printf("%d %d %f\n", a1+1, a2+1, r2) ;
	}
*/
}

/*
 * true if atoms i and j are bonded
 */
int bonded(int i, int j)
{
	struct atom *a ;
	int ib, nn ;

	a = &at[i] ;
	nn = a->nb ;
	for(ib = 0 ; ib < nn && a->b[ib] != j ; ib++) ;
	return ib != nn ;
}


/*
 * void gen3DwH() generates 3D coordinates with H
 */
void gen3DwH()
{
	int ia, nh, ih, k, nn ;
	struct atom *a, *aH ;
	
	gen3D(TRUE) ;
// checkcoo(na) ;
/* start with coordinate generation without H */
	nawH = na ;
/* nawH: number of atoims, with H included */
	for(ia = 0 ; ia < na ; ia++) {
		a = &at[ia] ;
		nh = a->mult ;
/* number of H to be added */
		for(ih = 0 ; ih < nh ; ih++) {
			if(nawH == MAX_AT) {
				printf("Too many atoms. Recompile with a higer MAX_AT.\n") ;
				exit(1) ;
			}
/* exit if too many atoms */
			aH = &at[nawH] ;
			aH->elt[0] = 'H' ;
			aH->elt[1] = '\0' ;
			aH->nb = 1 ;
			for(k = 0 ; k < 3 ; k++) {
				aH->coor[k] = a->coor[k] + myrand(1.0) ;
/* H is randomly placed around its parent atom */
			}
			aH->b[0] = ia ;
			aH->bo[0] = 1 ;
			for(k = 1 ; k < 4 ; k++) {
				aH->b[k] = 0 ;
			}
/* create the new H atom */
			nn = a->nb ;
			a->b[nn] = nawH ;
			a->bo[nn] = 1 ;
			a->nb++ ;
/* update atom ia, with its new H atom */
			nawH++ ;
/* ready for the next added H atom */
		}
	}
// checkcoo(nawH) ;
	minimiz(nawH, FALSE, 3) ;
/* optimize 3D coordinates of all atoms but without TT bonds */
}

/*
 * void gen2D() generates 2D coordinates without H
 */
void gen2D()
{
	gen3D(FALSE) ;
/* start with coordinate generation without H */
	smash() ;
/* make flat */
	minimiz(na, FALSE, 2) ;
/* optimize 2D coordinates without TT bonds */
}

/*
 * void smash() generates 2D coordinates from 3D coordinates using inertia matrix
 */
void smash()
{
	double gr[3], x, y, z ;
	double max2, max3, xn, yn ;
	int ia, k, l, iax1, iax2, iax3 ;
	double **iner, *mom, **evec, *coo ;
	unsigned int nrot ;

	iner = dmatrix(3, 3) ;
/* inertia matrix*/
	mom = dvector(3) ;
/* inertia moments */
	evec = dmatrix(3, 3) ;
/* inertia principal axes */
	if((!iner) || (!mom) || (!evec)) {
		printf("%s", "memory allocation failed.\n") ;
		exit(1) ;
	}
/* check memory allocation */
	for(k = 0 ; k < 3 ; k++) {
		gr[k] = 0.0 ;
	}
/* reset position of mass center (all atoms the same mass!) */
	for(ia = 0 ; ia < na ; ia++) {
		coo = at[ia].coor ;
		for(k = 0 ; k < 3 ; k++) {
			gr[k] += coo[k] ;
		}
	}
/* atomwise coordinate summation */
	for(ia = 0 ; ia < na ; ia++) {
		coo = at[ia].coor ;
		for(k = 0 ; k < 3 ; k++) {
			coo[k] -= gr[k]/na ;
		}
	}
/* molecule coordinate recentering */

	for(k = 0 ; k < 3 ; k ++) {
		for(l = 0 ; l < 3 ; l++) {
			iner[k][l] = 0 ;
		}
	}
/* reset inertia matrix */
	for(ia = 0 ; ia < na ; ia++) {
	 	coo = at[ia].coor ;
		x = coo[0] ;
		y = coo[1] ;
		z = coo[2] ;
		iner[0][0] += SQR(y) + SQR(z) ;
		iner[1][1] += SQR(z) + SQR(x) ;
		iner[2][2] += SQR(x) + SQR(y) ;
		iner[0][1] -= x * y ;
		iner[1][2] -= y * z ;
		iner[2][0] -= z * x ;
	}
/* upper part of inertia matrix */
	iner[1][0] = iner[0][1] ;
	iner[2][1] = iner[1][2] ;
	iner[0][2] = iner[2][0] ;
/* inertia matrix is symmetric */
/*
for(i=0 ; i<3 ; i++) {
	for(j=0 ; j<3 ; j++) {
		printf("%.6f ", iner[i][j]) ;
	}
	printf("%s", "\n") ;
}
*/

	switch (gsl_eigen_jacobi(3, iner, mom, evec, 1000, &nrot))
	{
		case 0 : break ;
		case 1 : printf("%s", "cannot converge.\n") ; exit(1) ;
		case 2 : printf("%s", "cannot allocate workspace.\n") ; exit(1) ;
	}
/* diagonalization of inertia matrix */


	max3 = 0.0 ;
/* looking for the axis of highest inertia moment */
	for(k = 0 ; k < 3 ; k++) {
		if(mom[k] > max3) {
			max3 = mom[k] ;
			iax3 = k ;
		}
	}
/* atom projections have to be made on the plane than is orthogonal to axis iax3 */
	max2 = 0.0 ;
/* looking for the axis of highest inertia moment that is orthogonal to iax3 */
	for(k = 0 ; k < 3 ; k++) {
		if(k != iax3 && mom[k] > max2) {
			max2 = mom[k] ;	
			iax2 = k ;
		}
	}
/* iax2 will be the y axis of the projection plane */
	iax1 = 3 - iax2 - iax3 ;
/* iax1 will be the x axis of the projection plane */

	for(ia = 0 ; ia < na ; ia++) {
	 	coo = at[ia].coor ;
		x = coo[0] ;
		y = coo[1] ;
		z = coo[2] ;
		xn = x*evec[0][iax1] + y*evec[1][iax1] + z*evec[2][iax1] ;
		yn = x*evec[0][iax2] + y*evec[1][iax2] + z*evec[2][iax2] ;
		coo[0] = xn ;
		coo[1] = yn ;
		coo[2] = 0 ;
	}
/* projection of atom 3D coordiantes in the xy plane */

	free_dmatrix(evec) ;
	free_dvector(mom) ;
	free_dmatrix(iner) ;
/* free dynamically allocated memory */
}

/*
 * void writesol() is the general solution writer
 */
void writesol()
{
	switch(p) {
		case 1 : writelist() ; break ;
		case 2 : printf("ARGOS format is obsolete") ; break ;
		case 3 : printf("JMNMOL format is obsolete") ; break ;
		case 4 : printf("MacroModel format is obsolete") ; break ;
		case 5 : writesmiles() ; break ;
		case 6 : writedraw() ; break ;
		case 7 :
		case 8 : writesdf(FALSE) ; break ;
		case 9 : writesdf(TRUE) ; break ;
		case 10 : writesdf(FALSE) ; break ;
		default : break ;
	}
}

/*
 * void writelist() writes molecule as a list of bonds
 */
void writelist()
{
	int ia, ib, in, nn, x, y ;
	struct atom *a ;

	printf("\n*** solution %d *** \n", isol) ;
	for(ia = 0 ; ia < na ; ia++) {
/* loop over bonds */
		a = &at[ia] ;
/* pointer to current atom */
		x = a->n2o ;
/* LSD number of current atom */
		nn = a->nb ;
/* number of neighbors of the current atom */
		for(ib = 0 ; ib < nn ; ib++) {
/* loop over bonds of current atom */
			in = a->b[ib] ;
/* in is the ib'th neighbor of ia */
			y = at[in].n2o ;
/* LSD number of current neighbor */
			if(x < y) {
/* do not duplicate bonds */
				printf("%4d%4d %c\n", x, y, bsym[a->bo[ib] - 1]) ;
/* print bond with symbol */
			}
		}
	}
}

/*
 * void writedraw() writes 2D molecules in the LSD draw format
 */
void writedraw()
{
	int ia, ib, k, nn ;
	struct atom *a ;
	double *coo ;

	printf("%3d  %3d\n", na, isol) ;
/* first line of atom description: number of atoms, solution index */
	for(ia = 0 ; ia < na ; ia++) {
		a = &at[ia] ;
		nn = a->nb ;
		coo = a->coor ;
		printf("%3d%3d%3s%3d%2d  ", ia, a->n2o, a->elt, a->charge, nn) ;
		for(ib = 0 ; ib < nn ; ib++) {
			printf("%3d%2d", a->b[ib], a->bo[ib]) ;
		}
		for(ib = nn ; ib < 6 ; ib++) {
			printf("%3d%2d", 0, 0) ;
		}
		for(k = 0 ; k < 2 ; k++) {
			printf("%10.5f ", coo[k]) ;
		}
		printf("\n") ;
	}
}

/*
 * void writesdf(int wH) writes 2D molecules in SDF format, with H if wH is true.
 */
void writesdf(int wH)
{
	int nat, ia, nbond, k, ib, nn, in, ncharge, c ;
	double *coo ;
	struct atom *a ;

	nat = wH ? nawH : na ;
/* number of atoms to consider */
	printf("\n     OUTLSD version 01/26/2022 \n\n");
/* MOL/SDF header */
	nbond = 0 ;
	ncharge = 0 ;
	for(ia = 0 ; ia < nat ; ia++) {
		nbond += wH ? at[ia].nb : at[ia].savnb ;
		if(at[ia].charge) {
			ncharge++ ;
		}
	}
	nbond /= 2 ;
/* bond number */
/*	printf("%3d%3d  0  0  0  0  0  0  0  0%3d V2000\n", nat, nbond, ncharge ? 2 : 1) ; */
	printf("%3d%3d  0  0  0  0  0  0  0  0999 V2000\n", nat, nbond) ;
/* molecule header */
	for(ia = 0 ; ia < nat ; ia++) {
		a = &at[ia] ;
		coo = a->coor ;
		for(k = 0 ; k < 3 ; k++) {
			printf("%10.4f", coo[k]) ;
		}
		printf(" %-3s 0%3d  0  0  0%3d  0  0  0  0  0  0\n", a->elt, a->charge ? 4 - a->charge : 0, ( (((!strcmp(a->elt, "N")) || (!strcmp(a->elt, "P"))) && (a->valence == 5)) || ((!strcmp(a->elt, "S")) && ((a->valence == 4) || (a->valence == 6))) ) ? a->valence : 0) ;
	}
/* write atom block */
	for(ia = 0 ; ia < nat ; ia++) {
		a = &at[ia] ;
		nn = wH ? at[ia].nb : at[ia].savnb ;
		for(ib = 0 ; ib < nn ; ib++) {
			in = a->b[ib] ;
			if(ia < in) {
				printf("%3d%3d%3d  0  0  0  0\n", ia + 1, in + 1, a->bo[ib]) ;
			}
		}
	}
/* write bond block */
	if(ncharge) {
		printf("M  CHG%3d", ncharge) ;
		for(ia = 0 ; ia < nat ; ia++) {
			if(c = at[ia].charge) {
				printf("%4d%4d", ia + 1, c) ;
			}
		}
		printf("\n") ;
	}
	printf("M  END\n") ;
	printf("$$$$\n") ;
/* end of molecule */
}

void checkcoo(int nat)
{
	double *coo ;
	int ia, k ;

	for(ia = 0 ; ia < nat ; ia++) {
		coo = at[ia].coor ;
		printf("%2d:", ia) ;
		for(k = 0 ; k < 3 ; k++) {
			printf("%25.20f", coo[k]) ; 
		}
		printf("\n") ; fflush(stdout) ;
	}
	printf("\n") ;	fflush(stdout) ;
}

