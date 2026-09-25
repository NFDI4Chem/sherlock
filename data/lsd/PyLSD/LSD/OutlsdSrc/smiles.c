/*
    file: OutlsdSrc/smiles.c
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
#define DEBUG 0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"

/*
 * SMILES chain generation
 * Each connected unit (CU) in the current solution is separately treated.
 */
 
extern struct atom at[MAX_AT] ; /* atoms */
extern int na ; /* number of OULSD atoms */
extern char bsym[3] ; /* bond symbols */

int cuna ;	/* number of atoms in a CU */
int cunb ;	/* number of bonds in a CU */
char buf[MAX_LEN] ;	/* buffer in which the SMILES chain is created */
int ibuf ;	/* buffer index */
int ue1[MAX_NB], ue2[MAX_NB] ;	/* atom indexes that define the bonds in a CU */
int ubtyp[MAX_NB] ;	/* bond order of the bonds in a CU */
int unn[MAX_AT] ;	/* number of neighbors of the atoms in a CU	*/
char label[MAX_AT][STRL] ;/* atom identifiers in a CU */
int d1[MAX_AT][MAX_AT] ;/* distance matrix */
int o_to_n[MAX_AT], n_to_o[MAX_AT] ; 
	/* conversions between OUTLSD atom numbers and CU atom numbers */
int sel[MAX_AT] ;	/* TRUE for already processed CU atoms */

void writesmiles() ;
void chain(int lev, int i, int j, int d, int type) ;
void w_char(char c) ;
void w_bond(int t) ;
void w_atom(int a) ;
int get_typ(int a1, int a2) ;

/*
 * void writesmiles() cuts the solution into connected units.
 * The ring closure bonds are searched for. Each CU is then analyzed in the order
 * of decreasing number of constituting atoms.
 * The distance matrix is calculated. The most distant atoms are selected
 * and used by chain() to recursively cut the molecule into non-branched pieces
 */

void writesmiles() {

	int nb ; 		/* number of bonds */
	int ia, in, ib, i, j, k, x ; 	/* indexes, current atoms */
	int btyp[MAX_NB], e1[MAX_NB], e2[MAX_NB], st[MAX_NB] ; 
			/* bond definition: bond type, extremities 1 and 2, closure status */
	int cua[MAX_AT], cub[MAX_NB] ; 	
			/* connected unit index for atoms and bonds */
	int t ; 		/* bond type, then current index */
	int nc1, nc2, nc3 ; 	/* CU numbers */
	int ncu, icu, cu, indxcu ;	
			/* nombre, indice, no. d'uc 			*/
	int cuindx[MAX_CU], cupop[MAX_CU] ;
		 	/* UC index and population */
	int d2[MAX_AT][MAX_AT] ;
		 	/* temporary distance matrix */
	int ust[MAX_NB] ; 	/* FALSE if closure bond in current CU */
	char *a ;		/* current string*/
	int a1, a2 ;		/* current atoms that define a current bond */
	int l1, l2 ;		/* used in distance calculation*/
	int im, jm, dm, dij ;	/* used in most distant atoms search */
	int c ;	/* atom charge */
	char bufstr[STRL] ;	/* string buffer */
#if DEBUG==1
	int z, zz ;		/* for debugging */
#endif

	ibuf = 0 ;
/* start writing SMILES chain */

	for(ia = 0 ; ia < na ; ia++) {
		if(at[ia].nb != at[ia].nbmax) {
			printf("outlsd: option 5: each atom must have all its neighbors\n") ;
			exit(1) ;
		}
	}
/* check that all atoms are complete */

	nb = 0 ;
	for(ia = 0 ; ia < na ; ia++) {
		nb += at[ia].nb ;
	}
	nb /= 2 ;
/* nb bonds found*/
	if(nb > MAX_NB) {
		printf("outlsd: too many bonds, recompile with higher MAX_NB.\n") ;
		exit(1) ;
	}
/* check limit on number of bonds */

#if DEBUG==1
	printf("%d atoms and %d bonds.\n", na, nb) ; 
#endif

	i = 0 ;
/* current bond index */
	for (t = 3 ; t > 0 ; t --) {
/* triple, then double, then single bonds */
		for(ia = 0 ; ia < (na-1) ; ia++) {
/* ia is the first atom in bond */
			for(in = 0 ; in < at[ia].nb ; in++) {
				ib = at[ia].b[in] ;
#if DEBUG==1
	printf("type %d:  ia: %2d  ib: %2d (%d)\n", t, ia, ib, at[ia].bo[in]) ;
#endif
/* loop through neighbors of current atom to find second atoms */
				if((ia < ib) && (at[ia].bo[in] == t)) {
#if DEBUG==1
	printf("OK i: %2d\n", i) ;
#endif
					e1[i] = ia ;
					e2[i] = ib ;
					btyp[i] = t ;
					i++ ;
/* new bond with extremities e1 and e2, and bond order btyp */
				}
			}
		}
	}
#if DEBUG==1	
	for(i = 0 ; i < nb ; i++) {
		printf("bond %d: %d - %d (type %d)\n", i+1, at[e1[i]].n2o, at[e2[i]].n2o, btyp[i]) ;
	}
#endif		

	for(ia = 0 ; ia < na ; ia++) {
		cua[ia] = ia + 1 ;
/* each atom gets a CU index equal to its index + 1(from 1 to na) */
	}
	nc3 = na + 1 ;
/* nc3 is the currently available CU index */
	for(i = 0 ; i < nb ; i++) {
		st[i] = FALSE ;
/* bond is a ring closure bond */
	}
	for(ib = 0 ; ib < nb ; ib++) {
/* loop through bonds */
		nc1 = cua[e1[ib]] ;
		nc2 = cua[e2[ib]] ;
/* CU indexes of the atoms of the current bond */
		if(nc1 == nc2) {
			continue ;
		}
/* if CU indexes are equal, then the current bond binds two atoms that already
   belong to the same CU and this bond is therefore a ring closure bond */
		for(ia = 0 ; ia < na ; ia++) {
/* loop through atoms */
			if(cua[ia] == nc1 || cua[ia] == nc2) {
				cua[ia] = nc3 ;
/* atoms e1[ib] and e2[ib] get nc3 as new CU index
   as well as and all atoms that have the same CU index as them  */
			}
		}
		nc3++ ;
/* next available CU index */
		st[ib] = TRUE ; 
/* this bond is not a ring closure */
	}
/* all atoms in the same CU have now the same cua */

	ncu = 0 ;
/* opposite of the number of CUs */
	for(;;) {
		for(ia = 0 ; ia < na && cua[ia] < 0 ; ia++) ;
/* loop over atoms stops when CU index is positive */
		if(ia == na) {
			break ;
/* if all atom CU index are negative, its finished */
		}
		cu = cua[ia] ;
/* cu is an index of still not reassigned UC index */
		ncu-- ;
/* new number of UC */
		for(ia = 0 ; ia < na ; ia++) {
			if(cua[ia] == cu) {
/* all atom with old CU index cu ... */
				cua[ia] = ncu ;
/* ... get ncu as new CU index */
			}
		}
	}
	ncu = -ncu ;
/* true number of CU */
	for(ia = 0 ; ia < na ; ia++) {
		cua[ia] = -cua[ia]-1 ;
/* UC index become positive, and 0 for the first one */
	}
    for(ib = 0 ; ib < nb ; ib++) {
    	cub[ib] = cua[e1[ib]] ;
/* the bonds also get a CU index, equal to the one of their extremities */
    }

#if DEBUG==1
	printf("BONDS\n") ;
	for(z = 0 ; z < nb ; z++) {
		printf("%4d(CU%2d)%c%4d(CU%2d) %s CLOSURE\n",
			at[e1[z]].n2o, cua[e1[z]], bsym[btyp[z]-1], at[e2[z]].n2o, cua[e2[z]],
			st[z] ? "NOT" : "") ;
	}
	printf("%d connected unit%s\n", ncu, (ncu>1) ? "s" : "") ;
#endif

	for(icu=0 ; icu < ncu ; icu++) {
		cupop[icu] = 0 ;
/* initial population of each CU is 0 */
	}
	for(icu = 0 ; icu < ncu ; icu++) {
		for(ia = 0 ; ia < na ; ia++) {
			if(cua[ia] == icu) {
/* ia is a new atom in CU icu */
				cupop[icu]++ ;
			}
		}
		cuindx[icu] = icu ;
/* initial index of current CU, for their sorting by order of decreasing population */
	}

	for(i = 0 ; i < ncu-1 ; i++) {
		for(j = i+1 ; j < ncu ; j++) {
			if(cupop[i] < cupop[j]) {
/* UC are sorted by decreasing population */
				t = cupop[i] ;
				cupop[i] = cupop[j] ;
				cupop[j] = t ;
/* population swap */
				t = cuindx[i] ;
				cuindx[i] = cuindx[j] ;
				cuindx[j] = t ;
/* CU index swap */
			}
		}
	}
	
#if DEBUG==1
	for(z = 0 ; z < ncu ; z++) {
		printf("CU %4d : index %4d : %4d atomes\n", z, cuindx[z], cupop[z]) ;
	}
#endif

	for(icu = 0 ; icu <ncu ; icu++) {
/* global CU processing */
		if(icu != 0) {
			w_char('.') ;
/* . is the CU separator in SMILES */
		}
		indxcu = cuindx[icu] ;
/* index of current CU, the biggest is processed first */
		cuna = cupop[icu] ;
/* number of atomes in the current CU */
		cunb = 0 ;
/* number of bonds in the current CU */
		for(ib = 0 ; ib < nb ; ib++) {
			if(cub[ib] == indxcu) {
/* bond ib has indxcu as CU index */
				cunb++ ;
/* one more bond in current CU */
			}
		}
		for(i = 0 ; i < cuna ; i++) {
			sel[i] = FALSE ;
		}
/* no atom selected so far */

		j = 0 ;
/* j is the index of the currently available atom number in current CU, "new" number */
		for(i=0 ; i < na ; i++) {
			if(cua[i] == indxcu) {
/* i is the atom index in the solution, "old" number */
				o_to_n[i] = j ;
				n_to_o[j] = i ;
/* old to now and new to old conversion */
				unn[j] = at[i].nb ;
/* number of neibors of CU atom j */
				j++ ;
/* ready for the next CU atom */
			}
		}


		j = 0 ;
/* j ist now the index of the currently available bond in current CU */
		for(i=0 ; i<nb ; i++) {
			if(cub[i] == indxcu) {
/* bond i belongs to the current CU */
				ue1[j] = o_to_n[e1[i]] ;
				ue2[j] = o_to_n[e2[i]] ;
/* new style bond */
				ust[j] = 1 - st[i] ;
/* bond status is now TRUE for ring closures */
				ubtyp[j] = btyp[i] ;
/* new style bond order */
				j++ ;
/* ready for the new CU bond */
			}
		}

		for(i = 0 ; i < cuna ; i++) {
			if(c = at[n_to_o[i]].charge) {
				sprintf(label[i], "[%s", at[n_to_o[i]].elt) ;
				if(at[n_to_o[i]].mult) {
					strcat(label[i], "H") ;
					if(at[n_to_o[i]].mult > 1) {
						sprintf(bufstr, "%d", at[n_to_o[i]].mult) ;
						strcat(label[i], bufstr) ;
					}
				}
				if(c < -1) {
					sprintf(bufstr, "%d", c) ;
					strcat(label[i], bufstr) ;
				} else if(c == -1) {
					strcat(label[i], "-") ;
				} else {
					strcat(label[i], "+") ;
					if(c > 1) {
						sprintf(bufstr, "%d", c) ;
						strcat(label[i], bufstr) ;
					}
				}
				strcat(label[i], "]") ;
/* atom label with charge */
			} else if( (at[n_to_o[i]].mult >= 2) && (((at[n_to_o[i]].valence == 5) && ((!strcmp(at[n_to_o[i]].elt, "N")) || (!strcmp(at[n_to_o[i]].elt, "P")))) || ((!strcmp(at[n_to_o[i]].elt, "S")) && ((at[n_to_o[i]].valence == 4) || (at[n_to_o[i]].valence == 6)))) ) {
				sprintf(label[i], "[%sH%d]", at[n_to_o[i]].elt, at[n_to_o[i]].mult) ;
			} else if(!strcmp(at[n_to_o[i]].elt, "Si")) {
				sprintf(label[i], "[%s]", at[n_to_o[i]].elt) ;
/* writing "Si" as "[Si]" */
			} else {
				strcpy(label[i], at[n_to_o[i]].elt) ;
			}         
/* atom label in SMILES start with the atomic symbol */
		}

		j = 1 ;
/* j is now the index of the next closure bond */
		for(i=0 ; i<cunb ; i++) {
			if(ust[i]) {
/* bond i is a ring closure bond */
				a = label[ue1[i]] ;
				a += strlen(a) ;
/* where to append the closure index for the label of the first atom of bond */
				if(j > 9) {
/* in this case the bond closure index starts with a % */
					sprintf(a, "%c", '%') ; a++ ;
				}
				sprintf(a, "%d", j) ;
/* writes the ring closure index */
				a = label[ue2[i]] ;
/* the same for the second atom of the current ring closure */
				a += strlen(a) ;
				if(j > 9) {
					sprintf(a, "%c", '%') ; a++ ;
				}
				sprintf(a, "%d", j) ;
/* writes the ring closure index */
				unn[ue1[i]]-- ;
				unn[ue2[i]]-- ;
/* the ring closure bond extremities get one less neighbor in the acyclic graph */
				j++ ;
/* ready for the next closure */
			}
		}

#if DEBUG==1
		for(z = 0 ; z < cuna ; z++) {
			printf ("atom %3d : old %3d : element %s : %d neighbors\n",
				z, at[n_to_o[z]].n2o, label[z], unn[z]) ;
		}
#endif

/* d1 est la matrice des distances */

		for(i = 0 ; i < cuna ; i++) {
			for(j = 0 ; j < cuna ; j++) {
				d1[i][j] = X ;
/* distance matrix initialization: all atoms are very far away */
			}
		}
		for(i = 0 ; i < cunb ; i++) {
			if(!ust[i]) {
/* non-closure bond i is now considered */
				a1 = ue1[i] ;
				a2 = ue2[i] ;
				d1[a1][a2] = d1[a2][a1] = 1 ;
/* distance of bonded atoms is 1 */
			}
		}
		for(k = 0 ; k < cuna ; k++) {
			for(i = 0 ; i < cuna ; i++) {
				for(j = 0 ; j < cuna ; j++) {
					l1 = d1[i][j] ;
					l2 = d1[i][k] + d1[k][j] ;
					d2[i][j] = (l1 < l2) ? l1 : l2 ;
				}
			}
/* d2 is a temporary distance matrix, updated at step k */
			for(i = 0 ; i < cuna ; i++) {
				for(j = 0 ; j < cuna ; j++) {
					d1[i][j] = d2[i][j] ;
				}
			}
/* copy d2 back in d1 */
		}
		for(i = 0 ; i < cuna ; i++) {
			d1[i][i] = 0 ;
/* fill diagonal with zeros */
		}

#if DEBUG==1
		for(z = 0 ; z < cuna ; z++) {
			for(zz = 0 ; zz < cuna ; zz++) {
				printf("%3d", d1[z][zz]) ;
			}
			printf("\n") ;
		}
#endif

		dm = im = jm = 0 ;
/* dm will be the longest distance between 2 atoms, ij and jm */
		for(i = 0 ; i < cuna-1 ; i++) if(unn[i] == 1) {
			for(j = i+1 ; j < cuna ; j++) if(unn[j] == 1) {
				dij = d1[i][j] ;
				if(dij > dm) {
/* a new maximum distance is found ... */
					dm = dij ;
					im = i ;
					jm = j ; 
/* remember where */
				}
			}
		}

		chain(0, im, jm, dm, 1) ;
/* call chain at ground level (0) and without attachment sign (1, as for a single bond) */
	} /* for iuc */
	printf("%s\n", buf) ;
} /* writesmiles */

/*
 * void chain(int lev, int i, int j, int d, int type)
 * creates a SMILES chain in buffer buf. 
 * lev is the recursive call level, i and j are the most distant (d) atoms of the chain.
 * a chain is connected (level differebt of 0) to its patent chain by a bond of type (or order) t.
 */
void chain(int lev, int i, int j, int d, int type) {
int k ;			/* atom on the main chain between i and j */
int x ;			/* so that x-1 is d(i, k) */
int kp ;		/* so that d(i, kp) = x	*/
int I ;			/* possible branching point at position k */
int t ; 		/* bond type (order) of the I-k bond */
int Jm ;		/* side chain extremity */
int Dm ;		/* Dm = d(I, Jm) */
int J ;			/* current side chain extremity	*/
int D ;			/* current D, for the search of Jm and Dm */

/*
 *
 * i---*---*---k---kp---*---*---*---j
 *             |
 *             I
 *             |
 *             *
 *             |
 *             Jm
 *
 * There are d=8 bonds between main chain extremity atoms i and j.
 * There is a branching point at k with x=4 (d(i,kp)).
 * kp immediately follows k when going from i to j.
 * There are Dm=2 bonds between the attach point I of the side chain and the Jm extremity.
 *
 */
#if DEBUG==1
	printf("Level %2d: the longest chain goes from %d to %d and has %d bonds\n",
		lev, at[n_to_o[i]].n2o, at[n_to_o[j]].n2o, d) ;
#endif
	if(lev) {
		w_char('(') ;
		w_bond(type) ;
	}
/* bind the current chain with its parent, if any (level being not 0) */
	w_atom(i) ;
/* writes first atom in chain */
#if DEBUG==1
	printf("***   %s\n", buf) ;
#endif
	if(d) { /* if this chain is not made of a single atom */
		k = i ;
/* k will be the index of intermediate atoms between i and j along the longest chain */
		for(x=1 ; x<=d ; x++) {
/* loop over the distance betwwen i an kp. d(i,k) is x-1 */
			sel[k] = TRUE ;
/* atom k has been visited */
			for(kp=0 ; !(d1[i][kp] == x && d1[kp][j] == d-x) ; kp++) ;
/* d(i,kp) is x and d(i,kp)+d(kp,j) is d. atom k is possibly a branching point */
#if DEBUG==1
			printf("between %d and %d : %d\n", at[n_to_o[i]].n2o, at[n_to_o[j]].n2o, at[n_to_o[k]].n2o) ;
#endif
			if(unn[k] > 2) {
/* k is a branching point if it has two neighbors along the main chain and others for branching */
				for(I = 0 ; I < cuna ; I++) {
/* I is any atom in the molecule ... */
#if DEBUG==1
					printf("trying I = %d \n", at[n_to_o[I]].n2o) ;
#endif
					if(d1[I][k] != 1) {
						continue ;
					}
/* ... that is a neighbor of k ... */
#if DEBUG
					printf(" %d is a neighbor of %d\n", at[n_to_o[I]].n2o, at[n_to_o[k]].n2o) ;
#endif
					if(I == kp) {
						continue ;
					}
/* ... that is not in the main chain in the direction of j ... */
#if DEBUG == 1
					printf(" %d different of kp : %d\n", at[n_to_o[I]].n2o, at[n_to_o[kp]].n2o) ;
#endif
					if(sel[I]) {
						continue ;
					}
/* ... that is not in the main chain in the direction of i (non selected). */
#if DEBUG==1
					printf(" %d not already visited\n", at[n_to_o[I]].n2o) ;
#endif
					t = get_typ(I, k) ;
/* bond type (order) for the link between the main chain and the side chain */
					if(unn[I] == 1) {
/* I is only bonded to k */
						chain(lev + 1, I, I, 0, t) ;
/* print that chain of length 0 */
					} else {
/* look for Jm, an atom that is as far as possible from I, at distance Dm of I */
						Dm = 0 ;
						for(J = 0 ; J < cuna ; J++) {
/* J is any atom in the CU ... */
							if(unn[J] != 1) {
								continue ;
							}
/* ... that is a terminal atom ... (one bond) */
							if(sel[J]) {
								continue ;
							}
/* ... that has never been visited ... */
							D = d1[I][J] ;
							if(d1[k][J] != D+1) {
								continue ;
							}
/* ... so that J is on the same side chain as I, so that I is between k and J */
							if(D > Dm) {
/* a candidate Jm has been found */
								Jm = J ;
								Dm = D ;
/* remember that J as Jm and that D as Dm */
							}
						}
						chain(lev + 1, I, Jm, Dm, t) ;
/* print that side chain */
					} /* if, on side chain length being 1 or longer */
				} /* for, search for I attached to k and located on a side chain */
			} /* if k may bear a side chain */
			t = get_typ(k, kp) ;
/* bond type (order) along the main chain, between k and kp */
			w_bond(t) ;
/* write that bond */
			w_atom(kp) ;
#if DEBUG==1
			printf("***   %s\n", buf) ;
#endif
			k = kp ;
/* moving one bon away along the main chain */
		} /* for, on distance betwwen i and kp */
	} /* if */
	sel[j] = TRUE ;
/* last atom in main chain has been visited */
	if(lev) w_char(')') ; 
/* end of chain, this chain has a parent chain */
#if DEBUG==1
	printf("***   %s\n", buf) ;
#endif
}

/*
 * void w_char(char c) writes a single character in the SMILES string buffer
 */
void w_char(char c)
{
	sprintf(buf + ibuf, "%c", c) ; ibuf++ ;
}

/*
 * void w_bond(int t) writes a bond symbol in the SMILES string buffer
 */
void w_bond(int t)
{
	if(t != 1) {
 		sprintf(buf + ibuf, "%c", (t == 2) ? '=' : '#') ; ibuf++ ;
 	}
}

/*
 * void w_atom(int a) writes the element label of atom a in the SMILES string buffer
 */
void w_atom(int a)
{
	sprintf(buf + ibuf, "%s", label[a]) ;
	ibuf += strlen(label[a])  ;
}

/*
 * returns the bond order (or bond type) of the (a1,a2) bond
 */
int get_typ(int a1, int a2) 
{
	int i;

	for(i=0 ; i<cunb ; i++) {
		if(ue1[i] == a1 && ue2[i] == a2 || ue1[i] == a2 && ue2[i] == a1) {
			return ubtyp[i] ;
		}
	}
}


