/*
    file: LsdSrc/corr.c
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
#include "defs.h"

int	start2(void) ;
int	start3(void) ;
int	pathincludes(struct lcorr *pc1, struct lcorr *pc2) ;
int	pathequals(struct lcorr *pc1, struct lcorr *pc2) ;
int	pathinter(struct lcorr *pc1, struct lcorr *pc2) ;
void	reducevarcorr(int i, int x) ;
void	newcorr(int p1, int p2, int ph1, int ph2, int p3, int p4, int type) ;
void	newcorr2(int newp1, struct lcorr *pc1) ;
void	updatecorr(int corr) ;
int	insideln(int ln1, int ln2) ;
int 	inln(int x1, int *v2, int sz2) ;
void typesetcorr(struct lcorr *pc) ;

extern void	non_nul(int errnum, int v, ...) ; /* lecture.c */
extern int	lnpos(void) ; /* lecture.c */
extern int	estlie(int x, int y) ; /* start.c */
extern int  estlielie(int p1, int p2) ; /* start.c */
extern void	wrln(int i) ; /* start.c */
extern void	s_lier(int p1, int p2, int fl) ; /* start.c */

char hlpstr[100] ;

extern int	valcorr[MAX_CORR] ;
extern struct lcorr propcorr[MAX_CORR] ;
extern int	icorr ;
extern int	flverb ;
extern int	maxj ;
extern int	lnsize[MAX_LN] ;
extern int	ln[MAX_LN][MAX_LNBUF] ;
extern int	flcorrvar ;
extern int	lnbuf[MAX_LNBUF] ; /* lecture.c */
extern int	ilnbuf ; /* lecture.c */
extern int	flelim ; /* lsd.c */


int		start2(void)
/*
 * treatment of non variant HMBC and COSY correlations
 * returns TRUE only if start2() has changed something to the correlation set
 */
{
	int	i, k, inter, x ;
	int work = 0 ;
	struct lcorr *pc1, *pc2 ;
	int c11, c12, c21, c22 ;
	
/* first part : comparison of non variant with non variant correlations */
	for(i = 1 ; i < icorr ; i++) {
		if(valcorr[i] && (pc1 = &propcorr[i], (c11 = pc1->at1) > 0)) {
/* correlation i is not variant */
			c12 = pc1->at2 ;
			for(k = 0 ; k < i ; k++) {
				if(valcorr[k] && (pc2 = &propcorr[k], (c21 = pc2->at1) > 0)) {
/* comparison of correlation i with non variant correlation k, k < i */
					c22 = pc2->at2 ;
					if(c11 == c21 && c12 == c22) {
/* the correlation is duplicated */
						if(flverb > 1) {
							printf("There are two correlations between %d and %d\n", c11, c12) ;
						}
						if(!pathequals(pc1, pc2)) {
							inter = pathinter(pc1, pc2) ;
/* calcul of intersection of the pair of min/max coupling path lengths */
							if(inter < 0) {
								typesetcorr(pc2) ;
								non_nul(258, FALSE, hlpstr) ;
							}
/* if intersection is negative the pair of correlations is conflicting */
							updatecorr(i) ;
/* the same correlation with different min and max coupling path length */
							work++ ;
						}
						if(flverb > 1) {
							printf("One of the two correlations between %d and %d is removed\n", c11, c12) ;
						}
						valcorr[k] = FALSE ;
					}
				}
			}
		}
	}

/* second part : bond formation */
	for(i = 0 ; i < icorr ; i++) {
		if(valcorr[i] && (pc1 = &propcorr[i], (c11 = pc1->at1) > 0) && !pc1->j2ok && !pc1->jnok) {
/* correlation i is not variant and is equivalent to a bond */
			c12 = pc1->at2 ;
			if(flverb > 1) {
				printf("Bond formation between %d and %d because this correlation is j1 only\n", c11, c12) ;
			}
			s_lier(c11, c12, FALSE) ;
			if(flverb > 1) {
				printf("Correlation between %d and %d [%d-%d] is removed\n",c11, c12, pc1->min, pc1->max) ;
			}
			work++ ;
			valcorr[i] = FALSE ;
/* bond formation when only 2J is allowed (j1ok only) */
		}
	}

/* third and last part : comparison of non variant correlations with bonds */
	for(i = 0 ; i < icorr ; i++) {
		if(valcorr[i] && (pc1 = &propcorr[i], (c11 = pc1->at1) > 0)) {
			c12 = pc1->at2 ;
			if (estlie(c11, c12)) {
				if(!pc1->j1ok) {
					typesetcorr(pc1) ;
					non_nul(259, FALSE, hlpstr, c11, c12) ;
				}
				if(flverb > 1) {
					printf("Correlation between %d and %d [%d-%d] is removed because %d and %d are already bound\n",
						c11, c12, pc1->min, pc1->max, c11, c12) ;
				}
				valcorr[i] = FALSE ;
/* the correlation c11 - c12 does not bring anything because c11 and c12 are already bound */
/* c11 --------- c12 */
			} else if(x = estlielie(c11, c12)) {
				if(!pc1->j1ok && !pc1->j2ok) {
					typesetcorr(pc1) ;
					non_nul(290, FALSE, hlpstr, c11, c12) ;
				}
				if(pc1->j2ok) {
					if(flverb > 1) {
						printf("Correlation between %d and %d [%d-%d] is removed due to %d - %d - %d bonds\n",
							c11, c12, pc1->min, pc1->max, c11, x, c12) ;
					}
					valcorr[i] = FALSE ;
/* the correlation does not bring anything because c11 and c12 have a common neighbor */
/* c11 --------- x --------- c12 */
				}
			}
		}
	}
	return(work) ;
}


int		start3(void)
/*
 * treatment of variant HMBC and COSY correlations
 * returns TRUE only if start3() has changed something to the correlation set
 */
{
	int i, j, k, s, p1, lnp1j, x ;
	int work = 0 ;
	struct lcorr *pc1, *pc2 ;
	int c11, c12, c21, c22 ;
	int tmp, order, icorr0 ;

	if(!flcorrvar) return FALSE ;
/* if there is no variant correlation then return */
	
/* first part : comparison of variant with non variant correlations */
	icorr0 = icorr ;
	for(i = 0 ; i < icorr0 ; i++) {
		if(valcorr[i] && (pc1 = &propcorr[i], (p1 = -pc1->at1) > 0)) {
/* correlation i is variant */
			tmp = pc1->at2 ;
			s = lnsize[p1] ;
			for (j = 0 ; j < s && valcorr[i] ; j++) {
/* enumerate the members of the group */
				lnp1j = ln[p1][j] ;
				order = lnp1j < tmp ;
				c11 = order ? lnp1j : tmp ;
				c12 = order ? tmp : lnp1j ;
/* ensure that c11 < c12, because c21 < c22 (below) */
				for(k = 0 ; k < icorr0 && valcorr[i] ; k++) {
					if(valcorr[k] && (pc2 = &propcorr[k], (c21 = pc2->at1) > 0) && (c11 == c21) && (c12 == (c22 = pc2->at2))) {
/* correlation between the element j of the list p1 and p2 already exists among non variant correlations */
						if(pathequals(pc1, pc2) || pathincludes(pc1, pc2)) {
							if(flverb > 1) {
								printf("Correlation ( ") ;
								wrln(p1) ;
								printf(") %d [%d-%d] is removed because it is already explained by correlation %d %d [%d-%d]\n",
									tmp, pc1->min, pc1->max, c21, c22, pc2->min, pc2->max) ;
							}
							flcorrvar-- ;
/* one less variant correlation */
							valcorr[i] = FALSE ;
/* variant correlation becomes invalid because pairs of min/max are equal or
the non variant correlation is completely included in the variant one */
						} else if (pathinter(pc1, pc2) < 0) {
/* calculate intersection of the min/max coupling path lengths */
/* if intersection is negative a member of the list is conflicting */
							if(flverb > 1) {
								printf("Atom %d from correlation ( ", lnp1j) ;
								wrln(p1) ;
								printf(") %d [%d-%d] is conflicting with correlation %d %d [%d-%d]\n",
									c12, pc1->min, pc1->max, c21, c22, pc2->min, pc2->max) ;
							}
							reducevarcorr(i, lnp1j) ;
/* creates a new correlation with lnp1j removed from the ln[p1] list */
							work++ ;
							if(flverb > 1) {
								printf("Correlation ( ") ;
								wrln(p1) ;
								printf(") %d [%d-%d] is removed\n", c22, pc1->min, pc1->max) ;
							}
							valcorr[i] = FALSE ;
/* old variant correlation becomes invalid */
						} 
					}
				}
			}
		}
	}
	
/* second part : comparison of variant with variant correlations */
	for(i = 0 ; i < icorr ; i++) {
		if(valcorr[i] && (pc1 = &propcorr[i], (c11 = -pc1->at1) > 0)) {
/* testing valid variable correlation 1 indexed i ... */
		c12 = pc1->at2 ;
			for(k = 0 ; k < icorr && valcorr[i] ; k++) {
				if(
				valcorr[k] && (i != k) &&
				(pc2 = &propcorr[k], (c21 = -pc2->at1) > 0) &&
/* ... against valid variable correlation 2 indexed k ... */
				((c22 = pc2->at2) == c12) &&
/* ... with the same non variable atom as correlation 1. */
				((c11 == c21) || insideln(c21, c11)) &&
/* if the two lists are the same or if the list in correlation 2 is inside of the list of correlation 1 ... */
				(pathequals(pc1, pc2) || pathincludes(pc1, pc2)) ) {
/* ... and if the coupling path of correlation 1 covers the one of correlation 2 ... */
					flcorrvar-- ;
/* one less variant correlation */
					if(flverb > 1) {
						printf("Correlation ( ") ;
						wrln(c11) ;
						printf(") %d [%d-%d] is removed because it is already explicated by correlation ( ",
							c12, pc1->min, pc1->max) ;
						wrln(c21) ;
						printf(") %d [%d-%d]\n", c22, pc2->min, pc2->max) ;
					}
					valcorr[i] = FALSE ;
/* ... then correlation 1 is useless */
				}
			}
		}
	}
	
	
/* third part : comparison of variant correlations with bonds */	
	for(i = 0 ; i < icorr ; i++) {
		if(valcorr[i] && (pc1 = &propcorr[i], (c11 = -pc1->at1) > 0)) {
			s = lnsize[c11] ;
			c12 = pc1->at2 ;
			for (j = 0 ; j < s && valcorr[i] ; j++) {
				lnp1j = ln[c11][j] ;
				if(estlie(lnp1j, c12)) {
/* test if two correlating atoms are bound */
					if(!pc1->j1ok) {
/* if j1 is not allowed the conflicting member must be removed from the list */
						if(flverb > 1) {
							printf("Correlation ( ") ;
							wrln(c11) ;
							printf(") %d [%d-%d] is conflicting with bond %d - %d\n", c12, pc1->min, pc1->max, lnp1j, c12) ;
						}
						reducevarcorr(i, lnp1j) ;
/* creates a new correlation with lnp1j removed from the ln[p1] list */
						work++ ;
						if(flverb > 1) {
							printf("Correlation ( ") ;
							wrln(c11) ;
							printf(") %d [%d-%d] is removed\n", c12, pc1->min, pc1->max) ;
						}
						valcorr[i] = FALSE ;
/* variant correlation becomes invalid */
					} else {
						flcorrvar-- ;
/* one less variant correlation */
						if(flverb > 1) {
							printf("Correlation ( ") ;
							wrln(c11) ;
							printf(") %d [%d-%d] is removed because %d and %d are already bound\n",
								c12, pc1->min, pc1->max, lnp1j, c12) ;	
						}
						valcorr[i] = FALSE ;
/* variant correlation becomes invalid because j1 is allowed */
					}
				} else if(x = estlielie(lnp1j, pc1->at2)) {
/* test if two correlating atoms are separated by a common neighbor */
					if(pc1->j2ok) {
						flcorrvar-- ;
/* one less variant correlation */
						if(flverb > 1) {
							printf("Correlation ( ") ;
							wrln(c11) ;
							printf(") %d [%d-%d] is removed due to %d - %d - %d bonds\n",
								c12, pc1->min, pc1->max, lnp1j, x, c12) ;	
						}
						valcorr[i] = FALSE ;
/* if j2 is allowed the correlation becomes invalid */
					} else if(!pc1->j1ok && !pc1->j2ok) {
/* if jn only is allowed the conflicting member must be removed from the list */
						if(flverb > 1) {
							printf("Correlation ( ") ;
							wrln(c11) ;
							printf(") %d [%d-%d] is conflicting with bonds %d - %d - %d\n",
								c12, pc1->min, pc1->max, lnp1j, x, c12) ;
						}
						reducevarcorr(i, lnp1j) ;
/* creates a new correlation with lnp1j removed from the ln[p1] list */
						work++ ;
						if(flverb > 1) {
							printf("Correlation ( ") ;
							wrln(c11) ;
							printf(") %d [%d-%d] is removed\n", c12, pc1->min, pc1->max) ;
						}
						valcorr[i] = FALSE ;
/* old variant correlation becomes invalid */
					}
				}	
			}
		}
	}
	
	return(work) ;
}


int	pathincludes(struct lcorr *pc1, struct lcorr *pc2)
/*
 * returns true if the coupling path of correlation 1 completely covers the
 * coupling path of correlation 2
 */
{
	int pc1min, pc1max, pc2min, pc2max ;
	
	pc1min = pc1->min ;
	pc1max = pc1->max ;
	pc2min = pc2->min ;
	pc2max = pc2->max ;
	return (
		((!pc1max) && (pc2min >= pc1min)) ||
		((pc1max && pc2max) && (pc2min >= pc1min) && (pc2max <= pc1max))
	) ;
}


int	pathequals(struct lcorr *pc1, struct lcorr *pc2)
/*
 * returns true if the coupling path of correlation 1 is identical to
 * coupling path of correlation 2
 */
{
	return (pc1->min == pc2->min) && (pc1->max == pc2->max) ;
}


int	pathinter(struct lcorr *pc1, struct lcorr *pc2)
/*
 * calculates the intersection of a pair of min/max coupling path lengths.
 * returns the max-min of the intersection.
 * if the correlation pointed by pc1 is not variant, then pc1 is modified
 * so that min and max are those of the intersection, if not empty (inter < 0) 
 */
{
	int	p3_1, p3_2, p4_1, p4_2, maxp3, minp4, inter ;

	p3_1 = pc1->min ;
	p4_1 = pc1->max ;
	p3_2 = pc2->min ;	
	p4_2 = pc2->max ;
	
	maxp3 = (p3_1 >= p3_2) ? p3_1 : p3_2 ;
/* Maximum of minimum coupling path lengths */
	if(!p4_1 || !p4_2) {
		minp4 = (p4_1 >= p4_2) ? p4_1 : p4_2 ;
	} else {
		minp4 = (p4_1 <= p4_2) ? p4_1 : p4_2 ;
	}
/* Minimum of maximum coupling path lengths */
	if(!p4_1 && !p4_2) {
		inter = 1 ;
	} else {
		inter = minp4 - maxp3 ;
	}
/* intersection is the difference of Min of max and Max of min */
	if((inter >= 0) && (pc1->at1 > 0)) {
		if(flverb > 1) {
			printf("Min and max coupling path lengths of correlation between %d and %d are redefined "
"as the intersection of [%d-%d] and [%d-%d]\n", pc1->at1, pc1->at2, p3_1, p4_1, p3_2, p4_2) ;
			printf("Correlation is redefined : %d %d [%d-%d]\n", pc1->at1, pc1->at2, maxp3, minp4) ;
		}
		pc1->min = maxp3 ;
		pc1->max = minp4 ;
/* puts new values in the correlation table only if it concerns a non variant correlation */
		pc1->origin = 1 ;
/* Min and Max coupling path lengths may have been modified */
	}
	return inter ;
}


void	reducevarcorr(int i, int x)
/*
 * considering variant correlation i, atom x is excluded from the list of alternatives
 * and a new correlation is created.
 * This is a variant correlation if more than one possibility remains
 * and a fixed correlation otherwise.
 */
{
	struct lcorr *pc1 ;
	int p1, s, l, newp1 ;
	
	pc1 = &propcorr[i] ;
	p1 = -pc1->at1 ;
	if(flverb > 1) {
		printf("Atom %d is removed from the list ( ", x) ;
		wrln(p1) ;
		printf(") correlating with %d [%d-%d]\n", pc1->at2, pc1->min, pc1->max) ;
	}
	s = lnsize[p1] ;
	for(l = 0, ilnbuf = 0 ; l < s ; l++) {
		if(ln[p1][l] != x) {
			lnbuf[ilnbuf] = ln[p1][l] ;
			ilnbuf++ ;
		}
	}
/* the conflicting member is removed from the list */
	if(ilnbuf == 1) {
		newp1 = lnbuf[0] ;
		if(flverb > 1) {
			printf("A new non variant correlation is defined : %d %d [%d-%d]\n", newp1, pc1->at2, pc1->min, pc1->max) ;
		}
/* the new correlation is non variant*/
		flcorrvar-- ;
/* one less variant correlation */
	} else {
		newp1 = -lnpos() ;
		if(flverb > 1) {
			printf("A new variant correlation is defined : ( ") ;
			wrln(-newp1) ;
			printf(") %d [%d-%d]\n", pc1->at2, pc1->min, pc1->max) ;
		}	
/* a new list is created if the list has more than one element */
	}
	
	newcorr2(newp1, pc1) ;
	updatecorr(icorr) ;
/* a new correlation is created */
	valcorr[icorr++] = TRUE ;
/* validity of the new correlation */
}


void	newcorr(int p1, int p2, int ph1, int ph2, int p3, int p4, int type)
/*
 * validates a new correlation and its min/max coupling path lengths
 */
{
	struct lcorr *pc ;
	int inp3, inp4 ;

	inp3 = p3 ;
	inp4 = p4 ;
	non_nul(291, icorr - MAX_CORR, MAX_CORR) ;
	non_nul(253, p3 != 1) ;
	non_nul(254, p4 != 1) ;
	if(type) {
/* type = 1 for COSY and type = 0 for HMBC */
		non_nul(242, p3 != 2) ;
		non_nul(243, p4 != 2) ;
	}
	non_nul(255, p3) ;
	if(p4 > 0) non_nul(256, p3 <= p4, p3, p4) ;
	if(p3 > 0) p3 = type ? p3 - 2 : p3 - 1 ;
	if(p4 > 0) p4 = type ? p4 - 2 : p4 - 1 ;
/* min and max coupling path lengths of COSY and HMBC correlations 
 * are stored as lengths between heavy atoms */
	if(p4 == -1) p4 = p3 ;
	if(p3 == -1 && p4 == -1) {
		p3 = 1 ;
		p4 = type ? 1 : (maxj ? maxj - 1 : 0) ;
	}
	if(!p4 && maxj) p4 = maxj - 1 ;
	if(maxj) non_nul(257, p4 <= maxj - 1, maxj + 1, maxj) ;
	pc = &propcorr[icorr] ; 
	pc->at1 = p1 < p2 ? p1 : p2 ;
	pc->at2 = p1 > p2 ? p1 : p2 ;
	pc->inputat1 = p1 ;
	pc->inputat2 = p2 ;
	pc->inputath1 = ph1 ;
	pc->inputath2 = ph2 ;
	pc->inputmin = inp3 ;
	pc->inputmax = inp4 ;
	pc->min = p3 ;
	pc->max = p4 ;
	pc->origin = 0 ;
	pc->type = type ;
}


void	newcorr2(int newp1, struct lcorr *pc1)
/*
 * validates a new variant correlation when an atom is removed from a list and a new list is created
 */
{
	struct lcorr *pc2 ;
	
	non_nul(291, icorr - MAX_CORR, MAX_CORR) ;
	pc2 = &propcorr[icorr] ; 
	pc2->at1 = newp1 ;
	pc2->at2 = pc1->at2 ;
	pc2->inputat1 = pc1->inputat1 ;
	pc2->inputat2 = pc1->inputat2 ;
	pc2->inputath1 = pc1->inputath1 ;
	pc2->inputath2 = pc1->inputath2 ;
	pc2->inputmin = pc1->inputmin ;
	pc2->inputmax = pc1->inputmax ;
	pc2->min = pc1->min ;
	pc2->max = pc1->max ;
	pc2->origin = 2 ;
	pc2->type = pc1->type ;
}


void	updatecorr(int corr)
/*
 * checks if a correlation is j1ok, j2ok or jnok from its min and max coupling path lengths
 */
{
	struct lcorr *pc ;
	
	pc = &propcorr[corr] ; 
	pc->j1ok = FALSE ;
	pc->j2ok = FALSE ;
	pc->jnok = FALSE ;
	if(pc->min == 1) {
		pc->j1ok = TRUE ;
		if(pc->max != 1) {
			pc->j2ok = TRUE ;
			if(pc->max != 2) pc->jnok = TRUE ;
		}
	} else if(pc->min == 2){
		pc->j2ok = TRUE ;
		if(pc->max != 2) pc->jnok = TRUE ; 
	} else pc->jnok = TRUE ;
	if(!pc->j1ok && !pc->j2ok && !flelim) {
		typesetcorr(pc) ;
		non_nul(293, FALSE, hlpstr) ;
	}
}


int 	insideln(int ln1, int ln2)
/*
 * returns true only if ln[ln1] is strictly inside of ln[ln2].
 * stricly inside means inside and not equal.
 * lists ln[ln1] and ln[ln2] must be sorted in ascending order 
 * and must have no duplicated elements.
 */
{
	int *v1, *v2 ;
	int sz1, sz2 ;
	int i1, x1 ;
	int found ;
	
	v1 = ln[ln1] ;
	v2 = ln[ln2] ;
/* is v1 strictly included in v2? */
	sz1 = lnsize[ln1] ;
	sz2 = lnsize[ln2] ;
	if (sz1 >= sz2) {
/* if v1 is strictly included in v2 then v1 must have strictly less elements than v2 */
		return FALSE ;
	}
	for (found = TRUE, i1 = 0 ; (i1 < sz1) && (x1 = v1[i1], found = inln(x1, v2, sz2)) ; i1++)
		;
/* x1 from v1 will be looked for in v2. All x1 have to be found in v2 */
	return found ;
}


int	inln(int x1, int *v2, int sz2)
/*
 * returns TRUE if x1 belongs to the strictly sorted numeric list v2 of size sz2
 */
{
	int i2, x2 ;
	int found ;
	int v2min, v2max ;
	
	v2min = v2[0] ;
	v2max = v2[sz2 -1] ;
	if ((x1 < v2min) || (x1 > v2max)) {
/* x1 is outside of the range of values in v2 */
		found = FALSE ;
	} else {
		for (i2 = 0 ; (i2 < sz2) && ((x2 = v2[i2]) < x1) ; i2++)
			;
/* the loop stops when x2 from v2 is equal or greater to x1 */
		found = (x1 == x2) ;
	}
	return found ;
}


void typesetcorr(struct lcorr *pc) {
/*
 * typeset
 */
	char *w ;
	int p1, i, s ;
	int nmin, nmax ;
	
	p1 = pc->inputat1 ;
	w = hlpstr ;
	if(pc->type) {
		nmin = 2 ;
		nmax = pc->max ? 2 : 0 ;
		w += sprintf(w, "H-") ;
		p1 = pc->inputath1 ;
	} else {
		nmin = 1 ;
		nmax = pc->max ? 1 : 0 ;
	}
	if(p1 > 0) {
/* it is a non variant correlation */
		w += sprintf(w, "%d H-%d [%d-%d]", p1, pc->inputath2, pc->min+nmin, pc->max+nmax) ;
	} else {
/* it is a variant correlation */
		p1 = -p1 ;
		s = lnsize[p1] ;
		w += sprintf(w, "(") ;
		for(i = 0 ; i < s ; i++) {
			w += sprintf(w, "%d%c", ln[p1][i], i == (s - 1) ? ')' : ' ') ;
		}
		w += sprintf(w, " H-%d [%d-%d]", pc->inputath2, pc->min+nmin, pc->max+nmax) ;
	}
}
