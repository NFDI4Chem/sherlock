/*
    file: LsdSrc/start.c
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

void	start(void) ;
void	start1(int ic, int o, int p1, int p2, int p3, int p4, int p5) ;
void	start4(void) ;
void	atomcorr(int corr) ;
void	atomcorrvar(int corr) ;
void	corrvar(int ic, int o) ;
void	s_lier(int p1, int p2, int fl) ;
int	s_carb(int x) ;
int	s_hete(int x) ;
int	s_sp3(int x) ;
int	s_sp2(int x) ;
int	s_sp(int x) ;
int s_char(int x) ;
int s_charpos(int x) ;
int s_charneg(int x) ;
int	s_quat(int x) ;
int	s_ch(int x) ;
int	s_ch2(int x) ;
int	s_ch3(int x) ;
int	s_full(int x) ;
int s_elem(int x, int y) ;
int	s_ge(int x, int y) ;
int	s_le(int x, int y) ;
int	s_gt(int x, int y) ;
int	s_lt(int x, int y) ;
int	s_ou(int x, int y) ;
int	s_et(int x, int y) ;
int	s_diff(int x, int y) ;
void	s_prop(int p, int (*f)(int)) ;
void	s_cmp(int p1, int p2, int (*f)(int, int)) ;
void	s_ens(int p1, int p2, int p3, int (*f)(int, int)) ;
int	*s_adr(int p, int *al) ;
void	s_useprop(int x, int n, int l, int c) ;
void	lier_prop(int p1, int p2) ;
int	estvalide(int x) ;
int	hestvalide(int x) ;
int	estlie(int x, int y) ;
int estlielie(int p1, int p2) ;
void	about(void) ;
void	wrln(int i) ;
void	wrln2(int y) ;
int	stoss(int i) ;

extern void	non_nul(int errnum, int v, ...) ; /* lecture.c */
extern void	check_list(int l)  ;
extern int	chklsd(void) ;
extern void	sub_about(void) ; /* sub.c */
extern void	test_open_path(int p) ; /* sub.c */
extern void	test_open_deff(int f, int p) ; /* sub.c */
extern void	test_open_skel(int f, int p) ; /* sub.c */
extern void	test_open_fexp(int p) ; /* sub.c */
extern void sub_start1(int ic) ; /* sub.c */
extern int	setcnx(int x, int y, int fl) ; /* connexe.c */
extern int	lnpos(void) ; /* lecture.c */
extern int	start2(void) ; /* corr.c */
extern int	start3(void) ; /* corr.c */
extern void	newcorr(int p1, int p2, int ph1, int ph2, int p3, int p4, int type) ; /* corr.c */
extern void	updatecorr(int corr) ; /* corr.c */

extern struct atom at[MAX_ATOMES] ;
extern int	list[MAX_LIST][MAX_ATOMES] ;
extern struct comtab coms[MAX_COM] ;
extern struct noy nucl[MAX_NUCL] ;
extern int	lignes_lues ;
extern struct crt cartes[MAX_CARTES] ;
extern int	valcorr[MAX_CORR] ;
extern struct lcorr propcorr[MAX_CORR] ;
extern int	icarte ;
extern int	max_atomes ;
extern int	icorr ;
extern int	flabout ;
extern int	flhistoire ;
extern int	fldisp ;
extern int	flverb ;
extern int	flpart ;
extern int	flstep ;
extern int	flwork ;
extern int	fldupl ;
extern int	flsubsgood ;
extern int	flsubsbad ;
extern int  flbredt ;
extern int	elim ;
extern int	flelim ;
extern int	maxj ;
extern int	flfilt ;
extern int	base[MAX_ATOMES] ;
extern int	lnsize[MAX_LN] ;
extern int	ln[MAX_LN][MAX_LNBUF] ;
extern int	flcorrvar ;
extern int	nivdisp ;
extern struct ss_atom sub[MAX_ATOMES] ;
extern int	nssl ;
extern int	nass ;
extern float	reels[MAX_REELS] ;
extern int	numeqdbbonds ; /* doubles.c */
extern int numsp2 ; /* doubles.c */
extern int	numrings ; /* chklsd.c */
extern int	fldeff ; /* lsd.c */
extern int	flfexp ; /* lsd.c */
extern int flmaxfrag ; /* lsd.c */
extern dyntext	*strtab ; /* lsd.c */
extern int	flcnx ; /* lsd.c */
extern int	ncnx ; /* connexe.c */
extern int	natomes ; /* connexe.c */
extern int	flmaxstruct ; /* lsd.c */
extern int	flmaxtime ; /* lsd.c */
extern int	flCclasses ; /* lsd.c */
extern int	lnbuf[MAX_LNBUF] ; /* lecture.c */
extern int	ilnbuf ; /* lecture.c */
extern dyntext	*pathtab ; /* lsd.c */
extern char *countfile ; /* lsd.c */
extern char *stopfile ; /* lsd.c */
extern float reels[MAX_REELS] ; /* lsd.c */


int	htoc[MAX_ATOMES] ; /* H index to C (or X) index from HMQC data */
int	ctoh[MAX_ATOMES][2] ; /* C (or X) index to H index from HMQC data */
int	deflist[MAX_LIST] ; /* true if the logical list is defined */

void	start(void)
/*
 * initialisations and analysis of the commands read by lecture()
 */
{
	int	i, ic, ip, w2, w3 ;
	int 	nval ;

	for (i = 0 ; i < MAX_ATOMES ; i++) {
		base[i] = htoc[i] = ctoh[i][0] = ctoh[i][1] = 0 ;
	}
/* resets the base (see lsd.c) and the C (or X) <-> H conversion tables */
	for (i = 0 ; i < MAX_ATOMES ; i++) {
		at[i].utile = FALSE ;
/* not initialized by a MULT command */
		at[i].nc = at[i].nlpr = at[i].nlia = 0 ;
/* no correlation, no property list, no neighbor */
		at[i].equstat = -1 ;
/* do not belong to an equivalence class */
		at[i].corrgrp = 0 ;
/* does not correlate with a group of unresolved X resonances */
		at[i].ingrp = FALSE ;
/* does not belong to a group of unresolved X resonances */
		at[i].icnx = 0 ;
/* no connected unit index assigned */
		at[i].ics = at[i].ihcs[0] = at[i].ihcs[1] = -1 ;
/* no chemical shifts */
	}
/* resets atoms characteristics */
	ncnx = 1 ;
/* ncnx is always the first available connected unit index */ 
	natomes = 0 ;
/* no atom yet */

	for (i = 0 ; i < MAX_CORR ; i++) {
		valcorr[i] = FALSE ;
	}
/* no valid correlation */
	for (i = 0 ; i < MAX_LIST ; i++) {
		deflist[i] = FALSE ;
	}
/* no defined logical list, but list[][] may already contains TRUE values */
	for (ip = 0 ; ip < MAX_PRIO ; ip++) {
/* analyzing commands following ascending priority order */
		for (ic = 0 ; ic < icarte ; ic++) {
			if (coms[cartes[ic].opnum].prior == ip) {
/* printf("# %d priority %d command %s\n", ic+1, ip, coms[cartes[ic].opnum].op) ; */
				lignes_lues = cartes[ic].ln ;
				cartes[ic].done = TRUE ;
				start1(	ic, 
					cartes[ic].opnum, 
					cartes[ic].pars[0], 
					cartes[ic].pars[1], 
					cartes[ic].pars[2],
					cartes[ic].pars[3],
					cartes[ic].pars[4]) ;
/* start1 does the real job */
			}
		}
	}
	non_nul(199, natomes) ;
	countfile = (countfile == NULL) ? SOLUTION_COUNT_FILE : countfile ;
	stopfile = (stopfile == NULL) ? LSD_STOP_FILE : stopfile ;
	if (flcorrvar) {
		for (ic = 0 ; ic < icarte ; ic++) {
			if (!cartes[ic].done) {
				corrvar(ic, cartes[ic].opnum) ;
			}
		}
/* if there are variants in HMBC and COSY correlation, the input analysis is completed by corrvar */
	}
	
	if(flverb > 1) {
		printf("Analysis of HMBC and COSY data.\n") ;
	}
	
	do {
		w2 = start2() ;
/* analysis of non variant correlations */
		w3 = start3() ;
/* analysis of variant correlations */
	} while (w2 || w3) ;
	
	if(flverb > 1) {
		printf("Analysis of HMBC and COSY data done.\n") ;
	}
	
	start4() ;
/* storage of correlations as atom properties */
	
 	non_nul(390, 1 - numsp2 % 2, numsp2) ;
/* the number of sp2 atoms cannot be odd */
	nval = chklsd() ;
	non_nul(391, 1 - nval % 2, nval) ;
/* the sum of all valencies cannot be odd, its half is the number of bonds */
	non_nul(392, (!flcnx) || (numrings >= 0), numrings) ;
/* only if a connected solution is required then numrings cannot be negative */ 
	if (flverb) {
		printf("start ok\n") ;
	/*	about() ; */
	}
}


void	start1(int ic, int o, int p1, int p2, int p3, int p4, int p5)
{
	struct atom *a ;
	int	i, j, cardinal, x ;
	int ph1, ph2 ;

	switch (o) {

	case 0 : /* VALE : valency A I R */
		non_nul(204, fldupl != 2) ;
		non_nul(200, !nucl[p1].val[1], nucl[p1].sym, nucl[p1].val[1]) ;
/* only concerns atoms for which valence is not defined */
		non_nul(201, p2, nucl[p1].sym) ;
/* valence cannot be 0 */
		non_nul(202, p2 <= MAX_NLIA, nucl[p1].sym, MAX_NLIA, p2) ;
/* valence cannot be higher than MAX_NLIA */
		nucl[p1].val[1] = p2 ;
/* redefines the valency */
		if(p3 >= 0) {
			nucl[p1].mass = reels[p3] ;
			non_nul(203, nucl[p1].mass >= 1.0, nucl[p1].sym, nucl[p1].mass) ;
		}
/* defines the mass of the atom, unless p3 < 0 */
		break ;

	case 1 : /* MULT : multiplicity I A I I Z */
		non_nul(210, p1 > 0) ;
		non_nul(211, p1 < MAX_ATOMES, MAX_ATOMES-1, p1) ;
		non_nul(212, (p3 >= 1) && (p3 <= 3), p1, p3) ;
		if((p2 == 3) || (p2 == 13)) {
			non_nul(215, p3 != 3) ;
/* atom N5(index 3 of the nuclei table) and P5(13) must not be sp3 because maximum number of neighbors is 4 */
		}
		if(p2 == 6) {
			non_nul(216, p3 == 1) ;
/* atom S6(6) must be sp */
		}
		a = &at[p1] ;
		a->utile = TRUE ;
		a->element = p2 ;
		a->charge = p5 ;
		a->valence = nucl[p2].val[p5+1] ;
/* valence is function of atomic symbol and charge */
		non_nul(213, a->valence > 0, p1, p5) ;
		a->hybrid = 3 - p3 ; /* 0 for sp3, 1 for sp2 and 2 for sp */
		a->mult = p4 ;
		a->nlmax = a->valence - a->hybrid - a->mult ; /* neighbours */
		non_nul(214, a->nlmax > 0, p1, a->valence, p3, p4) ;
		if (p1 > max_atomes) {
			max_atomes = p1 ;
/* atoms have not be to defined in the natural order. never tested */
		}
		numeqdbbonds += a->hybrid ;
		if(a->hybrid == 1) {
			numsp2++ ;
/* one more sp2 atom */
		}
		nucl[p2].count++ ;
		nucl[MAX_NUCL-1].count += p4 ;
/* hydrogen atom count is kept in the .count field of A-type atoms */
		a->icnx = ncnx++ ;
/* each atom initially receives a different connected unit index */
		natomes++ ;
		break ;

	case 2 : /* LIST : list L */
		cardinal = 0 ;
		for (x = 1 ; x <= max_atomes ; x++) {
			if (list[p1][x]) {
				non_nul(260, at[x].utile, x, p1) ;
				cardinal++;
			}
		}
		for (x = max_atomes + 1 ; x < MAX_ATOMES ; x++)
			non_nul(260, !list[p1][x], x, p1) ;
		non_nul(261, cardinal, p1) ;
		/*	printf("** liste %d **", p1) ; check_list(p1) ; printf("\n") ;*/
		deflist[p1] = TRUE ;
/* the list p1 is defined once all its elements are defined */
		break ;

	case 3 : /* CARB : carbons L */
		s_prop(p1, s_carb) ;
		break ;

	case 4 : /* HETE : heteroatoms L */
		s_prop(p1, s_hete) ;
		break ;

	case 5 : /* SP3  L*/
		s_prop(p1, s_sp3) ;
		break ;

	case 6 : /* SP2 L */
		s_prop(p1, s_sp2) ;
		break ;

	case 7 : /* QUAT : quaternary L */
		s_prop(p1, s_quat) ;
		break ;

	case 8 : /* CH   : L*/
		s_prop(p1, s_ch) ;
		break ;

	case 9 : /* CH2  : L */
		s_prop(p1, s_ch2) ;
		break ;

	case 10 : /* CH3  : L */
		s_prop(p1, s_ch3) ;
		break ;

	case 11 : /* FULL : all the defined atoms L */
		s_prop(p1, s_full) ; 
		break ;

	case 12 : /* GREQ : greater or equal L I */
		s_cmp(p1, p2, s_ge) ; 
		break ;

	case 13 : /* LEEQ : less or equal L I */
		s_cmp(p1, p2, s_le) ; 
		break ;

	case 14 : /* GRTH : strictly greater than L I */
		s_cmp(p1, p2, s_gt) ; 
		break ;

	case 15 : /* LETH : strictly small than L I */
		s_cmp(p1, p2, s_lt) ; 
		break ;

	case 16 : /* UNIO : union B B L */
		s_ens(p1, p2, p3, s_ou) ; 
		break ;

	case 17 : /* INTE : intersection B B L */
		s_ens(p1, p2, p3, s_et) ; 
		break ;

	case 18 : /* DIFF : difference B B L */
		s_ens(p1, p2, p3, s_diff) ; 
		break ;

	case 19 : /* PROP : properties B I L */
		if (p1 > 0) {
/* atom indexes start at 1 */
			non_nul(270, estvalide(p1), p1) ;
			s_useprop(p1, p2, p3, p4) ;
		} else {
/* p1 is -(logical list index), ok for p1 = 0 included */
			p1 = -p1 ;
/* true logical list index, 0 included */
			non_nul(271, deflist[p1], p1) ;
			for (x = 1 ; x <= max_atomes ; x++) {
				if (list[p1][x]) {
					s_useprop(x, p2, p3, p4) ;
/* applies properties to all the list's members */
				}
			}
		}
		break ;

	case 20 : /* BOND : bond I I */
		s_lier(p1, p2, 1) ;
		break ;
		
	case 21 : /* HSQC : hmqc I I */
	case 44 : /* HMQC : hsqc I I */
		non_nul(230, estvalide(p1), p1, p2, p1) ;
		non_nul(231, p2, p1, p2) ;
		non_nul(232, p2 < MAX_ATOMES, p1, p2, MAX_ATOMES-1) ;
   		for(i=0 ; i < 2 && ctoh[p1][i] ; i++)
   			;
   		non_nul(233, i != 2, p1, p2, p1, ctoh[p1][0], ctoh[p1][1]) ;
		non_nul(234, !htoc[p2], p1, p2, htoc[p2], p2) ;
		non_nul(235, i != 1 || at[p1].mult == 2, p1, p2, p1, ctoh[p1][0], nucl[at[p1].element].sym) ;
		ctoh[p1][i] = p2 ;
		htoc[p2] = p1 ;
		break ;

	case 22 : /* COSY : cosy V I O O */
		non_nul(240, hestvalide(p2), p1, p2, p2) ;
		if (p1 < 0) {
			flcorrvar++ ;
			cartes[ic].done = FALSE ;
			break ;
/* storage and treatment of fuzzy COSY correlations is delayed */
		} else {
			non_nul(240, hestvalide(p1), p1, p2, p1) ;		
			ph1 = p1 ;
			ph2 = p2 ;
			p1 = htoc[p1] ;
			p2 = htoc[p2] ;
			if(p1 == p2) {
				if(at[p1].mult != 2) {
					non_nul(241, FALSE, ph1, ph2, p1) ;
/* p1 must not correlate with itself */
				} else {
					break ;
/* correlation is valid but not useful, so it is discarded silently */
				}
			}
			
			newcorr(p1, p2, ph1, ph2, p3, p4, 1) ;
			updatecorr(icorr) ;
/* validates a new correlation and its min/max coupling path lengths */
			valcorr[icorr++] = TRUE ;

/* global validity of the new correlation */
/* warning: the validity of icorr is not checked */
		}
		break ;
	case 23 : /* HMBC : hmbc V I O O */
		non_nul(250, hestvalide(p2), p1, p2, p2) ;
		if (p1 < 0) {
			flcorrvar++ ;
			cartes[ic].done = FALSE ;
			break ;
/* storage and treatment of fuzzy HMBC correlations is delayed */
		} else {
			non_nul(251, estvalide(p1), p1, p2, p1) ;			
			ph1 = 0 ;
			ph2 = p2 ;
/* stores the hydrogen number of input file */ 
			p2 = htoc[p2] ;
/* transforms the proton reference into a carbon reference */
			non_nul(252, p2 != p1, p1, ph2, p1) ;
/* p1 must not correlate with p1 */
		
			newcorr(p1, p2, ph1, ph2, p3, p4, 0) ;
			updatecorr(icorr) ;
/* validates a new correlation and its min/max coupling path lengths */
			valcorr[icorr++] = TRUE ;

/* global validity of the new correlation */
/* warning: the validity of icorr is not checked */
		}
		break ;
	case 24 : /* ENTR : entry I */
		flabout = p1 ;
/* in order to display what has been understood from the input file */
		break ;
	case 25 : /* HIST : history I */
		flhistoire = p1 ;
/* in order to display a trace from the resolution mechanism */
		break ;
	case 26 : /* DISP : display I */
		fldisp = p1 ;
/* result as bonds or for JMNMOL */
		break ;
	case 27 : /* VERB : verbose I */
		flverb = p1 ;
		break ;
	case 28 : /* PART : partial I */
		flpart = p1 ;
/* in order to get only partial solutions */
		break ;
	case 29 : /* STEP : step I */
		flstep = p1 ;
/* for the single-step mode */
		break ;
	case 30 : /* WORK : work I */
		flwork = p1 ;
/* in order to do the job */
		break ;
	case 31 : /* MLEV maximum level I*/
		nivdisp = p1 ;
/* stop at a given heap level */
		break ;
	case 32 : /* SSTR : substructure S A V V */
		sub_start1(ic) ;
/* defines a subatom, see sub.c */
		break ;
	case 33 : /* ASGN : assign S I */
		sub_start1(ic) ;
/* updates substructure assignment information, see sub.c */
		break ;
	case 34 : /* LINK : link S S */
		sub_start1(ic) ;
/* set bonds between subatoms, see sub.c.
Validations could be improved */
		break ;
	case 35 : /* DUPL : duplicate solutions I */
		fldupl = p1 ;
		non_nul(375, (fldupl >= 0) && (fldupl <= 2), fldupl) ;
		break ;
	case 36 : /* SUBS : subs I */
		flsubsgood = (p1 > 0) ;
		flsubsbad = (p1 < 0) ;
		break ;
	case 37 : /* ELIM : eliminate correlations I I */
		elim = p1 ;
		flelim = (elim > 0) ;
		maxj = p2 ;
		if(maxj) {
			non_nul(380, maxj > 3, maxj) ;
		}
		break ;
	case 38 : /* FILT : filter substructure I */
		flfilt = p1 ;
		break ;
	case 39 : /* DEFF : define external substructure F C */
		non_nul(400, !dyntext_index_exists(strtab, p1), p1) ;
		test_open_deff(p1, p2) ; /* error 401 */
		dyntext_setindex(strtab, p2, p1) ;
		if(p1 > flmaxfrag) {
			flmaxfrag = p1 ;
		}
		fldeff++ ;
		break ;
	case 40 : /* FEXP : filter (logical) expression C */
		non_nul(402, !dyntext_index_exists(strtab, -1)) ;
		test_open_fexp(p1) ; /* errors 403 to 406 */
		break ;
	case 41 : /* CNTD : connected solutions I */
		flcnx = p1 ;
		break ;
	case 42 : /* MAXS : max number of solutions, 0 is no limit */
		flmaxstruct = p1 ;
		break ;
	case 43 : /* MAXT : max run time in seconds, 0 is no limit */
		flmaxtime = p1 ;
		break ;
/* next case is 45 because 44 is for HSQC and is grouped with 21: HMQC */
	case 45 : /* CCLA : equivalence classes for carbon atoms */
		flCclasses = p1 ;
		break ;
	case 46 : /* PATH :  directory for the finding of skeletons (see SKEL, below) C */
		test_open_path(p1) ; /* error 430 */
		break ;
	case 47 : /* SKEL : define external skeleton F C */
		non_nul(400, !dyntext_index_exists(strtab, p1), p1) ;
		dyntext_setindex(strtab, p2, p1) ;
		test_open_skel(p1, p2) ; /* errors 431 to 436 */
		if(p1 > flmaxfrag) {
			flmaxfrag = p1 ;
		}
		fldeff++ ;
		break ;
	case 48 : /* COUF : solution count file  P */
		non_nul(445, countfile == NULL, countfile) ;
		countfile = dyntext_getstr(pathtab, p1) ;
		break ;
	case 49 : /* STOF : lsd stop file  P */
		non_nul(446, stopfile == NULL, stopfile) ;
		stopfile = dyntext_getstr(pathtab, p1) ;
		break ;
	case 50 : /* SP   : list of sp atoms L */
		s_prop(p1, s_sp) ;
		break ;
	case 51 : /* CHAR : list of charged atoms L */
		s_prop(p1, s_char) ;
		break ;
	case 52 : /* CPOS : list of positively charged atoms L */
		s_prop(p1, s_charpos) ;
		break ;
	case 53 : /* CNEG : list of negatively charged atoms L */
		s_prop(p1, s_charneg) ;
		break ;
	case 54 : /* ELEM : list of element an type L A */
		s_cmp(p1, p2, s_elem) ;
		break ;
	case 55 : /* SHIX : chemical shift of non-hydrogen atoms I R */
		non_nul(300, estvalide(p1), p1) ;
		non_nul(301, at[p1].ics == -1, p1) ;
		at[p1].ics = p2 ;
		break ;
	case 56 : /* SHIH : chemical shift of hydrogen atoms I R */
		non_nul(302, hestvalide(p1), p1) ;
		for(i = 0 ; i < 2 && at[htoc[p1]].ihcs[i] != -1 ; i++)
			;
		non_nul(303, i != 2, p1) ;
		non_nul(304, i != 1 || at[htoc[p1]].mult == 2, p1, nucl[at[p1].element].symout) ;
		at[htoc[p1]].ihcs[i] = p2 ;
		break ;
	case 57 : /* BRUL : bredt rule I */
		flbredt = p1 ;
/* test for anti-Bredt structure */
		break ;
	}
} /* start1 */


void	corrvar(int ic, int o)
/*
 * pre-storage of variant HMBC and COSY correlations
 */
{
	int	p1, p2, ph1, ph2, p3, p4, type ;
	int s, i ;
	int *lnp1 ;
	int ok = TRUE ;
		
	type = 0 ;
	ph1 = 0 ;
	p1 = cartes[ic].pars[0] ;
/* carbon list reference */	
	ph2 = cartes[ic].pars[1] ;
	p2 = htoc[cartes[ic].pars[1]] ;
/* carbon from proton through HMQC data */
	p3 = cartes[ic].pars[2] ;
	p4 = cartes[ic].pars[3] ;
/* min and max coupling path lengths */
	
	lnp1 = ln[-p1] ;
	s = lnsize[-p1] ;
	
	if(o == 22) {
/* it is a COSY correlation */
		type = 1 ;
		ph1 = p1 ;
		for(i = 0, ilnbuf = 0 ; i < s ; i++) {
			non_nul(240, hestvalide(lnp1[i]), lnp1[i], ph2, lnp1[i]) ;
			lnbuf[ilnbuf] = htoc[lnp1[i]] ;
			ilnbuf++ ;
/* converts all hydrogens into carbons through HSQC data */
		}
		p1 = -lnpos() ;
/* creates a new atom list */
		lnp1 = ln[-p1] ;
		s = lnsize[-p1] ;
	} else {
		for(i = 0 ; i < s ; i++) {
			non_nul(251, estvalide(lnp1[i]), lnp1[i], p2, lnp1[i]) ;
		}
	}
	
	for(i = 0 ; (i < s) && (lnp1[i] != p2) ; i++)
		;
/* looking for p2 in the list indexed by p1 */
	if (i != s) {
/* p2 is in the p1 list and it must be removed from it.
If the list contains 2 elements, the correlation is a regular one,
assuming there is no list with a single value */
		if(o == 22) {
			if(at[lnp1[i]].mult != 2) {
				non_nul(241, FALSE, ctoh[lnp1[i]][0], ph2, p2) ;
			} else {
				ok = FALSE ;
			}
		}
		if (s > 2) {
/* p2 still correlates with a list */
			for(i = 0, ilnbuf = 0 ; i < s ; i++) {
				if (lnp1[i] != p2) {
					lnbuf[ilnbuf] = lnp1[i] ;
					ilnbuf++ ;
				}
			}
			p1 = -lnpos() ;
/* creates a new list */
		} else {
			p1 = (lnp1[0] == p2) ? lnp1[1] : lnp1[0] ;
/* p1 is the single remaining atom from the (p1 p2) list */
			flcorrvar-- ;
		}
	}
	
	if(ok) {
		newcorr(p1, p2, ph1, ph2, p3, p4, type) ;
		updatecorr(icorr) ;
		valcorr[icorr++] = TRUE ;
	}
}


void	start4(void)
/*
 * real storage of correlations in atom properties
 */
{
	int i ;
	
	for(i = 0 ; i < icorr ; i++) {
		if(valcorr[i] && propcorr[i].at1 > 0) {
			atomcorr(i) ;
/* storage of non variant correlations */
		}
	}
	for(i = 0 ; i < icorr ; i++) {
		if(valcorr[i] && propcorr[i].at1 < 0) {
			atomcorrvar(i) ;
/* storage of variant correlations */
		}
	}
}


void	atomcorr(int corr)
/*
 * updates atoms correlation data for non variant correlations
 */
{
	int	n, p1, p2 ;
	struct lcorr *pc ;
	struct atom *a ;
	
	pc = &propcorr[corr] ;
	p1 = pc->at1 ;
	p2 = pc->at2 ;
	a = &at[p1] ;
	n = a->nc ;
	non_nul(292, n - MAX_NC, MAX_NC) ;
	a->ex[n] = corr ;
	a->other[n] = p2 ;
	a->nctot = a->nc = n + 1 ;
	a = &at[p2] ;
	n = a->nc ;
	non_nul(292, n - MAX_NC, MAX_NC) ;
	a->ex[n] = corr ;
	a->other[n] = p1 ;
	a->nctot = a->nc = n + 1 ;
}


void	atomcorrvar(int corr)
/*
 * updates atoms correlation data for variant correlations
 */
{
	int	j, n, s, p1, p2 ;
	struct lcorr *pc ;
	struct atom *a ;
	
	pc = &propcorr[corr] ;
	p1 = pc->at1 ;
	p2 = pc->at2 ;
	
	s = lnsize[-p1] ;
	for (j = 0 ; j < s ; j++) {
		at[ln[-p1][j]].ingrp = TRUE ;
	}
/* atoms from p1 are in a group */
	a = &at[p2] ;
	n = a->nc ;
	non_nul(292, n - MAX_NC, MAX_NC) ;
	a->ex[n] = corr ;
	a->other[n] = p1 ;
	a->corrgrp++;
/* one more correlation of p2 with a group */
	a->nctot = a->nc = n + 1 ;
/* one more valid correlation of p2 */
}


void	s_lier(int p1, int p2, int fl)
/* 
 * initialises bonds. When called with fl being true, p1 and p2 must not
 * be previously bound 
 */
{
	int	n, paslie ;
	struct atom *a ;

	non_nul(280, estvalide(p1), p1, p2, p1) ;
	non_nul(280, estvalide(p2), p1, p2, p2) ;
	paslie = !estlie(p1, p2) ;
	if (fl) {
		non_nul(282, paslie, p1, p2) ;
	}
	if (paslie) {
/* does the job */
		a = at + p1 ;
		n = a->nlia ;
		a->lia[n] = p2 ;
		n++;
		non_nul(283, n <= a->nlmax, p1, p2, p1, a->nlmax) ;
		a->nlia = n ;

		a = at + p2 ;
		n = a->nlia ;
		a->lia[n] = p1 ;
		n++;
		non_nul(283, n <= a->nlmax, p1, p2, p2, a->nlmax) ;
		a->nlia = n ;

		lier_prop(p1, p2);
		lier_prop(p2, p1);
/* updates property information upon bonds formation */
		non_nul(420, setcnx(p1, p2, FALSE), p1, p2) ; 
/* updates connectivity indexing */
	}
}


int	s_carb(int x)
{ 
	return at[x].element == 0 ; 
}


int	s_hete(int x)
{ 
	return at[x].element != 0 ; 
}


int	s_sp3(int x)
{ 
	return at[x].hybrid == 0 ; 
}


int	s_sp2(int x)
{ 
	return at[x].hybrid == 1 ; 
}


int	s_sp(int x)
{ 
	return at[x].hybrid == 2 ; 
}


int s_char(int x)
{
	return at[x].charge != 0 ;
}


int s_charpos(int x)
{
	return at[x].charge > 0 ;
}


int s_charneg(int x)
{
	return at[x].charge < 0 ;
}


int	s_quat(int x)
{ 
	return s_carb(x) && at[x].mult == 0 ; 
}


int	s_ch(int x)
{ 
	return s_carb(x) && at[x].mult == 1 ; 
}


int	s_ch2(int x)
{ 
	return s_carb(x) && at[x].mult == 2 ; 
}


int	s_ch3(int x)
{ 
	return s_carb(x) && at[x].mult == 3 ; 
}


int	s_full(int x)
/* x is not used, I know... */
{ 
	return TRUE ; 
}


int s_elem(int x, int y)
{
	return at[x].element == y ;
}


int	s_ge(int x, int y)
{ 
	return x >= y ; 
}


int	s_le(int x, int y)
{ 
	return x <= y ; 
}


int	s_gt(int x, int y)
{ 
	return x > y ; 
}


int	s_lt(int x, int y)
{ 
	return x < y ; 
}


int	s_ou(int x, int y)
{ 
	return x || y ; 
}


int	s_et(int x, int y)
{ 
	return x && y ; 
}


int	s_diff(int x, int y)
{ 
	return x ? (!y) : x ; 
}


void	s_prop(int p, int (*f)(int))
/*
 * builds the list number p according to the atom's property f.
 */
{
	int	cardinal, x, r ;

	cardinal = 0 ;
	for (x = 1 ; x <= max_atomes ; x++) {
		r = at[x].utile && (*f)(x) ;
		list[p][x] = r ;
		if (r) {
			cardinal++;
		}
	}
	non_nul(261, cardinal, p) ;
	/*	printf("** liste %d **", p) ; check_list(p) ; printf("\n") ; */
	deflist[p] = TRUE ;
}


void	s_cmp(int p1, int p2, int (*f)(int, int))
/* 
 * selects atoms according to a comparison function f 
 */
{
	int	cardinal, x, r ;

	cardinal = 0 ;
	for (x = 1 ; x <= max_atomes ; x++) {
		r = at[x].utile && (*f)(x, p2) ;
		list[p1][x] = r ;
		if (r) {
			cardinal++;
		}
	}
	non_nul(261, cardinal, p1) ;
	/*  printf("** liste %d **", p1) ; check_list(p1) ; printf("\n") ; */
	deflist[p1] = TRUE ;
}


void	s_ens(int p1, int p2, int p3, int (*f)(int, int))
/*
 * operates sets function f. l1 and l2 are 2 lists used in order
 * to store sets made of a single element
 */
{
	int	l1[MAX_ATOMES], l2[MAX_ATOMES] ;
	int	*a1, *a2, x, r, cardinal ;

	a1 = s_adr(p1, l1) ;
	a2 = s_adr(p2, l2) ;
	cardinal = 0 ;
	for (x = 1 ; x <= max_atomes ; x++) {
		r = at[x].utile && (*f)(a1[x], a2[x]) ;
		list[p3][x] = r ;
		if (r) 
			cardinal++;
	}
	non_nul(261, cardinal, p3) ;
	/*	printf("** liste %d **", p3) ; check_list(p3) ; printf("\n") ; */
	deflist[p3] = TRUE ;
}


int	*s_adr(int p, int *al)
/* returns the address of a set that can be a single element, too */
{
	int	i ;

	if (p > 0) {
/* generates a set with a single element : p */
		non_nul(320, estvalide(p), p) ;
		for (i = 1 ; i <= max_atomes ; i++) 
			al[i] = FALSE ;
		al[p] = TRUE ;
		return al ;
	} else {
		p = -p ;
/* p is a list reference */
		non_nul(321, deflist[p], p) ;
		return list[p] ;
	}
}


void	s_useprop(int x, int n, int l, int c)
/*
 * atom x has exactly n neighbours in list l if n is not 0 and c equals -1
 * atom x has all its neighbours in list if n equals 0 and c equals -1
 * atom x has n (n > 0) or less neighbours if c equals '-'
 * atom x has n (n > 0) or more neighbours if c equals '+'
 * atom x has no neighbour in list l if n equals 0 and c equals '0'
 */
{
	int	i, max, o, eo ;
	struct atom *pa ;

/* printf("x=%d  n=%d  l=%d  c=%d\n", x, n, l, c) ; */
	pa = &at[x] ;
	max = pa->nlmax ;
//	non_nul(331, max >= n, n, x, l, max) ;
	if ((c == -1) || (c == '+')) {
		non_nul(331, max >= n, n, x, l, max) ;
	}
/* n is in the 0..max range */
	o = FALSE ;
/* becomes TRUE when at[x].occ is forced to be >= 0 */
	eo = FALSE ;
/* becomes TRUE when at[x].extocc is forced to be >= 0 */
	switch (c) {
	case -1 :
/* as before 3.3.10 */
		if (n == 0) {
/* 0 neighbours means max neighbours when the PROP command has only 3 parameters */
			n = max ;
			eo = TRUE ;
/* n is equal to max, no need to impose occ >= 0 */
		} else if (n == max) {
/* n is equal to max, no need to impose occ >= 0 */
			eo = TRUE ;
		} else {
/* n is the exact number of neighbours, requires both occ >= 0 and extocc >= 0 */
			o = TRUE ;
			eo = TRUE ;
		}
		break ;

//	case '0' :
//		non_nul(332, n == 0, n, x, l) ;
/* only possible if n equals 0, really means 0 neighbours */
//		*po = TRUE ;
/* requires only occ >= 0 */
//		break ;

	case '+' :
		if (n != 0) {
			eo = TRUE ;
/* requires only extocc >= 0 */
		}
		break ;
	case '-' :
		if (n != max) {
/* max of less neighbors does not require any control */
			o = TRUE ;
/* requires only occ >= 0 */
		}
		break ;
	}
	
	if (o || eo) {
/* create new property of atom x */
		i = pa->nlpr ;
/* i is the present number of property lists */
		non_nul(330, i < MAX_NLPR, x, MAX_NLPR) ;
/* check that maximum i has not already been reached */
		pa->lpr[i] = l ;
		pa->occ[i] = n ;
/* number of atoms inside of the property list */
		pa->extocc[i] = max - n ;
/* number of atoms out of the property list */
		pa->flocc[i] = o ;
/* true if occ must be checked when property is used */
		pa->flextocc[i] = eo ;
/* true if extocc must be checked when property is used */
		pa->nlpr++;
/* ready for the next (if any) property list of atom x */
	}
}


void	lier_prop(int p1, int p2)
/*
 * property data management upon bond formation between p1 and p2
 * looks for p2 in all property lists of p1
 */
{
	int	i, n, *a, inside ;
	struct atom *pa ;

	pa = &at[p1] ;
	n = pa->nlpr ;
	for (i = 0 ; i < n ; i++) {
		inside = list[pa->lpr[i]][p2] ;
		a = NULL ;
		if (pa->flocc[i] && inside) {
/* was "if (pa->flocc && inside) {" before 3.4.9 and was a bug */
			a = &pa->occ[i] ;
		}
		if (pa->flextocc[i] && !inside) {
/* was "if (pa->flextocc && inside) {" before 3.4.9 and was a bug */
			a = &pa->extocc[i] ;
		}
		if(a != NULL) {
			non_nul(340, (*a) > 0, p1, p2, p1) ;
			(*a)-- ;
		}
	}
}


int	estvalide(int x)
/*
 * returns true if x ( x >= 0 )is a defined atom
 */
{
	return x < MAX_ATOMES && at[x].utile ;
}


int	hestvalide(int x)
/*
 * returns true if x ( x >= 0 ) is a defined hydrogen atom
 */
{
	return x < MAX_ATOMES && htoc[x] ;
}


int	estlie(int x, int y)
/*
 * returns true if x and y are bound
 */
{
	int	i, n ;

	n = at[x].nlia ;
	if (!n || !at[y].nlia) {
		return FALSE ;
	}
	for (i = 0 ; (i < n) && (at[x].lia[i] != y) ; i++) 
		;
/* y is x's ith neighbour */
	return (i != n) ;
}


int	estlielie(int p1, int p2)
/*
 * returns the index of the common neighbor if p1 and p2 share
 * a common neighbour, returns FALSE otherwise
 */
{
	int	i, n, x ;

	n = at[p1].nlia ;
	if (!n || !at[p2].nlia) {
		return FALSE ;
	}
/* p1 or p2 has no bond */
	for (i = 0 ; (i < n) && ((x = at[p1].lia[i]), !estlie(x, p2)) ; i++) 
		;
/* x is the hypothetic common neighbour */
	return i == n ? FALSE : x ;
}


void	about()
/*
 * displays what has been understood from the input data.
 * If bonds are established using either the BOND, the COSY or the HMBC commands,
 * correlation and property information are modified accordingly.
 */
{
	struct atom *a ;
	int	i, x, n ;
	struct lcorr *pc ;
	int c11 ;
	int nmin, nmax ;
	
	if (nass && fldeff) {
		printf("Warning: printing of atom assignment is disabled\n") ;
	}
	for (x = 1 ; x <= max_atomes ; x++) {
		a = &at[x] ;
		if (!a->utile) {
			continue ;
		}
		printf("atom %d, element %s, hybridization sp%d, charge %d, ",
			x, nucl[a->element].sym, 3 - a->hybrid, a->charge) ;
		n = a->mult ;
		printf("bears %d hydrogen%s \n",  n, (n > 1) ? "s" : " ") ;
		n = a->ics ;
		if(n != -1) {
			printf("\tchemical shift: %.2f ppm \n", reels[a->ics]) ;
		}
		n = a->nlpr ;
		if (n) {
			printf("\towns %d property list%s \n", n, (n > 1) ? "s" : " ") ;
			for (i = 0 ; i < n ; i++) {
				int o, eo, fo, feo ;
				o = a->occ[i] ;
				eo = a->extocc[i] ;
				fo = a->flocc[i] ;
				feo = a->flextocc[i] ;
				printf("\t\tList : ") ;
				check_list(a->lpr[i]) ;
				if(fo) {
						printf(" (In: %d) ", o) ;
				}
				if(feo) {
						printf(" (Out: %d) ", eo) ;
				}
				printf("\n") ;
			}
		}
		n = a->nlia ;
		if (n) {
			printf("\towns %d bond%s with", n, (n > 1) ? "s" : " ") ;
			for (i = 0 ; i < n ; i++) {
				printf("%3d", a->lia[i]) ;
			}
			printf("\n") ;
		}
		n = ctoh[x][0] ;
		if (n) {
			printf("\tbears hydrogen %d\n", n) ;
		}
		n = ctoh[x][1] ;
		if (n) {
			printf("\tand bears hydrogen %d\n", n) ;
		}
		n = a->ihcs[0] ;
		if(n != -1) {
			printf("\thydrogen chemical shift: %.2f ppm", reels[n]) ;
			n = a->ihcs[1] ;
			if(n != -1) {
				printf(" and %.2f ppm\n", reels[n]) ;
			} else {
				printf("\n") ;
			}
		}
		n = a->nctot ;
		if (n) {
			int	j ;

			for (j = 0 ; j < n &&  !((a->other[j] > 0) && valcorr[a->ex[j]]) ; j++) 
				;
			if (j != n) {
				printf("\tcorrelates with : ") ;
				for (i = 0 ; i < n ; i++) {
					if ((a->other[i] > 0) && valcorr[a->ex[i]]) {
						pc = &propcorr[a->ex[i]] ;
						printf("%3d [%d-%d]", a->other[i], pc->min, pc->max) ;
					}
				}
				printf("\n") ;
			}
			if (a->corrgrp) {
				for (i = 0 ; i < n ; i++) {
					if ((a->other[i] < 0) && valcorr[a->ex[i]]) {
						pc = &propcorr[a->ex[i]] ;
						printf("\tcorrelates with an atom among : ") ;
						wrln(-a->other[i]) ;
						printf(" [%d-%d]", pc->min, pc->max) ;
						printf("\n") ;
					}
				}
			}
		}
		n = a->equstat ;
		if (n >= 0) {
			printf("\telement %d from equivalence class %d\n", a->iequ, n) ;
		}
		if (!fldeff) {
			n = a->asgn ;
			if (n > 0) {
				printf("\tassignment of sub-atom %d\n",  sub[n].ss_ref) ;
			}
		}
	}
	
	printf("Input file HMBC and COSY correlations:\n") ;
	for(i = 0 ; i < icorr ; i++) {
		pc = &propcorr[i] ;
		if((pc->origin < 2) && (pc->at1 > 0)) {
			printf("\t%2d %2d [%d-%d]\n", pc->inputat1, pc->inputat2, pc->inputmin, pc->inputmax) ;
		}
	}
	for(i = 0 ; i < icorr ; i++) {
		pc = &propcorr[i] ;
		if((pc->origin < 2) && (pc->at1 < 0)) {
			printf("\t( ") ;
			wrln(-pc->inputat1) ;
			printf(") %2d [%d-%d]\n", pc->inputat2, pc->inputmin, pc->inputmax) ;
		}
	}
	printf("Effective HMBC and COSY correlations:\n") ;
	for(i = 0 ; i < icorr ; i++) {
		pc = &propcorr[i] ;
		if(valcorr[i] && (pc->at1 > 0)) {
			nmin = pc->type ? 2 : 1 ;
			nmax = pc->type ? (pc->max ? 2 : 0) : (pc->max ? 1 : 0) ;
			printf("\t%2d %2d [%d-%d]\n", pc->at1, pc->at2, pc->min + nmin, pc->max + nmax) ;
		}
	}
	for(i = 0 ; i < icorr ; i++) {
		pc = &propcorr[i] ;
		if(valcorr[i] && (pc->at1 < 0)) {
			nmin = pc->type ? 2 : 1 ;
			nmax = pc->type ? (pc->max ? 2 : 0) : (pc->max ? 1 : 0) ;
			printf("\t( ") ;
			wrln(-pc->at1) ;
			printf(") %2d [%d-%d]\n", pc->at2, pc->min + nmin, pc->max + nmax) ;
		}
	}
	if (nass || fldeff) { 
		printf("sub-structure search ") ;
		if (flsubsgood) {
			printf("active, keeps successful results\n") ;
		} else if (flsubsbad) {
			printf("active, keeps unsuccessful results\n") ;
		} else {
			printf("inactive\n") ;
		}
		sub_about() ;
	}
}


void	wrln(int i)
/*
 * writes the ith numerical list. No \n issued
 */
{
	int	j ;

	for (j = 0 ; j < lnsize[i] ; j++) {
		printf("%d ", ln[i][j]) ;
	}
}


void	wrln2(int y)
/*
 * writes either a numerical list or an integer. 
 */
{
	if (y >= 0) {
		printf("%d\n", y) ;
	}
	else {
		wrln(-y) ;
		printf("\n") ;
	}
}


int	stoss(int i)
/*
 * gives subatom's internal reference number from its external reference
 * number ( the one following [Ss] in the input file ). Returns 0
 * if the subatom is unknown.
 */
{
	int	j ;

	for (j = 1 ; j <= nass && sub[j].ss_ref != i ; j++) 
		;
	return (j > nass) ? 0 : j ;
}

