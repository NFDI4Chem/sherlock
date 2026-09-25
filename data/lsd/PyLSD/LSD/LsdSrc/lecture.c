/*
    file: LsdSrc/lecture.c
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
#include <stdarg.h>
#include "defs.h"

void	lecture(void) ;
void	remcomment(void) ;
void	lirecom(int ic) ;
void	lireliste(void) ;
void	non_nul(int errnum, int v, ...) ;
int	z(char *noyau) ;
void	check_lecture(void) ;
void	check_list(int l)  ;
int	lnpos(void) ;

extern void	classe(int n, int t[]) ; /* lsd.c */
extern void	myexit(int i) ; /* lsd.c */

extern int	lnsize[MAX_LN] ;
extern int	ln[MAX_LN][MAX_LNBUF] ;
extern int	list[MAX_LIST][MAX_ATOMES] ;
extern FILE 	*fic1 ;
extern dynstring	*curstring ; /* lsd.c */
extern dyntext	*strtab ; /* lsd.c */
extern dyntext	*pathtab ; /* lsd.c */

struct comtab coms[MAX_COM] = {
/*
 * field 1 : mnemonic
 * field 2 : number of arguments, -1 if it is a list (of unknown length)
 * field 3 : type of the arguments, see lirecom
 * field 4 : analysis priority order
 */
/*	{ "VALE", 4, "AIIR", 0 }, */ /* 0 */
	{ "VALE", 3, "AIR", 0 }, /* 0 */
	{ "MULT", 5, "IAIIZ", 1 }, /* 1 */
	{ "LIST", -1, "L", 3 }, /* 2 */
	{ "CARB", 1, "L", 3 }, /* 3 */
	{ "HETE", 1, "L", 3 }, /* 4 */
	{ "SP3 ", 1, "L", 3 }, /* 5 */
	{ "SP2 ", 1, "L", 3 }, /* 6 */
	{ "QUAT", 1, "L", 3 }, /* 7 */
	{ "CH  ", 1, "L", 3 }, /* 8 */
	{ "CH2 ", 1, "L", 3 }, /* 9 */
	{ "CH3 ", 1, "L", 3 }, /* 10 */
	{ "FULL", 1, "L", 3 }, /* 11 */
	{ "GREQ", 2, "LI", 3 }, /* 12 */
	{ "LEEQ", 2, "LI", 3 }, /* 13 */
	{ "GRTH", 2, "LI", 3 }, /* 14 */
	{ "LETH", 2, "LI", 3 }, /* 15 */
	{ "UNIO", 3, "BBL", 3 }, /* 16 */
	{ "INTE", 3, "BBL", 3 }, /* 17 */
	{ "DIFF", 3, "BBL", 3 }, /* 18 */
	{ "PROP", 4, "BILH", 4 }, /* 19 */
	{ "BOND", 2, "II", 5 }, /* 20 */
	{ "HMQC", 2, "II", 6 }, /* 21 */
	{ "COSY", 4, "VIOO", 8 }, /* 22 */
	{ "HMBC", 4, "VIOO", 9 }, /* 23 */
	{ "ENTR", 1, "I", 0 }, /* 24 */
	{ "HIST", 1, "I", 0 }, /* 25 */
	{ "DISP", 1, "I", 0 }, /* 26 */
	{ "VERB", 1, "I", 0 }, /* 27 */
	{ "PART", 1, "I", 0 }, /* 28 */
	{ "STEP", 1, "I", 0 }, /* 29 */
	{ "WORK", 1, "I", 0 }, /* 30 */
	{ "MLEV", 1, "I", 0 }, /* 31 */
	{ "SSTR", 4, "SAVV", 10 }, /* 32 */
	{ "ASGN", 2, "SI", 11 }, /* 33 */
	{ "LINK", 2, "SS", 12 }, /* 34 */
	{ "DUPL", 1, "I", 0 }, /* 35 */
	{ "SUBS", 1, "T", 0 }, /* 36 */
	{ "ELIM", 2, "II", 0 }, /* 37 */
	{ "FILT", 1, "I", 0 }, /* 38 */
	{ "DEFF", 2, "FC", 14 }, /* 39 */
	{ "FEXP", 1, "C", 15 }, /* 40 */
	{ "CNTD", 1, "I", 0 }, /* 41 */
	{ "MAXS", 1, "I", 0 }, /* 42 */
	{ "MAXT", 1, "I", 0 }, /* 43 */
	{ "HSQC", 2, "II", 6 }, /* 44 */
	{ "CCLA", 1, "I", 0 }, /* 45 */
	{ "PATH", 1, "C", 13 }, /* 46 */
	{ "SKEL", 2, "FC", 14 }, /* 47 */
	{ "COUF", 1, "P", 0 }, /* 48 */
	{ "STOF", 1, "P", 0 }, /* 49 */
	{ "SP  ", 1, "L", 3 }, /* 50 */
	{ "CHAR", 1, "L", 3 }, /* 51 */
	{ "CPOS", 1, "L", 3 }, /* 52 */
	{ "CNEG", 1, "L", 3 }, /* 53 */
	{ "ELEM", 2, "LA", 3 }, /* 54 */
	{ "SHIX", 2, "IR", 2 }, /* 55 */
	{ "SHIH", 2, "IR", 7 }, /* 56 */
	{ "BRUL", 1, "I", 0 } /* 57 */
/* 58 */
};


struct noy nucl[MAX_NUCL] = {
	{ "C", {3, 4, 3, 0}, 12.0, 0, "C" },
	{ "O", {1, 2, 3, 0}, 16.0, 0, "O" },
	{ "N", {2, 3, 4, 3}, 14.0, 0, "N" },
	{ "N5", {0, 5, 0, 0}, 14.0, 0, "N" },
	{ "S", {1, 2, 3, 4}, 32.1, 0, "S" },
	{ "S4", {0, 4, 0, 0}, 32.1, 0, "S" },
	{ "S6", {0, 6, 0, 0}, 32.1, 0, "S" },
	{ "F", {0, 1, 0, 0}, 19.0, 0, "F" },
	{ "Cl", {0, 1, 0, 0}, 35.5, 0, "Cl" },
	{ "Br", {0, 1, 0, 0}, 79.9, 0, "Br" },
	{ "I", {0, 1, 0, 0}, 126.9, 0, "I" },
	{ "Si", {3, 4, 3, 0}, 28.1, 0, "Si" },
	{ "P", {2, 3, 4, 0}, 31.0, 0, "P" },
	{ "P5", {0, 5, 0, 0}, 31.0, 0, "P" },
	{ "B", {4, 3, 0, 0}, 10.8, 0, "B" },
	{ "X", {0, 0, 0, 0}, -1.0, 0, "X" },
	{ "A", {-1, -1, -1, -1}, 1.0, 0, "A" }
};


struct crt cartes[MAX_CARTES] ;

extern float reels[MAX_REELS] ; /* lsd.c */
int	ireels = 0 ;
float 	*preels = reels ;

int	lignes_lues ;
int	icarte ;
int	lnbuf[MAX_LNBUF] ;

int	ilnbuf ;
int	iln = 1;

void	lecture(void)
{
	char	c[5] ;
	int	ic, go, n ;

	lignes_lues = 0 ;
	icarte = 0 ;
	go = TRUE ;
	while (go) {
		remcomment() ;
/* removes comment(s) if any */
		go = !feof(fic1) ;
/* eof without an EXIT command should be possible */
		if (!go) {
			break ;
		}
		non_nul(100, fscanf(fic1, " %4[A-Z 23]", c)) ;
		non_nul(101, strlen(c) == 4, c) ;
/* read mnemonic */
		go = strcmp(c, "EXIT") ;
		if (!go) {
			break ;
		}
/* EXIT not encountered */
		lignes_lues++;
/* counts command lines */
		for (ic = 0 ; ic < MAX_COM && strcmp(coms[ic].op, c) ; ic++)
			;
		non_nul(102, MAX_COM - ic, c) ;
/* the mnemonic exists */
		lirecom(ic) ;
/* reads arguments following the mnemomic */
	}
	/*check_lecture() ; */
}


void	remcomment()
/*
 * a comment starts at a ; and ends up at the end of the line.
 */
{
	char	s[2] ;

	while (!feof(fic1) && fscanf(fic1, " %1[;]", s)) {
		while (!feof(fic1) && (getc(fic1)) != '\n') 
			;
	}
}


void	lirecom(int ic)
/* reads arguments of current command */
{
	int	i, nf, fl, v, h ;
	char	f, noyau[3], l[2] ;

	cartes[icarte].opnum = ic ;
	cartes[icarte].ln = lignes_lues ;
	nf = coms[ic].fields ;
/* expected arguments number */
	fl = nf < 0 ;
	if (fl) {
		nf = -nf ;
	}
/* special case of enumerated lists */
	for (i = 0 ; i < nf ; i++) {
		f = coms[ic].forme[i] ;
/* expecting argument of type f. its value will be stored in v */
		switch (f) {
		case 'A' :
/* symbol of atom requested */
			non_nul(110, fscanf(fic1, "%s", noyau)) ;
			v = z(noyau) ;
/* v is the nucleus number */
			break ;
		case 'B' :
/* either an atom index or [Ll]stuck with an integer
 * that is either an atom index or a logical list index */
			h = fscanf(fic1, " %1[Ll]", l) ;
/* h is true if it is a list index */
			non_nul(111, fscanf(fic1, "%d", &v)) ;
/* v is either an atom index or a list index */
			if (h) {
				non_nul(112, v < MAX_LIST, v, MAX_LIST-1) ;
				non_nul(115, v >= 0, v) ;
/* list indexes start at 0 */
				v = -v ; 
/* the list index is stored as its opposite */
			} else {
				non_nul(123, v > 0, v) ;
/* atom indexes are stored as themselves, start at 1 */
			}
			break ;
		case 'L' :
/* a list reference */
			non_nul(113, fscanf(fic1, " %1[Ll]", l)) ;
			non_nul(114, fscanf(fic1, "%d", &v)) ;
			non_nul(112, v < MAX_LIST, v, MAX_LIST-1) ;
			non_nul(115, v >= 0, v) ;
			break ;
		case 'I' :
/* an integer that is positive or zero */
			non_nul(116, fscanf(fic1, "%d", &v)) ;
			non_nul(124, v >= 0, v, i) ;
			break ;
		case 'T' :
/* an integer that can only be -1 0 1 */
			non_nul(116, fscanf(fic1, "%d", &v)) ;
			non_nul(128, v >= -1 && v <= 1, v, i) ;
			break ;
		case 'V' :
/* a numeric list or a single integer. a numeric list is enclosed 
 * between parentheses */
			h = fscanf(fic1, " %1[(]", l) ;
			if (h) {
/* this is a list */
				ilnbuf = 0 ;
/* index on numeric list buffer */
				do {
					non_nul(118, fscanf(fic1, "%d", &v)) ;
					lnbuf[ilnbuf++] = v ;
					non_nul(125, v >= 0, v) ;

/* put integer into the buffer */
				} while ((!fscanf(fic1, " %1[)]", l)) && (ilnbuf < MAX_LNBUF)) ;
/* while no ) is encountered */
				non_nul(119, MAX_LNBUF - ilnbuf, MAX_LNBUF) ;
				v = -lnpos() ;
/* lnpos stores the content of the buffer to its definive location
 * lnpos returns the numeric list number, stored as its opposite */
			} else {
				non_nul(120, fscanf(fic1, "%d", &v)) ;
				non_nul(124, v >= 0, v, i) ;
			}
/* single integer */
			break ;
		case 'S' :
/* subatome reference, starts with [Ss] */
			non_nul(121, fscanf(fic1, " %1[Ss]", l)) ;
			non_nul(122, fscanf(fic1, "%d", &v)) ;
/* range checking for v is achieved in sub.c */
			break ;
		case 'R' :
/* read an optional real number, either atom mass or chemical shift */
			non_nul(146, fscanf(fic1, "%f", preels)) ;
			v = ireels ;
			preels++ ;
			ireels++ ;
			non_nul(145, MAX_REELS - ireels, MAX_REELS) ;
			break ;
		case 'F' :
/* fragment (substructure) reference, starts with [F1] */
			non_nul(150, fscanf(fic1, " %1[Ff]", l)) ;
			non_nul(151, fscanf(fic1, "%d", &v)) ;
			non_nul(152, v > 0, v) ;
			break ;
		case 'C' :
/* character string */
		case 'P' :
/* file path */
			non_nul(160, fscanf(fic1, " %1[\"]", l)) ;
			dynstring_reset(curstring) ;
			v = 0 ;
			while(!feof(fic1) && (h = fgetc(fic1)) != '\"') {
				dynstring_pushc(curstring, h) ;
				v++ ;
			}
			non_nul(161, h == '\"') ;
			non_nul(162, v) ;
			if (f == 'C') {
				v = dyntext_push(strtab, curstring, 0, NULL) ;
			} else {
				v = dyntext_push(pathtab, curstring, 0, NULL) ;
			}
			non_nul(163, v >= 0) ;
			break ;
		case 'H' :
/* + or - */
			h = fscanf(fic1, " %[+-]", l) ;
			v = (h == 1) ? (int)l[0] : -1 ;
			break ;
		case 'O' :
/* optional parameter */
			v = -1 ;
			h = fscanf(fic1, " %d", &v) ;
			if(h == 1) {
				non_nul(170, v >= 0, v) ;
			}
			break ;
		case 'Z' :
/* optional integer that can only be -1, 0, 1 or 2 */
			v = 0 ;
			if (fscanf(fic1, "%d", &v)) {
				non_nul(129, v >= -1 && v <= 2, v) ;
			}
			break ;
		}
		cartes[icarte].pars[i] = v ;
/* v is stored */
	}
	if (fl) {
/* there is still the content of the list to read */
		lireliste() ;
	}
	icarte++;
/* next command */
	non_nul(117, MAX_CARTES - icarte, MAX_CARTES) ;
}


int	lnpos(void)
{
	int	i, j ;

	classe(ilnbuf, lnbuf) ;
/* ordering numeric values in lnbuf */
	for (j = 0 ; j < ilnbuf-1 && lnbuf[j] != lnbuf[j+1] ; j++)
		;
	non_nul(140, j == ilnbuf-1, lnbuf[j]) ;
/* searching for duplicated values */
	for (i = 1 ; i < iln ; i++) {
/* looking for the equality of lnbuf and a previously entered list */
/* remember that iln starts at 1 */
		if (ilnbuf != lnsize[i]) {
			continue ;
/* they do not have the same number of elements */
		}
		for (j = 0 ; j < ilnbuf && lnbuf[j] == ln[i][j] ; j++) 
			;
		if (j == ilnbuf) {
			return i ;
/* ln[i] is identical to lnbuf, return i */
		}
	}
	non_nul(141, MAX_LN - iln, MAX_LN) ;
	for (j = 0 ; j < ilnbuf ; j++) {
		ln[iln][j] = lnbuf[j] ;
	}
	lnsize[iln] = ilnbuf ;
/* stores lnbuf as ln[iln] and stores its size */
	i = iln ;
/* current list is indexed by iln */
	iln++ ;
	return i ;
}


void	lireliste(void)
{
	int	i, il, a ;

	il = cartes[icarte].pars[0] ;
	for (i = 0 ; i < MAX_ATOMES ; i++) {
		list[il][i] = FALSE ;
	}
/* resets the list */
	while (fscanf(fic1, "%d", &a)) {
		non_nul(126, a < MAX_ATOMES, a, MAX_ATOMES-1) ;
		non_nul(127, a > 0, a) ;
		list[il][a] = TRUE ;
/* reads the list */
	}
}


void	non_nul(int errnum, int v, ...)
/* 
 * prints a message if the value v is 0 
 */
{
#define MAX_MSG 136
	va_list l ;
	int	i ;
	static struct msg {
		int	errco ;
		char	*errmsg ;
	} msgs[MAX_MSG] = {
		{ 100, "Bad first character in command name." }, /*  */
		{ 101, "Bad character in command name. "
"Reading stopped after \"%s\"." }, /*  */
		{ 102, "Unknown command name: %s" }, /*  */
		{ 110, "Cannot read atom symbol." }, /*  */
		{ 111, "Cannot read atom or atom list index." }, /*  */
		{ 112, "List index %d is too high.\n"
"Maximum list index is %d.\n"
"Try with a lower list index or recompile with a higher MAX_LIST." }, /*  */
		{ 115, "List index cannot be negative (is %d)." }, /*  */
		{ 113, "L expected." }, /*  */
		{ 114, "Cannot read list index." }, /*  */
		{ 123, "Atom indexes start at 1 (is %d)." }, /*  */
		{ 116, "Cannot read integer value." }, /*  */
		{ 124, "Integer value cannot be negative (is %d) field %d." }, /*  */
		{ 118, "Cannot read integer value in a (...) list." }, /*  */
		{ 125, "Value in a (...) list cannot be negative (is %d)." }, /*  */
		{ 119, "The size of a (...) list is limited to %d elements.\n"
"Carefully consider recompiling with a higher MAX_LNBUF." }, /*  */
		{ 120, "Cannot read integer value." }, /*  */
		{ 117, "The number of commands is limited to %d.\n"
"Recompile with a higher MAX_CARTES." }, /*  */
		{ 121, "A sub-structure atom reference starts with \"S\"." }, /*  */
		{ 122, "Cannot read sub-structure atom index after \"S\"." }, /*  */
		{ 140, "Duplicated value (is %d) in a (...) list." }, /*  */
		{ 141, "The number of (...) lists is limited to %d.\n"
"Carefully consider recompiling with a higher MAX_LN." }, /*  */
		{ 126, "Atom index %d is higher than the allowed maximum (%d)." }, /*  */
		{ 127, "Atom indexes start at 1 (is %d)." }, /*  */
		{ 128, "Integer value should be -1, 0, or 1 (is %d)." }, /*  */
		{ 129, "Integer value should be -1, 0, 1 or 2 (is %d)." }, /*  */
		{ 131, "Unknown atom symbol: %s." }, /*  */
		{ 145, "The number of real values is limited to %d.\n"
"Recompile with a higher MAX_REELS." }, /*  */
		{ 146, "Cannot read real value." }, /*  */
		{ 150, "F expected." }, /*  */
		{ 151, "Cannot read fragment (sub-structure) index after F." }, /*  */
		{ 152, "Fragment (sub-structure) indexes start at 1 (is %d)." }, /*  */
		{ 160, "A string must start with a \"." }, /*  */
		{ 161, "End Of File was reached during string reading." }, /*  */
		{ 162, "String cannot be empty." }, /*  */
		{ 163, "Memory allocation error.\n"
"Cannot store new string." }, /*  */
		{ 170, "Optional argument cannot be negative (is %d)." }, /*  */
		{ 199, "Try again, but this time with a molecule!" }, /*  */
		{ 200, "Valence of \"%s\" is already defined (is %d)." }, /*  */
		{ 201, "Valence of \"%s\" cannot be 0." }, /*  */
		{ 202, "Valence of \"%s\" cannot be greater than %d (is %d).\n"
"Recompile with a higher MAX_NLIA." }, /*  */
		{ 203, "Mass of \"%s\" cannot be less than 1.0 (is %.5f)." }, /*  */
		{ 204, "VALE is not compatible with \"DUPL 2\"." }, /*  */
		{ 210, "Atom index cannot be 0." }, /*  */
		{ 211, "The atom indexes are limited to %d (is %d).\n"
"Recompile with a higher MAX_ATOMES." }, /*  */
		{ 212, "Hybridization state of atom %d is %d.\n"
"It must be 1, 2 or 3 respectively for sp, sp2 and sp3 atoms." }, /*  */
		{ 213, "The charge of atom %d cannot be %d." }, /*  */
		{ 214, "Atom %d cannot be of valence %d, be sp%d hybridized, "
"and be bonded to %d hydrogen atoms." }, /*  */
		{ 215, "Atoms N5 and P5 must not be sp3 hybridized because maximum number of neighbors is 4." }, /*  */
		{ 216, "Atom S6 must be sp hybridized because maximum number of neighbors is 4." }, /*  */
		{ 260, "Atom %d is undefined but in list L%d." }, /*  */
		{ 261, "List L%d cannot be empty." }, /*  */
		{ 270, "Cannot define property of unknown atom %d." }, /*  */
		{ 271, "Cannot define property of atoms in unknown list L%d." }, /*  */
		{ 230, "Cannot set an HMQC correlation of atom %d with H-%d "
"because atom %d is not defined by a MULT command." }, /*  */
		{ 231, "Cannot set an HMQC correlation of atom %d with H-%d "
"because hydrogen atom indexes start at 1."}, /*  */
		{ 232, "Cannot set an HMQC correlation of atom %d with H-%d "
"because the atom indexes are limited to %d.\n"
"Recompile with a higher MAX_ATOMES." }, /*  */
		{ 233, "Cannot set an HMQC correlation of atom %d with H-%d "
"because atom %d already correlates with H-%d and H-%d." }, /*  */
		{ 234, "Cannot set an HMQC correlation of atom %d with H-%d "
"because atom %d already correlates with H-%d." }, /*  */
		{ 235, "Cannot set an HMQC correlation of atom %d with H-%d "
"because atom %d already correlates with H-%d and does not belong to a "
"%sH2 group." }, /*  */
		{ 240, "Cannot set a COSY correlation between H-%d and H-%d "
"because H-%d is not defined by an HMQC command." }, /*  */
		{ 241, "Correlation between H-%d and H-%d makes %d correlating with itself." }, /*  */
		{ 242, "Minimal coupling path length of a COSY correlation cannot be 2." }, /*  */
		{ 243, "Maximal coupling path length of a COSY correlation cannot be 2." }, /*  */
		{ 250, "Cannot set an HMBC correlation between %d and H-%d "
"because H-%d is not defined by an HMQC command." }, /*  */
		{ 251, "Cannot set an HMBC correlation between %d and H-%d "
"because atom %d is not defined by a MULT command." }, /*  */
		{ 252, "Correlation between %d and H-%d makes %d correlating with itself." }, /*  */
		{ 253, "Minimal coupling path length cannot be 1." }, /*  */
		{ 254, "Maximal coupling path length cannot be 1." }, /*  */
		{ 255, "Minimal coupling path length cannot be 0." }, /*  */
		{ 256, "Minimal coupling path length (%d) cannot be greater than the maximum one (%d)." }, /*  */
		{ 257, "Limit for maximal coupling path length is %d for COSY correlations "
"and %d for HMBC correlations (value provided by the ELIM command)." }, /*  */
		{ 258, "Correlation %s is already defined with conflicting coupling path lengths." }, /*  */
		{ 259, "Correlation %s is conflicting with bond %d-%d." }, /*  */
		{ 290, "Correlation %s is not valid because %d and %d are separated by two bonds." }, /*  */
		{ 291, "The number of correlations is limited to %d.\n"
"Recompile with a higher MAX_CORR." }, /*  */
		{ 292, "The number of correlations of an atom is limited to %d.\n"
"Carefully consider recompiling with a higher MAX_NC." }, /*  */
		{ 293, "The very long-range correlation %s cannot be eliminated. "
"Consider adding an ELIM command to the data file." }, /*  */
		{ 280, "Attempting to bind atoms %d and %d but %d is not defined." }, /*  */
		{ 282, "Attempting to bind already bonded atoms %d and %d." }, /*  */
		{ 283, "Attempting to bind atoms %d and %d but %d has already "
"reached the maximum number of neighbors (is %d)." }, /*  */
		{ 300, "Cannot set chemical shift of atom %d because it is not defined by a MULT command." }, /*  */
		{ 301, "Chemical shift of atom %d is already defined." }, /*  */
		{ 302, "Cannot set chemical shift of hydrogen H-%d because it is not defined by an HSQC command." }, /*  */
		{ 303, "Chemical shift of hydrogen H-%d is already defined." }, /*  */
		{ 304, "Cannot set a second chemical shift to H-%d because it does not belong to a %sH2 group." }, /*  */
		{ 320, "Attempting to use the undefined atom %d in a list "
"construction command." }, /*  */
		{ 321, "Attempting to use the undefined atom list L%d in a list "
"construction command." }, /*  */
		{ 330, "The maximum number of property lists of atom %d "
"is %d.\nRecompile with a higher MAX_NLPR." }, /*  */
//		{ 331, "Number of neighbors (%d) in PROP command of atom %d "
//"and list L%d is limited to %d." }, /*  */
		{ 331, "The requested number of neighbors (%d) in the PROP command of atom %d "
"and list L%d is not compatible with the exact number of neighbors (%d)." }, /*  */
		{ 340, "Cannot bind atoms %d and %d.\n"
"Property list violation of atom %d, from either a BOND or a COSY command." }, /*  */
		{ 350, "Sub-atom indexes start at 1 (is %d)." }, /*  */
		{ 351, "Attempting to redefine sub-atom S%d." }, /*  */
		{ 352, "Sub-atom indexes are limted to %d.\n"
"Recompile with a higher MAX_NASS." }, /*  */
		{ 353, "Hybridization state of sub-atom S%d can only be "
"sp, sp2 or sp3 (is sp%d)." }, /*  */
		{ 354, "Bad multiplicity of sub-atom S%d (is %d)." }, /*  */
		{ 360, "Sub-atom S%d is not defined by a SSTR command." }, /*  */
		{ 361, "Cannot identify sub-atom S%d to unknown atom %d." }, /*  */
		{ 362, "Cannot identify sub-atom S%d to atom %d "
"because it is already identified to atom %d." }, /*  */
		{ 370, "Too many bonds in sub-structure.\n"
"Recompile with a higher MAX_NSSL." }, /*  */
		{ 371, "Attempting to bind sub-atoms S%d and S%d "
"but S%d is not defined by a SSTR command." }, /*  */
		{ 372, "Attempting to bind sub-atoms S%d and S%d "
"but S%d has too many (is %d) bonds.\nRecompile with a higher MAX_NLIA." }, /*  */
		{ 375, "DUPL is %d but can only be 0, 1 or 2"}, /*  */
		{ 380, "Number of bonds in ELIM must be either 0 (no limit) "
"or strictly greater than 3 (is %d)." }, /*  */
		{ 390, "Odd number of sp2 atoms (is %d)." }, /*  */
		{ 391, "Odd total sum of valences (is %d)." }, /*  */
		{ 392, "Number of rings uncompatible with connected solutions (is %d)." }, /*  */
		{ 400, "Attempting to redefine fragment F%d." }, /*  */
		{ 401, "Cannot open file %s that defines fragment F%d." }, /*  */
		{ 402, "FEXP: Attempting to redefine the filter expression." }, /*  */
		{ 403, "FEXP: Cannot open expression file: %s." }, /*  */
		{ 404, "FEXP: Internal error.\n" 
"failed to create the logical expression compiler." }, /*  */
		{ 405, "FEXP: Failed to compile the logical expression: \"%s\"" }, /*  */
		{ 406, "FEXP: F%d is an unknown fragment." }, /*  */
		{ 410, "Invalid command in this fragment file: %s." }, /*  */
		{ 420, "Attempting to bind atoms %d and %d leads to "
"non-connected solutions.\nAdd a \"CNTD 0\" command if this is what is wanted." }, /*  */
		{ 430, "No directory named: %s." }, /*  */
		{ 431, "Directory name is too long:\n%s\n%d characters, maximum is %d." }, /*  */
		{ 432, "Directory cannot be opened: %s." }, /*  */
		{ 433, "File name is too long:\n%s\n%d characters, maximum is %d." }, /*  */
		{ 434, "File cannot be opened: %s." }, /*  */
		{ 435, "Fragment F%d not found: %s." }, /*  */
		{ 436, "Fragment F%d is ambiguous: %s." }, /*  */
		{ 437, "PATH directory duplicated: %s." }, /*  */
		{ 438, "PATH directories embedded: %s in %s." }, /*  */
		{ 439, "PATH directory name cannot be empty string" }, /*  */
		{ 440, "PATH directory name cannot end with /: %s" }, /*  */
		{ 441, "PATH directory name cannot contain a \\: %s\nUse / to separate directory names." }, /*  */
		{ 445, "Solution count file name already defined: %s" }, /*  */
		{ 446, "Lsd stop file name already defined: %s" }, /*  */
		{ 500, "Internal error.\n"
"Cannot create string buffer." }, /*  */
		{ 501, "Internal error.\n"
"Cannot create text buffer." }, /*  */
		{ 502, "Internal error.\n"
"Cannot store native substructure." }, /*  */
		{ 503, "Internal error.\n"
"Cannot bind %d with %d." }, /*  */
		{ 504, "Internal error.\n"
"Cannot bind %d with %d again." }, /*  */
		{ 505, "Internal error.\n"
"Atoms %d and %d should not correlate." }, /*  */
		{ 506, "Error during InChI creation: %s" } /*  */
	};


	if (v) {
		return ;
/* nothing to do */
	}
/* fprintf(stderr, "**************** error %d was detected.\n", errnum) ; */
	printf("error %d - %d commands read\n", errnum, lignes_lues) ;
	for (i = 0 ; i < MAX_MSG && msgs[i].errco != errnum ; i++) 
		;
	if (i != MAX_MSG) {
		va_start(l, v) ;
		vprintf(msgs[i].errmsg, l) ;
		va_end(l) ;
		printf("\n") ;
/* prints the message */
	} else {
		printf("unknown error!!!!\n") ;
/* should not happen */
	}
	check_lecture() ;
/*	histoire() ; */
	myexit(-1) ;
}


int	z(char *noyau)
/* 
 * looks up in the nucl table in order to get the nucleus reference 
 */
{
	int	i ;

	for (i = 0 ; (i < MAX_NUCL) && strcmp(noyau, nucl[i].sym) ; i++)  
		;
	non_nul(131, MAX_NUCL - i, noyau) ;
	return i ;
}


void	check_lecture(void)
/* lists for each read command the command number,
 * the mnemonic, the arguments as the v values determined in lirecom */
{
	int	i, j, o, n, fl ;

	for (i = 0 ; i < icarte ; i++) {
		o = cartes[i].opnum ;
		n = coms[o].fields ;
		fl = n < 0 ;
/* for the LIST command */
		if (fl) {
			n = -n ;
		}
		printf("%4d : %4s  ", i + 1, coms[o].op) ;
		for (j = 0 ; j < n ; j++) {
			printf("%4d", cartes[i].pars[j]) ;
		}
		if (fl) {
			check_list(cartes[i].pars[0]) ;
		}
		printf("\n") ;
	}
/* could be done in a better way, but would be more expensive */
}


void	check_list(int l)
/*
 * displays a (logical) list 
 */
{
	int	i ;

	for (i = 0 ; i < MAX_ATOMES ; i++) {
		if (list[l][i]) {
			printf("%3d", i) ;
		}
	}
}


