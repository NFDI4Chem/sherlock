/*
    file: LsdSrc/sub.c
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

#include "defs.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef	WIN_NT
#include <windows.h>
#else
#include <sys/param.h>
#include <dirent.h>
#endif

extern int	flverb ; /* lsd.c */
extern int	flsubsgood ; /* lsd.c */
extern int	flsubsbad ; /* lsd.c */
extern struct atom 	at[MAX_ATOMES] ; /* lsd.c */
extern int	max_atomes ; /* lsd.c */
extern int	lnsize[MAX_LN] ; /* lsd.c */
extern int	ln[MAX_LN][MAX_LNBUF] ; /* lsd.c */
extern dynstring	*curstring ; /* lsd.c */
extern dyntext	*strtab ; /* lsd.c */
extern FILE	*fic1 ; /* lsd.c */
extern int	flclose ; /* lsd.c */
extern int	fldeff ; /* lsd.c */
extern int	flfexp ; /* lsd.c */
extern int	flnative ; /* lsd.c */
extern struct noy 	nucl[MAX_NUCL] ; /* lecture.c */
extern struct comtab 	coms[MAX_COM] ; /* lecture.c */
extern struct crt 	cartes[MAX_CARTES] ; /* lecture.c */
extern int	icarte ; /* lecture.c */
extern int	lignes_lues ; /* lecture.c */

extern void	marque(void) ; /* heap.c */
extern void	empile0(int *a, int v) ; /* heap.c */
extern void	empile(int t, int a, int v) ; /* heap.c */
extern void	relache(void) ; /* heap.c */
extern int	ok(void) ; /* lsd.c */
extern int	estlie(int x, int y) ; /* start.c */
extern int	estvalide(int x) ; /* start.c */
extern void	lecture(void) ; /* lecture.c */
extern void	non_nul(int errnum, int v, ...) ; /* lecture.c */
extern void	wrln2(int y) ; /* start.c */
extern int	stoss(int i) ; /* start.c */

struct ss_atom 	sub[MAX_ATOMES] ; /* sub atoms storage */
int	nssl = 0 ; /* number of bonds in the substructure - see main */
int	nass = 0 ; /* number of subatoms */
int	s1[MAX_NSSL] ; /* ordered subbonds starting point */
int	s2[MAX_NSSL] ; /* ordered subbonds ending point */
int	ss_ok ; /* true if the substructure has been found in the structure */

int numpath = 0 ; /* number of PATH command */

void	native_sstr(void) ;
int 	isdir(char *dirfilename) ;
void 	exploredir(char *dirIn, char *target, char *(*ppwhere)[2]) ;
int 	nam2file(int noreset, char *initDir, char *skel, char *(*ppwhere)[2]) ;
void	test_open_path(int p) ;
void	test_open_deff(int f, int p) ;
void	test_open_skel(int f, int p) ;
void	test_open_fexp(int p) ;
void	readfragments(void) ;
int	*readfragment(int f, char *filename) ;
void	sub_start(void) ;
void 	sub_start1(int ic) ;
void	unused_fragments(int n) ;
int	init_sstr(int f) ;
int	print_sub_index(int f) ;
int	*serialize() ;
int	critere(int x, int ssx) ;
int	crit_e(int x, int ssx) ;
int	crit_h(int x, int ssx) ;
int	crit_m(int x, int ssx) ;
int	attribuer(void) ;
int	fragment(int f) ;
void	explode(int f) ;
void	reset_atom_asgn(void) ;
void	attrib(int i) ;
void	criter2(int i, int y, int ssy) ;
void	sub_about(void) ;
void	sub_about2(int f, char *s) ;


void	native_sstr(void)
/*
 * stores native substructure information so that it will be possible
 * to refer to it as F0 (fragment 0) in a complexe substructure search.
 * The native substructure is the one that is not defined by a
 * DEFF command.
 */
{
	int n ;
	int *data ;

	nssl = init_sstr(0) ;
	flnative = (nass > 0) ? 1 : 0 ;
	dynstring_reset(curstring) ;
	n = dynstring_pushs(curstring, "NATIVE") ;
	non_nul(502, n) ;
	data = serialize() ;
	non_nul(502, data != NULL) ;
	n = dyntext_push(strtab, curstring, 0, data) ;
	non_nul(502, n >= 0) ;
}


int isdir(char *dirfilename)
/* returns -1 if file not accessible, 0 if not a directory, 1 if a directory */
{
#ifdef	WIN_NT
	DWORD FileAttributes = GetFileAttributes(dirfilename);

	if (FileAttributes == INVALID_FILE_ATTRIBUTES)
		return -1;

	if (FileAttributes & FILE_ATTRIBUTE_DIRECTORY)
		return 1;

	return 0;
#else
	struct stat dirfilebuf ;
	
	return (stat(dirfilename, &dirfilebuf) < 0) ? -1 : S_ISDIR(dirfilebuf.st_mode) ;
#endif
}


#ifdef	WIN_NT
#ifndef	MAXPATHLEN
#define	MAXPATHLEN MAX_PATH
#endif
#endif

void exploredir(char *dirIn, char *target, char *(*ppwhere)[2])
/* recursively searches in directories, starting for dirIn for a toc file that contains a ligne equal to target.
The path to the desidered structure file is in (*ppwhere)[0], if found. If another target is found its name is stored
in (*ppwhere)[1]. */
{
	char *where ;
	char *whereAlt ;
	int isDir ;
	char *stName ;
	char szFullName[MAXPATHLEN*2];
	int len ;
	FILE *toc ;
	int hastoc ;
	char skelName[MAXPATHLEN*2] ;
	char *whereFound ;
	int iline ;
	int quit ;
#ifdef	WIN_NT
	char szDirCopy[MAX_PATH];
	WIN32_FIND_DATA wfd;
	HANDLE hFindFile;
#else
	DIR *stDirIn ;
	struct dirent *stFiles ;
#endif
	
	len = strlen(dirIn) ;
	non_nul(431, len <= MAXPATHLEN, dirIn, len, MAXPATHLEN) ;
	where = (*ppwhere)[0] ;
	whereAlt = (*ppwhere)[1] ;
	isDir = isdir(dirIn) ;
/* fprintf(stderr, "%s is %sa directory\n", dirIn, isDir>0 ? "" : "not ") ; */
	non_nul(430, isDir == 1, dirIn, isDir) ;
#ifdef	WIN_NT
	strcpy(szDirCopy, dirIn);
	strcpy(szDirCopy + len, "/*");
	hFindFile = FindFirstFile(szDirCopy, &wfd);
	non_nul(432, hFindFile != INVALID_HANDLE_VALUE, dirIn) ;
#else
	stDirIn = opendir(dirIn) ;
	non_nul(432, stDirIn != NULL, dirIn) ;
#endif
	hastoc = 0 ;
	quit = 0 ;
#ifdef	WIN_NT
	while (!hastoc) {
		hastoc = !strcmp((stName = wfd.cFileName), "toc") ;
#else
	while (	((stFiles = readdir(stDirIn)) != NULL) && (!hastoc)) {
		hastoc = !strcmp((stName = stFiles->d_name), "toc") ;
#endif
		if(hastoc) {
			sprintf(szFullName, "%s/%s", dirIn, stName) ;
			len = strlen(szFullName) ;
			non_nul(433, len <= MAXPATHLEN, szFullName, len, MAXPATHLEN) ;
/* fprintf(stderr, "Analyzing %s\n", szFullName) ; */
			toc = fopen(szFullName, "r") ;
			non_nul(434, toc != NULL, szFullName) ;
			iline = 0 ;
			while( (fscanf(toc, "%s\n", skelName) == 1) && (!quit) ) {
				iline++ ;
				if(!strcmp(skelName, target)) {
/*
					whereFound = (strlen(where) == 0) ? where :
						((strlen(whereAlt) == 0) ? whereAlt : NULL) ;
					if(whereFound != NULL) {
						sprintf(whereFound, "%s/file%d", dirIn, iline) ;
/*
/* fprintf(stderr, "found %s in %s at line %d.\n", target, szFullName, iline) ;
fprintf(stderr, "\tstored at address %0x\n", whereFound) ; */
/*
					}
*/
					whereFound = (strlen(where) == 0) ? where : whereAlt ;
					sprintf(whereFound, "%s/file%d", dirIn, iline) ;
					len = strlen(whereFound) ;
					non_nul(433, len <= MAXPATHLEN, whereFound, len, MAXPATHLEN) ;
					quit = (strlen(whereAlt) != 0) ;
/* printf("quit = %d.\n", quit) ; */
				}
			}
			fclose(toc) ;
		}
#ifdef	WIN_NT
		if (FindNextFile(hFindFile, &wfd) == 0)
			break;
#endif
	}
	if(!quit) {
#ifdef	WIN_NT
		FindClose(hFindFile);
		hFindFile = FindFirstFile(szDirCopy, &wfd);
		while (!quit) {
			stName = wfd.cFileName;
#else
		rewinddir(stDirIn) ;
		while ( ((stFiles = readdir(stDirIn)) != NULL) && (!quit) ) {
			stName = stFiles->d_name ;
#endif
			if(stName[0] != '.') {
				sprintf(szFullName, "%s/%s", dirIn, stName) ;
				len = strlen(szFullName) ;
				non_nul(431, len <= MAXPATHLEN, szFullName, len, MAXPATHLEN) ;
				if(isdir(szFullName)) {
					exploredir(szFullName, target, ppwhere) ;
					quit = (strlen(whereAlt) != 0) ;
				}
			}
#ifdef	WIN_NT
		if (FindNextFile(hFindFile, &wfd) == 0)
			break;
#endif
		}
	}
#ifdef	WIN_NT
	FindClose(hFindFile);
#else
	closedir(stDirIn) ;
#endif
}


int nam2file(int noreset, char *initDir, char *skel, char *(*ppwhere)[2])
/* recursively searches from initDir for a toc file that contains a line with target. 
The path to structure file(s), if any, are stored in ppwhere (see exploreredir(). 
Must be called first with noreset=0 and then with noreset=1 so that directories are sucessively explored.
Returns 0 is target is not found, 1 if found once, 2 if found at least twice. */
{
	static char where[MAXPATHLEN+1] ;
	static char whereAlt[MAXPATHLEN+1] ;
	
	if(noreset == 0) {
		strcpy(where, "") ; strcpy(whereAlt, "") ;
		(*ppwhere)[0] = where ; (*ppwhere)[1] = whereAlt ;
	}
	exploredir(initDir, skel, ppwhere) ;
	
	return (strlen(where) == 0) ? 0 : ((strlen(whereAlt) == 0) ? 1 : 2) ;
}


void	test_open_path(int p)
/* checks for name at line p in strtab is the one of a directory */
{
	char *newpath, *shortpath, *longpath, *oldpath ;
	int ok ;
	int i, newlen, oldlen, shortlen, embedded ;
	
	newpath = dyntext_getstr(strtab, p) ;
	newlen = strlen(newpath) ;
	non_nul(439, newlen) ;
	non_nul(440, newpath[newlen-1] != '/', newpath) ;
	non_nul(441, strchr(newpath, '\\') == NULL, newpath) ;
	numpath++ ;
	dyntext_setindex(strtab, p, -1-numpath) ;
/* The successive paths get the -2, -3, -4... index. index -1 is for the FEXP argument, 0 for the native substructure, 1, 2, 3 for the others */
	for(i=(-numpath) ; i<(-1) ; i++) {
		non_nul(437, strcmp(newpath, dyntext_str_from_index(strtab, i)), newpath) ;
	}
/* duplication test */
	for(i=(-numpath) ; i<(-1) ; i++) {
		oldpath = dyntext_str_from_index(strtab, i) ;
		oldlen = strlen(oldpath) ;
		if(oldlen != newlen) {
			shortlen = oldlen < newlen ? oldlen : newlen ;
			embedded = !strncmp(oldpath, newpath, shortlen) ;
			if(embedded) {
				shortpath = oldlen < newlen ? oldpath : newpath ;
				longpath = oldlen < newlen ? newpath : oldpath ;
				non_nul(438, longpath[shortlen] != '/', shortpath, longpath) ;
			}
		}	
	}
/* inclusion test */
	ok = isdir(newpath) ;
	non_nul(430, ok == 1, newpath) ;
/* test if path exists and corresponds to a directory */
}


void	test_open_deff(int f, int p)
/* tries to open file whose name is at line p in strtab
   that is relative to fragment Ff */
{
	char *filename ;
	FILE *test ;
	
	filename = dyntext_getstr(strtab, p) ;
	test = fopen(filename, "r") ;
	non_nul(401, test != NULL, filename, f) ;
	fclose(test) ;
}

void	test_open_skel(int f, int p)
/* translates the  fragment Ff name at line p in strtab to a fragment file name. check the latter for opening */
{
	char *fragname, *dirIn ;
	char *pwhere[2] ;
	int i, imax, j ;
	line **ppl, *pl ;
	int numpatheff ;
	
	fragname = dyntext_getstr(strtab, p) ;
/* string at line p */

/* the case of numpath = 0 must be separately treated */
	numpatheff = numpath ? numpath : 1 ;
	for(i=(-1-numpatheff), j=0, imax=0; i<(-1) && (imax<2); i++, j++) {
		dirIn = numpath ? dyntext_str_from_index(strtab, i) : DEFAULT_PATH ;
/* fprintf(stderr, "Looking for %s in %s noreset is %d.\n", fragname, dirIn, j) ; */
		imax = nam2file(j, dirIn, fragname, &pwhere) ;
	}
/* error handling is missing */
	for(i = 0 ; i<imax ; i++) {
		fprintf(stderr, "Skeleton %s is described in %s.\n", fragname, pwhere[i]) ;
	}
	switch (imax) {
		case 0 : non_nul(435, 0, f, fragname) ; break ;
		case 1 : break ;
		case 2 : non_nul(436, 0, f, fragname) ; break ;
	}
	ppl = strtab->lines ;
	pl = dyntext_getline(strtab, p) ;
/* fprintf(stderr, "line was %s\n", pl->str) ;*/
	line_free(&pl) ;
	pl = line_new(pwhere[0], f, NULL) ;
/* fprintf(stderr, "line in now %s\n", pl->str) ; */
	ppl[p] = pl ;
	test_open_deff(f, p) ;
}

void	test_open_fexp(int p)
/* formats the expression string, either from LSD input file
 * or from external file, initializes the expression compiler,
 * checks if the fragments are defined.
 * The logical expression or the >filename is at line p in strtab.
 */
{
	char *exprinfo, *expr, *filename ;
	FILE *in ;
	char c, *pc ;
	plvm logic ;
	int fragc, *fragv, i ;

	exprinfo = dyntext_getstr(strtab, p) ;
/* exprinfo may either directly be a logical expression
or the name of a file that contains it */
	if(exprinfo[0] == '>') {
/* the file name starts immediately after the > sign */
		filename = exprinfo + 1 ;
		in = fopen(filename, "r") ;
		non_nul(403, in != NULL, filename) ;
		dynstring_reset(curstring) ;
/* the actual expression will be stored in curstring */
		while(!feof(in)) {
			dynstring_pushc(curstring, fgetc(in)) ;
		}
		fclose(in) ;
		expr = dynstring_getstring(curstring) ;
/* expression in expr is the string field of curstring */
	} else {
		expr = exprinfo ;
/* expression in expr is directly the one in strtab at line p */ 
	}
	for(pc = expr ; (c = *pc) != 0 ; pc++) {
		if(isspace(c)) {
			*pc = ' ' ;
		}
	}
/* changes all spaces to blanks for the lvm compiler. 
This operation does not change the length
of the string and may therefore be carried out in-place. */
	
	if(expr != exprinfo) {
		dyntext_setstr(strtab, p, expr) ;
/* this is not strictly speaking useful but may be, one day, a many
expression LSD input file will be implemented, so it's better to keep
the expression from external file in strtab */
	}
	logic = lvm_new(expr) ;
	non_nul(404, logic != NULL) ;
	non_nul(405, !lvm_error(logic), expr) ;
	lvm_vars(logic, &fragc, &fragv) ;
	for(i=0 ; i<fragc ; i++) {
		non_nul(406, dyntext_index_exists(strtab, fragv[i]), fragv[i]) ;
	}
	flfexp = fragc ;
	dyntext_setserial(strtab, p, (int *)logic) ;
/* I know, this is not elegant */
	dyntext_setindex(strtab, p, -1) ;
}


void	readfragments(void)
/* processes non-native fragment files */
{
	int i, l, f ;
	char *filename ;

	if(flclose) {
		fclose(fic1) ;
		fic1 = NULL ;
	}
/* closes the main input file so that fragment files can be read.
Setting fic1 to NULL avoids re-closing this file in myexit(). */
	l = dyntext_getlength(strtab) ;
	for(i=0 ; i<l ; i++) {
		f = dyntext_getindex(strtab, i) ;
		if(f > 0) {
			filename = dyntext_getstr(strtab, i) ;
			dyntext_setserial(strtab, i, readfragment(f, filename)) ;
		}
	}
	unused_fragments((flfexp == 0) ? fldeff : fldeff + flnative - flfexp) ;
}


int	*readfragment(int f, char *filename)
/* returns serialized sub-structure from file whose name is filename */
{
	fprintf(stderr, "Reading fragment file: %s\n", filename) ;
	fic1 = fopen(filename, "r") ;
/* should not fail because it was already tested in test_open_deff */
	lecture() ;
/* reads commands from fragment file */
	fclose(fic1) ;
	fic1 = NULL ;
/* no need any more for this input file */
	sub_start() ;
/* decodes commands from fragment file */
	nssl = init_sstr(f) ;
/* builds fragment assignment strategy */
	return serialize() ;
}


void	sub_start(void)
/* 
 * decodes commands from fragment file. This is a simplified copy
 * of start() in start.c
 */
{
	int ic, ip ;
	char *o ;
	
	nass = 0 ;
	nssl = 0 ;
	reset_atom_asgn() ;
/* resets substructure information */
	for (ic = 0 ; ic < icarte ; ic++) {
		o = coms[cartes[ic].opnum].op ;
		non_nul(410, !(
			strcmp(o, "SSTR") && 
			strcmp(o, "LINK") && 
			strcmp(o, "ASGN")), o) ;
	}
/* checks for substructure-only commands */
	for (ip = 10 ; ip <= 12 ; ip++) {
/* analyzing commands following ascending priority order,
8 is priority for SSTR, 9 for ASGN, 10 for LINK */
		for (ic = 0 ; ic < icarte ; ic++) {
			if (coms[cartes[ic].opnum].prior == ip) {
/* printf("# %d priority %d command %s\n", ic+1, ip, coms[cartes[ic].opnum].op) ; */
				lignes_lues = cartes[ic].ln ;
				sub_start1(ic) ;
			}
		}
	}
}


void 	sub_start1(int ic)
{
	int o, p1, p2, p3, p4, *pp ;
	struct crt *pcrt ;
	int i, j ;
	
	pcrt = &(cartes[ic]) ;
	o = pcrt->opnum ;
	pp = pcrt->pars ;
	p1 = pp[0] ;
	p2 = pp[1] ;
	p3 = pp[2] ;
	p4 = pp[3] ;
	switch (o) {
	case 32 : /* SSTR : substructure S A V V */
		non_nul(350, p1 > 0, p1) ;
		non_nul(351, !stoss(p1), p1) ;
/* p1 is a new subatom reference */
		nass++;
		non_nul(352, nass < MAX_NASS, MAX_NASS-1) ;
		sub[nass].ss_ref = p1 ;
		sub[nass].ss_element = p2 ;
		sub[nass].ss_hybrid = p3 ;
		if (p3 > 0) 
			non_nul(353, (p3 >= 1) && (p3 <= 3), p1, p3) ;
		else 
			for (i = 0 ; i < lnsize[-p3] ; i++) {
				j = ln[-p3][i] ;
				non_nul(353, (j >= 1) && (j <= 3), p1, p3) ;
			}
		sub[nass].ss_mult = p4 ;
		if (p4 >= 0) 
			non_nul(354, (p4 >= 0) && (p4 <= 3), p1, p4) ;
		else 
			for (i = 0 ; i < lnsize[-p4] ; i++) {
				j = ln[-p4][i] ;
				non_nul(354, (j >= 0) && (j <= 3), p1, p4) ;
			}
		sub[nass].ss_nlia = 0 ;
		sub[nass].ss_asgn = 0 ;
		sub[nass].ss_init = 0 ;
/* defines a subatom */
		break ;
	case 33 : /* ASGN : assign S I */
		i = stoss(p1) ;
		non_nul(360, i, p1) ;
		non_nul(361, estvalide(p2), p1, p2) ;
		non_nul(362, !sub[i].ss_asgn, p1, p2, sub[i].ss_asgn) ;
		sub[i].ss_asgn = p2 ;
		at[p2].asgn = i ;
/* updates substructure assignment information */
		break ;
	case 34 : /* LINK : link S S */
		non_nul(370, nssl < MAX_NSSL) ;
		i = stoss(p1) ;
		j = stoss(p2) ;
		non_nul(371, i, p1, p2, p1) ;
		non_nul(371, j, p1, p2, p2) ;
		non_nul(372, sub[i].ss_nlia < MAX_NLIA, p1, p2, p1, MAX_NLIA) ;
		non_nul(372, sub[j].ss_nlia < MAX_NLIA, p1, p2, p1, MAX_NLIA) ;
		sub[i].ss_lia[sub[i].ss_nlia++] = j ;
		sub[j].ss_lia[sub[j].ss_nlia++] = i ;
		nssl++;
/* set bonds between subatoms - validations could be improved */
		break ;
	}
}

void	unused_fragments(int n)
{
	if(n > 1) {
		fprintf(stderr, "Warning: there are unused fragments\n") ;
	}
	if(n == 1) {
		fprintf(stderr, "Warning: there is an unused fragment\n") ;
	}
}


int	init_sstr(int f)
/*
 * defines the strategy for subatoms assignment. The subbonds are reordered
 * so that the substructure will be exploited in order to built a
 * connected group of subatoms (when possible). s1[i] gives the start
 * of the ith subbond and s2[i] gives its end. s2[i] will be assigned according
 * to the assignment of s1[i]. If s1[i] is zero there is no preceeding
 * assignment to rely on. s2[j] is then the current starting subatom.
 */
{
	int	i, j, k, x, y, ilss, e ;

	if (!nass) {
		return 0 ;
	}
/* no subatom : nothing to do */
	ilss = nssl ;
/* nssl is the number of subbonds */
	i = j = 0 ;
	while (ilss) {
/* as long as there are subbonds */
		if (i == j) {
/* cold start */
			for (x = 1 ; (x <= nass) && !( sub[x].ss_asgn && !sub[x].ss_init) ; x++)
				;
/* x is assigned but still not used */
			if (x > nass) {
/* no x available */
				for (x = 1 ; (x <= nass) && sub[x].ss_init ; x++)
					;
/* x is not used */
			}
			s1[i] = 0 ;
			s2[i] = x ;
			sub[x].ss_init = 1 ;
			i++;
/* x is now used. subbond 0-x stored */
		} else {
			x = s2[j] ;
			if (sub[x].ss_init == 1) {
/* x is an ending point that can be used as a new starting point */
				for (k = 0 ; k < sub[x].ss_nlia ; k++) {
					y = sub[x].ss_lia[k] ;
					e = sub[y].ss_init ;
/* y is the kth subneighbour of x. Its state is e */
					if (e != 2) {
/* y has never been used as starting point */
						s1[i] = x ;
						s2[i] = y ;
						i++;
						ilss--;
/* the bond x-y is now stored */
						if (!e) {
							sub[y].ss_init = 1 ;
/* if y has never been an ending point, it becomes an ending point */
						}
					}
				}
				sub[x].ss_init = 2 ;
/* x has been used as a starting point */
			}
			j++;
/* ready for the next starting point */
		}
	}
	if (flverb) {
/* listing of the assignment order */
		print_sub_index(f) ;
		printf("assignment order\nx : ") ;
		for (j = 0 ; j < i ; j++) {
			printf("%4d", sub[s1[j]].ss_ref) ;
		}
		printf("\n") ;
		printf("y : ") ;
		for (j = 0 ; j < i ; j++) {
			printf("%4d", sub[s2[j]].ss_ref) ;
		}
		printf("\n") ;
	}
	return i ;
/* sets the number of recursive steps in the assignment process */
}


int	print_sub_index(int f)
{
	if (f == 0) {
		printf("Native substructure: ") ;
	} else {
		printf("Substructure F%d: ", f) ;
	}
}


int	*serialize()
/*
 * builds a vector of integers that represents the current 
 * substructure (or fragment)
 */
{
	int sz, szs1s2, szsub ;
	int *data ;
	int *p ;

	szs1s2 = nssl * sizeof(int) ;
	szsub = nass * sizeof(struct ss_atom) ;
	sz = 2 * sizeof(int) + 2 * szs1s2 + szsub ;
/* for nass and nssl (2 ints), s1 and s2 (2*nssl ints) and nass ss_atoms */
	data = (int *)malloc(sz) ;
	if(data == NULL) return NULL ;
	p = data ;
	*p = nass ;
	p++ ;
	*p = nssl ;
	p++ ;
	if(nssl) {
		memcpy(p, s1, szs1s2) ;
		p += nssl ;
		memcpy(p, s2, szs1s2) ;
		p += nssl ;
	}
	if(nass) {
		memcpy(p, sub+1, szsub) ;
	}
	return data ;
}


int	critere(int x, int ssx)
/*
 * checks the compatibility between the atom x and the subatom ssx
 */
{
	return crit_e(x, ssx) &&  crit_h(x, ssx) &&  crit_m(x, ssx) ;
}


int	crit_e(int x, int ssx)
/* 
 * checks the identity of elements
 */
{
	int e ;

	e = sub[ssx].ss_element ;
	if (nucl[e].val[0] < 0) {
		return 1 ;
/* atom type A, with valency -1, stands for any atom type - new with 3.1.4 */
	}
	return at[x].element == e ;
}


int	crit_h(int x, int ssx)
/* 
 * checks the compatibility of the hybridization states 
 */
{
	int	i, h, s ;

	h = sub[ssx].ss_hybrid ;
	if (h >= 0)
		return h == 3 - at[x].hybrid ;
	else {
		for (i = 0 ; i < (s = lnsize[-h]) &&  ln[-h][i] != 3 - at[x].hybrid ; i++) 
			;
		return i != s ;
	}
}

int	crit_m(int x, int ssx)
/* 
 * checks the compatibility of the multiplicities 
 */
{
	int	i, m, s ;

	m = sub[ssx].ss_mult ;
	if (m >= 0)
		return m == at[x].mult ;
	else {
		for (i = 0 ; i < (s = lnsize[-m]) &&  ln[-m][i] != at[x].mult ; i++)
			;
		return i != s ;
	}
}

int	attribuer(void)
/*
 * top level of the recursive assignment process
 */
{
	plvm logic ;

	if ((!flsubsgood) && (!flsubsbad)) {
		ss_ok = TRUE ;
	} else {
		if (flfexp == 0) {
			if (flnative == 0) {
				ss_ok = TRUE ;
			} else {
				ss_ok = fragment(0) ;
			}
		} else {
			logic = (plvm)dyntext_serial_from_index(strtab, -1) ;
			ss_ok = lvm_run(logic, fragment) ;
		}
		if (flsubsbad) {
			ss_ok = !ss_ok ;
		}
	}
	return ss_ok ;
}

int	fragment(int f)
/*  returns true if fragment Ff is included inside of the current structure.
this function is of lvm_eval type (see lvm.h) */
{
	if (flverb) {
		print_sub_index(f) ;
		printf("search started.\n") ;
	}
	explode(f) ;
	ss_ok = FALSE ;
	attrib(0) ;
	if (flverb) {
		print_sub_index(f) ;
		printf("%s found.\n", ss_ok ? "" : "not ") ;
	}
	return ss_ok ;
}

void	explode(int f)
/*
 * reverses the action of serializes: restores substructure data
 * as if it alsways had been here. f is a substructure index,
 * like in Ff.
 */
{
	int *p, szs1s2, szsub ;
	int x, i ;

	reset_atom_asgn() ;
/* resets atom assignment */
	p = dyntext_serial_from_index(strtab, f) ;
	nass = *p ;
	p++ ;
	nssl = *p ;
	p++ ;
	szs1s2 = nssl * sizeof(int) ;
	szsub = nass * sizeof(struct ss_atom) ;
	if(nssl) {
		memcpy(s1, p, szs1s2) ;
		p += nssl ;
		memcpy(s2, p, szs1s2) ;
		p += nssl ;
	}
	if(nass) {
		memcpy(sub+1, p, szsub) ;
	}
/* restores s1, s2 and the subatoms */
	for(i = 1 ; i <= nass ; i++) {
		x = sub[i].ss_asgn ;
		if(x) {
			at[x].asgn = i ;
		}
	}
/* restores atom assignment(s) */
}


void	reset_atom_asgn(void)
/*
 * removes any previous assignment of an atom to a subatom
 */
{
	int x ;

	for(x = 1 ; x <= max_atomes ; x++) {
		if(at[x].utile) {
			at[x].asgn = 0 ;
		}
	}
}


void	attrib(int i)
/*
 * assignment of the substructure. nssl is now the number of assignment
 * steps (see main).
 */
{
	int	x, y, ssx, ssy, j ;

	if (i == nssl) {
		ss_ok = TRUE ;
		return ;
/* end condition of recursivity. if this point is reached the substructure
 * is fully assigned. That means it is present in the structure */
	}
	ssx = s1[i] ;
	ssy = s2[i] ;
	y = sub[ssy].ss_asgn ;
	if (ssx) {
		x = sub[ssx].ss_asgn ;
/* the starting subatom ssx is already assigned to x */
		if (y) {
/* if y is already assigned, too */
			if (estlie(x, y)) {
				attrib(i + 1) ;
/* and if x and y are bound, everything is ok. goes on. */
			} else {
				return ;
/* something is wrong somewhere. goes back. */
			}
		} else {
			for (j = 0 ; j < at[x].nlia ; j++) {
/* x <--> ssx. ssx is bound to ssy. therefore if y <--> ssy, x and y
 * have to be bound. y is x's jth neighbour. */
				y = at[x].lia[j] ;
				criter2(i, y, ssy) ;
/* checks if y and ssy have compatible status and call attrib at the next level.
   performs then the actual assignment job */
				if (ss_ok) {
					break ;
/* once all the assignments done there is no need to look for an other
 * possible assignments set. Assignment are never used in the solution output */ 
				}
			}
		}
	} else {
/* ssy is assigned from scratch */
		if (y) {
			attrib(i + 1) ;
/* if ssy is already assigned there is nothing to do */
		} else {
			for (y = 1 ; y <= max_atomes ; y++) {
				if (!at[y].utile) {
					continue ;
				}
				criter2(i, y, ssy) ;
/* checks if y and ssy have compatible status and call attrib at the next level.
   performs then the actual assignment job */
				if (ss_ok) {
					break ;
				}
			}
		}
	}
}


void	criter2(int i, int y, int ssy)
{
	if ((!at[y].asgn) && critere(y, ssy)) {
/* ssy has never been assigned, therefore y cannot be already assigned.
 * y and ssy must be compatible */
		if (flverb) 
			printf("assigns S%d to %d\n", sub[ssy].ss_ref, y) ;
		if (ok()) {
			marque () ;
			empile0(&(sub[ssy].ss_asgn), y) ;
			empile0(&(at[y].asgn), ssy) ;
			attrib(i + 1) ;
			relache() ;
/* updates assignement information, to be released after backtracking */
		}
	}
}


void	sub_about(void)
/*
 * prints everything about substructures
 */
{
	int i, l, f ;
	char *s ;

	l = dyntext_getlength(strtab) ;
	for(i=0 ; i<l ; i++) {
		f = dyntext_getindex(strtab, i) ;
		if(f >= 0) {
			s = dyntext_getstr(strtab, i) ;
			sub_about2(f, s) ;
		}
	}
}


void	sub_about2(int f, char *s)
/*
 * prints details about substructure (fragment) Ff from file s
 */
{
	int x, n, i ;

	if (strcmp(s, "NATIVE")) {
		printf("substructure from file: %s\n", s) ;
	} else {
		printf("native substructure\n") ;
	}
	explode(f) ;
	if(!nass) {
		printf("\tis empty\n") ;
		return ;
	}
	for (x = 1 ; x <= nass ; x++) {
		printf("\tsub-atom %d\n", sub[x].ss_ref) ;
		printf("\t\telement %s\n",  nucl[sub[x].ss_element].sym) ;
		printf("\t\tmultiplicity : ") ;
		wrln2(sub[x].ss_mult) ;
		printf("\t\thybridization : ") ;
		wrln2(sub[x].ss_hybrid) ;
		n = sub[x].ss_asgn ;
		if (n) {
			printf("\t\tassignment of atom : %d\n", n) ;
		}
		printf("\t\tbound to : ") ;
		for (i = 0 ; i < sub[x].ss_nlia ; i++) {
			printf("%d ", sub[sub[x].ss_lia[i]].ss_ref) ;
		}
		printf("\n") ;
	}
}

