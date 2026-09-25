/*
    file: LsdSrc/lsd.c
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
#include <time.h>

#include <signal.h>
#include <stdlib.h>
/* new with 3.1.7, 9/16/2005, template of exit() */

extern void	lecture(void) ; /* lecture.c */
extern void	non_nul(int errnum, int v, ...) ; /* lecture.c */
extern void	start(void) ; /* start.c */
extern int	estlie(int x, int y) ; /* start.c */
extern void	about(void) ; /* start.c */
extern void	wrln2(int y) ; /* start1.c */
extern void	display(int p) ; /* display.c */
extern int	multiplelia(void) ; /* double.c */
extern int	antibredt(void) ; /* double.c */
extern void	native_sstr(void) ; /* sub.c */
extern void	readfragments(void) ; /* sub.c */
extern int	attribuer(void) ; /* sub.c */
extern void	initequ(void) ; /* equiv.c */
extern int	calcstat(struct atom *a) ; /* equiv.c */
extern int	testequ(int x) ; /* equiv.c */
extern void	modequ(int x) ; /* equiv.c */
extern void	marque(void) ; /* heap.c */
extern void	empile0(int *a, int v) ; /* heap.c */
extern void	empile(int t, int a, int v) ; /* heap.c */
extern void	empile6(int a, int v, int c) ; /* heap.c */
extern void	relache(void) ; /* heap.c */
extern void	histoire(void) ; /* heap.c */
extern void	remcorr3(int fl, int n, int x, int y) ; /*remcorr.c */
extern int	inutil1(int x, int y) ; /*remcorr.c */
extern int  remcorr2(int x, int z, int d) ; /*remcorr.c */
extern int	inutil2(int x, int z, int y) ; /*remcorr.c */
extern void	remcorrvar(int x, int z, int d) ; /*remcorr.c */
extern int	maxjok(void) ; /* distance.c */
extern int	setcnx(int x, int y, int fl) ; /* connexe.c */
extern int	isom(void) ; /* isom.c */
extern int	identical(void) ; /* isom.c */
extern void	isom_free(void) ; /* isom.c */

extern struct 	noy nucl[MAX_NUCL] ; /* lecture.c */
extern int	nssl ; /* sub.c */
extern int	nass ; /* sub.c */

void	ouvrir(int argc, char *argv[]) ;
void	dynvars(void) ;
void	stopit(void) ;
void	updatesolncounter(void) ;
void	myexit(int i) ;
void	controlC(int i) ;
void	timelimit() ;
void	work(void) ;
void	chlies(void) ;
int	insatur(int x) ;
void	queue(int x, int y, int n) ;
void	j1hyp(int x, int y) ;
void	j2hyp(int x, int y) ;
void	j2_1(int x, int y) ;
void	j2_2(int x, int y) ;
int	useprop(int x, int y) ;
int	usep(int x, int y) ;
int	valide(int x, int y) ;
int	goodsp2(int x, int flx) ;
int	goodsp(int x, int flx) ;
void	ajbase(int x) ;
void	wrbase(void) ;
void	paires(void) ;
void	pair1(void) ;
void	finir(void) ;
int	dupsol(void) ;
void	classe(int n, int t[]) ;
int	ok(void) ;
void	piege(void) ;

FILE *fic1, *fic2 ; /* input, output file */
int flclose = FALSE ; /* TRUE if fic1 and fic2 need to be closed */

struct atom	at[MAX_ATOMES] ; /* all the informations about atoms */
struct pile	heap[MAX_HEAP] ; /* heap for the optimization of recursivity */
struct lcorr propcorr[MAX_CORR] ; /* informations about HMBC and COSY correlations */
float reels[MAX_REELS] ; /* contains all real numbers */
int	list[MAX_LIST][MAX_ATOMES] ; /* logical lists */
int	valcorr[MAX_CORR] ; /* to know if correlations have been used */
int	max_atomes = 0 ; /* highest atom number */
int	icorr = 0 ; /* correlations counter, starts at 0 */
int	phase = 1 ; /* 1, 2 or 3, see work() */
int	nivheap = 0 ; /* heap level */
int	flpaires = FALSE ; /* for final atom pairing */
int	iessai ; /* for final atom pairing */
int	icomp ; /* for final atom pairing */
int	comp[MAX_ATOMES] ; /* for final atom pairing */
int	essai[MAX_BOND][2] ; /* for final atom pairing */
int	isol = 0 ; /* solution number */
int	iph2 = 0 ; /* number of entries into the pairing process */
int	flabout = FALSE ; /* to display input final interpretation */
int	flhistoire = FALSE ; /* to display the resolution path */
int	fldisp = 1 ; /* 0 : displays soutions as bond, 1: for JMNMOL */
int	flverb = FALSE ; /* verbosity, 0, 1, 2 */
int	flpart = FALSE ; /* gives partial solutions */
int	flbackph2 ; /* true upon coming back successfully from atoms pairing */
int	flstep = FALSE ; /* step by step operation */
int	flwork = TRUE ; /* true if the solution(s) has to be found */
int	fldupl = 2 ; /* 0: no check, 1: no duplicates, 2: no isomorphs */
int	flsubsgood = TRUE ; /* true if substructure must be present */
int	flsubsbad = FALSE ; /* true if substructure must be absent */
int	flcnx = TRUE ; /* true if the solution structure must be connected */
int flbredt = TRUE ; /* test for anti-Bredt structure */
int	base[MAX_ATOMES] ; /* determines the order of the use of correlations */
int	nbase = 0 ; /* number of atoms in the base */
int	lnsize[MAX_LN] ; /* sizes of numerical lists */
int	ln[MAX_LN][MAX_LNBUF] ; /* numerical lists */
int	flcorrvar = 0 ; /* hmbc and cosy variant correlation(s) counter */
int	supniv = 0 ; /* highest heap level */
int	nivdisp = 0 ; /* heap level at which full program stop is required */
int	elim = 0 ; /* number of HMBC/COSY correlations to eliminate */
int	flelim = FALSE ; /* HMBC/COSY correlations can be eliminated */
int	maxj = 0 ; /* max number of bonds between atoms (ELIM command) */
int	flfilt = 0 ; /* use LSD as substructure filter */
dynstring *curstring = NULL ; /* F type (filters, fragments) argument storage */
dyntext *strtab = NULL ; /* container of all strings and related information */
int	fldeff = 0 ; /* number of substructure(s) defined by DEFF and/or SKEL command(s) */
int	flfexp = 0 ; /* number of substructure(s) defined by the FEXP command */
int flmaxfrag = 0 ; /* maximum index of substructures defined by DEFF and FEXP commands */
int	flnative = 0 ; /* number of native substructure (0 or 1) */
int	flmaxstruct = 0 ; /* max number of structures on output, 0 is no limit */
int	flmaxtime = 0 ; /* max run time in seconds, 0 is no limit */
char	*reason = "Ctrl-C" ; /* reason for which execution was interrupted */
time_t	tinit ; /* initial time */
int	numskelatoms = 0 ; /* number of skeletal atoms */
int	numskelbonds = 0 ; /* number of skeletal bonds */
int	flCclasses = FALSE ; /* equiv. classes for C are forbidden/allowed 0/1 */
dyntext *pathtab = NULL ; /* container of solution count and lsd stop path names */
char *countfile = NULL ; /* name of solution count file name */
char *stopfile = NULL ; /* name of lsd stop file name */

int main(int argc, char *argv[])
/*
 * input, output file set-up and resolution top-level
 */
{
/* store time at which execution has started */
	tinit = time(NULL) ;
	
/* header */
	fprintf(stderr, "LSD-3.5.2\n") ;
	fprintf(stderr, "Copyright(C)2000 CNRS - UMR 7312 - Jean-Marc Nuzillard.\n") ;
	fprintf(stderr, "Please read the LICENCE file.\n") ;
	fprintf(stderr, "jm.nuzillard@univ-reims.fr.\n") ;

	ouvrir(argc, argv) ;    /* defines input, output files */
	dynvars() ; 	/* initializes curstring and strtab */
	lecture() ;/* reads input file */
	start() ;  /* input file analysis */
	stopit() ; /* stops everything is the stop file is there */
	updatesolncounter() ; /* creates the file that contains the number of solutions */
	native_sstr() ;
/* immediately after calling start(), nssl is the true number of subbonds.
 * after native_sstr it contains the number of recursive steps in the
 * assignment process for the native substructure. The latter is the one
 * that is defined in the lsd input file, by opposition to the ones that
 * defined by DEFF commands. */
  	readfragments() ;
 	heap[0].niv = 0 ;  /* heap initialisation */
	if (flfilt) {
		myexit(attribuer()) ;
/* new with 3.1.2. LSD for subs-structure filtering only */
	}
	initequ() ;/* equivalence classes initialisation */
	if (flabout) {
		about() ;
	}
/* displays the starting state of the resolution process */
	if (flwork) {
		work() ;
/* solution(s) search */
		if ((fldisp == 1) && (isol)) {
			fprintf(fic2, "0\n") ;
/* means EOF for OUTLSD */
		}
	}
	if (flwork && (isol == 0) && (nivdisp == 0)) {
/* no solution has been found. Prints the highest reached heap level (nivdisp)
 * rerunning lsd with the "MLEV nivdisp" command will dump the resolution
 * state at its most advanced point */
		fprintf(stderr, "No solution found.\n"
"You should try to:\n"
"  - Check all you MULT commands\n"
"  - Remove commands for small intensity COSY and HMBC correlations\n"
"  - Let LSD remove HMBC and COSY correlations using the ELIM command\n"
"  - Remove BOND commands\n"
"  - Remove sub-structure constraints\n"
"Max stack level : %d\n", supniv) ;
		myexit(0) ;
	}
	if (isol) {
		fprintf(stderr, "%d solution%s found.\n",
			isol, (isol>1) ? "s" : "") ;
	}
        myexit(1) ;
} /* main */


void	ouvrir(int argc, char *argv[])
/*
 * ouvrir means "to open". sets FILE* variables according to lsd call line
 */
{
	char	*nom1, *nom2 ;
	dynstring	*dnom2 ;
	int	len ;
	char	*ext = ".lsd" ;
	int	c ;
	int	i, j ;


	switch (argc) {
	case 1 :
/* lsd without argument : uses standard devices */
		fic1 = stdin ;
		fic2 = stdout ;
/* data are not copied to output file */
		break ;
	case 2 :
/* lsd with one argument. the file "file" is opened for reading,
 * the file "file.sol" is opened for writing */
 		nom1 = argv[1] ;
		fic1 = fopen(nom1, "r");
		if (fic1 == NULL) {
			printf("Cannot open file: %s\n", nom1) ;
			exit(-1) ;
		}
/* if the file name ends with .lsd, this extension must be removed first */
		len = strlen(nom1) ;
		i = len - 1 ;
		j = strlen(ext) - 1 ;
		for ( ; (i >= 0) && (j >= 0) && (nom1[i] == ext[j]) ; i--, j-- )
			;
		dnom2 = dynstring_new(DYNSTRING_BLOCKSIZE) ;
		non_nul(500, dnom2 != NULL) ;
		dynstring_pushsn(dnom2, nom1, ((j < 0) && (i >= 0)) ? i+1 : len) ;
		dynstring_pushs(dnom2, ".sol") ;
		nom2 = dynstring_getstring(dnom2) ;
		fic2 = fopen(nom2, "w");
		if (fic2 == NULL) {
			printf("Cannot create file: %s\n", nom2) ;
			dynstring_free(&dnom2) ;
			fclose(fic1) ;
			exit(-1) ;
		}
		dynstring_free(&dnom2) ;
		flclose = TRUE ;
		fprintf(fic2, "# From file: %s.\n", nom1) ;
/* header, line 1, new with 3.1.2 */
		while ( (c = fgetc(fic1)) != EOF ) {
			fputc(c, fic2) ;
		}
/* header, data are copied to output file */
		fputs("#\n", fic2) ;
/* header, last line */
		rewind(fic1) ;
/* input file ready for analysis */
		break ;
	default :
		printf("Usage : lsd [input file name]\n") ;
		exit(-1) ;
		break ;
	}
	signal(SIGINT, controlC) ;
} /* ouvrir */


void	dynvars(void)
/*
 * sets curstring and dyntext
 */
{
	curstring = dynstring_new(DYNSTRING_BLOCKSIZE) ;
	non_nul(500, curstring != NULL) ;
	strtab = dyntext_new(DYNTEXT_BLOCKSIZE) ;
	non_nul(501, strtab != NULL) ;
	pathtab = dyntext_new(DYNTEXT_BLOCKSIZE) ;
	non_nul(501, pathtab != NULL) ;
}


void stopit()
/*
 * stops execution if lsd stop file exists. 
 */
{
	static int first = TRUE ;
	
	if (fopen(stopfile, "r") != NULL) {
		if (first) {
			fprintf(stderr, "Please remove the %s file.\n", stopfile) ;
			exit(0) ;
		} else {
			reason = "STOP file detected" ;
			controlC(0) ;
		}
	}
	first = FALSE ;
}


void	updatesolncounter(void)
/*
 * writes the number of found solutions to the solution count file 
 */
{
	FILE *f ;

	if ((f = fopen(countfile, "w")) != NULL) {
		fprintf(f, "%d\n", isol) ;
		fclose(f) ;
	}
}


void	myexit(int i)
/*
 * exits from LSD with appropriate return code
 * closes input and output files if needed
 * frees memory that is occupied by dynamic objects
 */
{
	fflush(fic2) ;
	if (flclose == TRUE) {
		if(fic1 != NULL) fclose(fic1) ;
		fclose(fic2) ;
	}
	dynstring_free(&curstring) ;
	dyntext_free(&strtab) ;
	isom_free() ;
	updatesolncounter() ;
	dyntext_free(&pathtab) ;
	exit(i) ;
}

void	controlC(int i)
/*
 * cleanly closes files when the user types Ctrl-C
 * or when the maximum number of structures is reached
 * or when run time is too long.
 */
{
	fprintf(stderr, "Aborted by user (%s) after %d solution%s.\n",
		reason, isol, (isol>1) ? "s" : "") ;
	fprintf(fic2, "0\n") ;
	myexit(0) ;
}

void	timelimit()
/*
 * terminates execution if run time is too long.
 */
{
 	time_t tcurrent ;
 	
 	tcurrent = time(NULL) ;
 	if(difftime(tcurrent, tinit) >= (double)flmaxtime) {
 		reason = "time limit" ;
		controlC(0) ;
	}
}

void	work(void)
/*
 * Warning : serious stuff begins here !
 * the variable phase controls the resolution process. It has been initialised
 * to 1. Work is called recursively from chlies() (means : make bonds from
 * C-H correlations) and from paires()
 */
{
 	static long n = 0 ;

	n++ ;
	if(n == CALL_TIME_EVERY) {
		n = 0 ;
/* reset call count */
		updatesolncounter() ;
/* reports the number of found solutions, so far */
		if(flmaxtime) {
			timelimit() ;
/* handling time limit, new with 3.2.4 */
		}
		stopit() ;
/* looks for the presence of the file that stops execution */
	}
	switch (phase) {
/* algorithm execution according to current resolution phase */
	case 1 : 
		if (flverb > 1) {
			printf("%d : entering phase 1\n", nivheap) ;
		}
		chlies() ;
		if (flverb > 1) {
			printf("%d : leaving phase 1\n", nivheap) ;
		}
/* uses correlations */
		break ;
	case 2 : 
		if (flverb > 1) {
			printf("%d : entering phase 2\n", nivheap) ;
		}
		paires() ;
		if (flverb > 1) {
			printf("%d : leaving phase 2\n", nivheap) ;
		}
/* uncomplete atoms pairing */
		break ;
	case 3 : 
		if (flverb > 1) {
			printf("%d : entering phase 3\n", nivheap) ;
		}
		finir() ;
		if (flverb > 1) {
			printf("%d : leaving phase 3\n", nivheap) ;
		}
/* validity checks and printout */
		break ;
	}
} /* work */

void	chlies(void)
/*
 * phase 1 : uses correlations
 */
{
	int	go, i, x, y, n ;

	go = TRUE ;
/* will be false when a suitable correlating atom x will be found */
	if (flverb > 1) {
		wrbase() ;
	}
/* the base is an initially empty set that is progressively filled  
 * with the atoms most recently involved in bonds. If the base is not empty
 * correlations of atoms present in the base are considered in priority.
 * When all the correlations of an atom from the base are used this
 * atom is removed from the base. */
	if (nbase) {
/* the base is not empty */
		for (i = 0 ; (i < 5) && go ; i++) {
			for (x = 1 ; (x <= max_atomes) && (go = !( base[x] && (insatur(x) == i))) ; x++) 
				;
/* atoms from the base have all at least one usable correlation.
 * x has i free positions */
		}
	} else {
/* the base is empty */
		for (i = 0 ; (i < 5) && go ; i++) {
			for (x = 1 ; (x <= max_atomes) && (go = !( at[x].utile && (insatur(x) == i) && (at[x].nc != 0))) ; x++)
				;
/* x has a usable correlation and i free positions */
		}
	}
	if (i == 5) {
/* x has not been found. All the correlation data have been exploited */
		iph2++;
/* counts the number of times the program has entered into the phase #2
 * this has been useful for debugging */
		marque() ;
		empile0(&phase, 2) ;
/* now phase equals 2 */
		work() ;
/* recursive call */
		relache() ;
/* back to the original state */
	} else {	
/* x is found. now a correlation of x with y has to be found */
		for (i = 0 ; (n = at[x].ex[i], !valcorr[n]) ; i++) 
			;	
/* the first one which is valid will be the good one. This is x's ith 
 * correlation which is also the nth correlation */
		y = at[x].other[i] ;	
/* x correlates with y */
		if (flverb) {	/* message */
			printf("%d : analysis of correlation %d - ",  nivheap, x) ;
			wrln2(y) ;
		}
		marque() ;
		remcorr3(FALSE, n, x, y) ;
/* correlation number n between x and y is removed. the first argument (FALSE)
 * indicates that this removal is due to the choice of a correlation */
		empile(1, x, y) ;
/* the choice of this correlation will be mentionned in the history listing */
		if(propcorr[n].j1ok || propcorr[n].j2ok) {
			queue(x, y, n) ;
/* analyses the possible consequences of the x-y correlation */
		}
/* wrong in 3.1.0 */
		if(elim && propcorr[n].jnok) {
/* wrong in 3.1.0 */
/* corrected 24/12/2003 in v3.1.1 */
/* the HMBC correlation can now be rejected */
/* COSY correlations can be rejected too in v3.4.1 */
			if (flverb) {	/* message */
				printf("%d : drop correlation %d - ",  nivheap, x) ;
				wrln2(y) ;
			}
			empile6(x, y, n) ;
/* trace of correlation elimination in the history listing */
/*			elim-- ; */
			empile0(&elim, elim - 1) ;
/* one less correlation to be possibly eliminated */
			work() ;
/* does as if the correlation does not exist */
/*			elim++ ; */
/* restores the number of eliminable correlations */
		}
		relache() ;
	}
}


int	insatur(int x)	
/*
 * returns x's number of free positions
 */
{
	return at[x].nlmax - at[x].nlia ;
}


void	queue(int x, int y0, int n)
/* 
 * makes bond(s) from correlation of x with y0, indexed by n
 */
{
	int	i, y ;

	if (y0 > 0) {
		y = y0 ;
/* x correlates with an atom, not with a group of atoms */
		if(propcorr[n].j1ok) {
			j1hyp(x, y) ; 
/* x and y can directly be bound */
		}
		if(propcorr[n].j2ok) {
			j2hyp(x, y) ;
/* x and y can be bound to a common atom */
		}
	} else {
/* y0 stands for a group of atoms */
		for (i = 0 ; i < lnsize[-y0] ; i++) {
			y = ln[-y0][i] ;
/* y correlates possibly with x */
			/*if (testequ(y)) {*/
/* test considered useless since version 3.4.2 */
/* y has to the available according to equivalence status */
				marque() ;
				modequ(y) ;
/* updates y's equivalence class (if any) */
				empile(5, x, y) ;
/* the choice of y as a substitute for y0 is registered for the
 * history listing */
				if (flverb) {
					printf("%d : hypothesis : correlation %d - %d\n", nivheap, x, y) ;
				}
				if(propcorr[n].j1ok) {
					marque() ;
/* considering the j1 hypothesis, other j1ok correlation will be removed
and have to be restaured before considering the j2 hypothesis */
					remcorrvar(x, y, 1) ; 
/* if x correlates with one (or more) lists that contains y, that is (are) different from y0
and that partially overlap(s) with y0, then this (these) lists must be removed */
					remcorrvar(y, x, 1) ; 
/* if y correlates with one (or more) lists that contains x,
then this (these) lists must be removed */
					j1hyp(x, y) ;
					relache() ;
				}
				if(propcorr[n].j2ok && !estlie(x, y)) {
/* If x and y are bound the correlation cannot be j1ok because
 * it would have been removed before. The current y is not
 * a good choice */
					marque() ;
					remcorrvar(x, y, 2) ; 
					remcorrvar(y, x, 2) ; 
					j2hyp(x, y) ;
					relache() ;
				}
/* x and y are possibly either bound together or bound to a common atom. */
				relache() ;
			/*}*/
		}
	}
}


void	j1hyp(int x, int y)
/*
 * binds x and y, performs housekeeping
 */
{
	if (insatur(x) && insatur(y)) {
/* both x and y have to own more then 0 free position */
		marque() ;
		if (useprop(x, y) && valide(x, y) && inutil1(x, y) && inutil1(y, x)) {
/* the x-y bond has to be compatible with x's and y's property information.
 * housekeeping for x and y 
 */
			if (flverb > 1) {
			printf("%d : tries j1 %d - %d\n", nivheap, x, y) ;
			}
			if (flverb) {
				printf("%d : *** OK for j1 %d - %d\n", nivheap, x, y) ;
			}
			if (ok()) { 
/* eventual user's agreement */
				ajbase(x) ;
				ajbase(y) ;
/* base management */
				empile(2, x, y) ;
/*				inutil1(x, y) ; */
/*				inutil1(y, x) ; */
/* removes eventual confirmation correlations arising from the new x-y bond */
				work() ;
/* goes on */
			}
		}
		relache() ;
	}
}


void	j2hyp(int x, int y)
/*
 * in j2_1 and j2_2 the program will look for z to be bound to x and y.
 * if x is complete then z is obviously one of x's neighbours. This
 * situation is handled by j2_1(x, y). x and y roles are symmetrical.
 */
{
	if (insatur(x)) {
		if (insatur(y)) {
			j2_2(x, y) ;
		} else {
			j2_1(y, x) ;
		}
	} else {
		if (insatur(y)) {
			j2_1(x, y) ;
		} else {
			return ;
/*
 * this cannot happen when x and y are bound to a common atom because
 * x and y would be visible, and because correlations between visible
 * atoms are systematically removed before this point is reached. 
 * therefore x and y must have non-overlapping neighbours sets, and
 * this is in contradiction with the existence if a correlation between them.
 * The program must backtrack if this point is reached. Subtle isn't it?
 */
		}
	}
}


void	j2_1(int x, int y)
/*
 * finds the intermediate between x and y as one of x's neighbours when
 * x is already complete.
 */
{
	int	i, z ;

	for (i = 0 ; i < at[x].nlia ; i++) {
		z = at[x].lia[i] ;
		if (insatur(z)) {
/* z is x's ith neighbour, it can at least be bound to a new atom
 * and complies to equivalence constraints */ 
			marque() ;
			if (useprop(y, z) && valide(y, z) && inutil2(y, z, x)) {
/* the y-z bond is compatible with y and z properties. uptades data about
 * x and y */
				if (flverb > 1) {
					printf("%d : tries j2 %d - %d - %d\n", nivheap, x, z, y) ;
				}
				if (flverb) {
					printf("%d : *** OK for j2 %d - %d - %d\n", nivheap, x, z, y) ;
				}
				if (ok()) {
/* upon eventual user's agreement */
					ajbase(y) ;
					ajbase(z) ;
/* updates the content of the base */
					empile(2, y, z) ;
/* the y-z bond is stored for history listing */
/*					inutil2(y, z, x) ; */
/* removes eventual confirmation correlations. New bond is y-z. Bond x-z already exists */
					work() ;
/* goes on */
				}
			}
			relache() ;
		}
	}
}


void	j2_2(int x, int y)
/*
 * finds the intermediate z between x and y when neither x nor y is complete
 */
{
	int	z, flx, fly, m, valx, valy ;
	int val ;

	for (z = 1 ; z <= max_atomes ; z++) {
		if (at[z].utile && (z != x) && (z != y) && testequ(z)) {
/* z is different from x and from y and complies to equivelence contraints */
			flx = estlie(x, z) ;
/* if true there is no need to bind x and z */
			fly = estlie(y, z) ;
/* if true there is no need to bind y and z */
			non_nul(505, !(flx && fly), x, y) ;
/* x and y would be visible. it is impossible because the pending
 * correlation shoud have been removed previously */
 /*
			if (flx && fly) {
				printf("%3d and %3d should not correlate\n", x, y) ;
				histoire() ;
				about() ;
				myexit(-1) ;
			}
*/
			m = (flx || fly) ? 1 : 2 ;
/* number of bonds to provide to z */
			if (insatur(z) >= m) {
/* z has enough free positions */
				marque() ;
				valx = (flx) ? TRUE : (useprop(x, z) && valide(x, z)) ;
/* checks for the x-z bond */
				valy = valx &&  ((fly) ? TRUE : (useprop(y, z) && valide(y, z))) ;
/* if OK checks for z-y bond */
				val = valy && ((flx) ? TRUE : inutil2(x, z, y)) && ((fly) ? TRUE : inutil2(y, z, x)) ;
				if (val) {
/* everything is OK, proceeds to housekeeping */
					if (flverb > 1) {
						printf("%d : tries j2 %d - %d - %d\n", nivheap, x, z, y) ;
					}
					if (flverb) {
						printf("%d : *** OK for j2 %d - %d - %d\n", nivheap, x, z, y) ;
					}
					if (ok()) {
						if(m == 2) {
							modequ(z) ;
						}
						if (!flx) {
/* binds x to z */
							ajbase(x) ;
							ajbase(z) ;
							empile(2, x, z) ;
/*							inutil2(x, z, y) ; */
						}
						if (!fly) {
/* binds y to z */
							ajbase(y) ;
							ajbase(z) ;
							empile(2, y, z) ;
/*							inutil2(y, z, x) ; */
						}
						work() ;
/* goes on */
					}
				}
				relache() ;
			}
		}
	}
}


int	useprop(int x, int y)
/*
 * returns true if x and y can be neighbours according to their properties
 */
{
	return usep(x, y) && usep(y, x) ;
}


int	usep(int x, int y)
/*
 * return true if y can be bound to x according to x's properties 
 */
{
	int	i, n, ok, *a, inside ;
	struct atom *pa ;

	pa = &at[x] ;
	n = pa->nlpr ;
/* n is the number of x's property lists */
	if (!n) {
		return TRUE ;
/* no property */
	}
	
	ok = TRUE ;
/* becomes FALSE when a property cannot be satisfied */
	for(i = 0 ; (i < n) && ok ; i++) {
/* loop over properties of x and stop if one is not satisfied */
		inside = list[pa->lpr[i]][y] ;
/* true if y is inside of the ith property list of x */
		a = NULL ;
/* a adress of the number of remaining occurences, in occ if y is in the list, in extocc otherwise */
		if (pa->flocc[i] && inside) {
/* was "if (pa->flocc && inside) {" before 3.4.9 and was a bug */
			a = &pa->occ[i] ;
		}
		if (pa->flextocc[i] && !inside) {
/* was "if (pa->flextocc && inside) {" before 3.4.9 and was a bug */
			a = &pa->extocc[i] ;
		}
		if (a != NULL) {
/* the validity of the property must be checked */
			ok = *a ;
/* number of occurences, cannot be 0 because it has to be decremented */
			if (ok) {
/* change the number of remaining occurences */
				empile0(a, ok - 1) ;
			}
		}
	}
	return ok ;
}


int	valide(int x, int y)
/*
 * bond formation housekeeping
 */
{
	int	n, flx, fly, *a, nm ;

	non_nul(503, x != y, x, y) ;
	non_nul(504, !estlie(x, y), x, y) ;
	a = &at[x].nlia ;
/* address of the actual number of atoms bound to x */
	n = (*a) + 1 ;
/* what this number will be */
	nm = at[x].nlmax ;
/* its maximum value */
	if (n > nm) {
		return FALSE ;
/* of course */
	}
	flx = (n == nm) ;
/* x becomes complete */
	empile0(a, n) ;
/* updates the actual number of atoms bound to x */
	at[x].lia[n-1] = y ;
/* y is the last x's neighbour */

	a = &at[y].nlia ;
	n = (*a) + 1 ;
	nm = at[y].nlmax ;
	if (n > nm) 
		return FALSE ;
	fly = (n == nm) ;
	empile0(a, n) ;
	at[y].lia[n-1] = x ;
/* as well for y */

	return setcnx(x, y, TRUE) && goodsp2(x, flx) && goodsp2(y, fly) && goodsp(x, flx) && goodsp(y, fly) ;
/* binding x et y must satisfy the connectivity constraint, if any.
   a complete sp2 atom must have at least one sp2 neighbour */
}


int	goodsp2(int x, int flx)
/*
 * a complete sp2 atom must have at least one sp2 or sp neighbour
 */
{
	int	i, s ;

	if ((!flx) || (at[x].hybrid != 1)) {
		return TRUE ;
/* no point if either x is not complete or x is not sp2 */
	}
	s = 0 ;
	for (i = 0 ; i < at[x].nlmax ; i++) {
		s += at[at[x].lia[i]].hybrid ;
	}
/* s is the number of x's sp2 equivalent neighbours */
	return s ;
}


int goodsp(int x, int flx)
/*
 * a complete sp atom must have at least one sp or exactly two sp2 neighbours
 */
{
	int i, s ;
	 
	if ((!flx) || (at[x].hybrid != 2)) {
		return TRUE ;
/* no point if either x is not complete or x is not sp */
	}
	s = 0 ;
	for (i = 0 ; i < at[x].nlmax ; i++) {
		s += at[at[x].lia[i]].hybrid ;
	}
/* s is the number of x's sp2 equivalent neighbours */
	return (s >= 2) ;
}


void	ajbase(int x)
/*
 * attempts to add x into the base 
 */
{
	if (!base[x] && at[x].nc) {
/* x must not be already in the base and must have at least one valid 
 * correlation */
		empile0(&(base[x]), TRUE) ;
		empile0(&nbase, nbase + 1) ;
	}
}


void	wrbase(void)
/*
 * writes the base's content
 */
{
	int	x ;

	if (!nbase)  {
		printf("empty connected set\n") ;
	} else {
		printf("\t*%d : ", nbase);
		for (x = 1 ; x <= max_atomes ; x++) {
			if (base[x]) {
				printf("%d ", x) ;
			}
		}
		printf("\n") ;
	}
}


void	paires(void)
{
	int	x ;
	int	i ;
	int	ninsatur ;

	if (flpaires) {
/* when pairing initialisation has been performed */
		pair1() ;
	}
	else {
/* pairing initialisation */
		icomp = 0 ;
		for (x = 1 ; x <= max_atomes ; x++) {
			if (at[x].utile && insatur(x)) {
				comp[icomp++] = x ;
			}
		}
		if(flverb > 1) {
			printf("%d : Pairing initialisation\n", nivheap) ;
			for (i = 0 ; i < icomp ; i++) {
				x = comp[i] ;
				ninsatur = insatur(x) ;
				printf("atom %d : %d free position%s\n", x, ninsatur,
					(ninsatur > 1) ? "s" : "") ;
			}
		}
/* the table comp contains the icomp uncomplete atoms */
		marque() ;
		flbackph2 = FALSE ;
/* will become true only in partial mode and if pairing and 
 * phase #3 are successfull */
		empile0(&flpaires, TRUE) ;
		iessai = 0 ;
/* will be the number of times pair1 will be called recursively. it happens
 * each time a new bond is established, and the bound atoms will
 * be stored in essai[1][iessai] and essai[2][iessai] */
		work() ;
/* goes on */
		if (flbackph2 && !dupsol()) {
			isol++;
			display(fldisp) ;
/* in partial mode the solution is displayed here because only bonds
 * established during phase #1 are present */
		}
		relache() ;
	}
}


void	pair1(void)
{
	int	i, j, k, x, y, stoppe ;

	for (i = 0 ; (i < icomp) && (x = comp[i], !(insatur(x) && testequ(x))) ; i++)
		;
/* x the first available uncomplete atom from comp, at index i */
	if (i == icomp) {
/* no such x, goes to phase #3 */
		marque() ;
		empile0(&phase, 3) ;
		work() ;
		relache() ;
	} else {
/* ok for x */
		marque() ;
		modequ(x) ;
		for (j = i + 1 ; j < icomp ; j++) {
			y = comp[j] ;
/* y will be tentatively bound to x */
			if (x == y) {
				continue ;
/* how could it be possible ? */
			}
			if (!insatur(y)) {
				continue ;
/* y must still be uncomplete */
			}
			if (!testequ(y)) {
				continue ;
/* looks for equivalence */
			}
			if (estlie(x, y)) {
				continue ;
/* cannot bind x and y if it is already done */
			}
			if (flverb > 1) {
				printf("%d : tries to bind %d - %d\n", nivheap, x, y) ;
			}
			stoppe = FALSE ;
			for (k = iessai ; k > 0 &&  (essai[k][1] == x) &&  (stoppe = (essai[k][2] > y)) ; k--)
				;
/* goes back through the previously established bonds where x is the
 * starting atom. If any, the previous y must not be greater than the actual
 * one */
			if (stoppe) continue ;
			marque() ;
			if (useprop(x, y) && valide(x, y)) {
/* bonds housekeeping */
				if (flverb)
					printf("%d : *** OK to bind %d - %d\n", nivheap, x, y) ;
				if (ok()) {
					modequ(y) ;
/* looks for equivalence */
					empile0(&iessai, iessai + 1) ;
					essai[iessai][1] = x ;
					essai[iessai][2] = y ;
					empile(4, x, y) ;
/* stores x and y for the history listing */
					work() ;
/* attempts no establish one more bond */
				}
			}
			relache() ;
			if (flbackph2) {
/* in partial mode, if phase #3 has been successfully there is no need
 * to look for other pairing possibilities */
				break ;
			}
		}
		relache() ;
	}
}


void	finir(void)
/*
 * final tests before solution delivery
 */
{
	if(flverb > 1) {
		printf("Entering test for very long-range correlations\n") ;
	}
	if (flelim && (!maxjok())) {
		if(flverb > 1) {
			printf("solution rejected due to invalid long-range correlation(s)\n") ;
		}
		return ;
/* correlations could have been eliminated, distance control is set
 * distance constraint for atoms in eliminated correlations are checked */
	}
	if(flverb > 1) {
		printf("test for very long-range correlations passed\n") ;
	}
	if(flverb > 1) {
		printf("Entering test for double and/or triple bond placement\n") ;
	}
	if (!multiplelia()) {
		if(flverb > 1) {
			printf("solution rejected due to bad placement of double and/or triple bonds\n") ;
		}
		return ;
/* impossible to set double and/or triple bonds : backtrack */
	}
	if(flverb > 1) {
		printf("test for double and/or triple bond placement passed\n") ;
	}
	if(flbredt) {
		if(flverb > 1) {
			printf("Entering test for anti-Bredt solution\n") ;
		}
		if (antibredt()) {
			if(flverb > 1) {
				printf("anti-Bredt solution rejected\n") ;
			}
			return ;
/* antibredt solution : backtrack */
		}
		if(flverb > 1) {
			printf("test for anti-Bredt solution passed\n") ;
		}
	}
	if (!attribuer()) {
		if(flverb > 1) {
			if(flsubsgood) {
				printf("substructure expected but not found\n") ;
			} else {
				printf("substructure not expected but found\n") ;
			}
		}
		return ;
/* impossible to find the substructure : backtrack */
	}
	if(flverb > 1) {
		printf("test for sub-structure placement passed\n") ;
	}
	if (flpart) {
		flbackph2 = TRUE ;
/* the solution is a good one. in partial mode the test for unicity has
 * to be performed within the partial structures set */
		return ;
	} else if (dupsol()) {
/* only original solutions are printed */
		return ;
	}
	if(flverb > 1) {
		printf("test for uniqueness passed\n") ;
	}
	isol++;
	if(flverb > 1) {
		printf("Solution %d found\n", isol) ;
	}
	display(fldisp) ;
/* hurray !!!! */
	if (flhistoire) {
		histoire() ;
	}
	if(isol == flmaxstruct) {
/* handling the limit of structure number, new with 3.2.4 */
		reason = "too many solutions" ;
		controlC(0) ;
	}
}


int	dupsol(void)
/*
 * return true when the solution is duplicated.
 * with DUPL 1 duplicated means identical.
 * with DUPL 2 duplicated means isomorph.
 */
{
	if (fldupl == 0) {
		return FALSE ;
	}
/* nothing to do when DUPL 0 */
	return identical() || ((fldupl == 2) && isom()) ;
/* discard if solution is identical to a previous one.
discard isomorphs only if DUPL is 2 */
}


void	classe(int n, int t[])
/*
 * the most stupid sorting algorithm in the world
 */
{
	int	i, j, x ;

	for (i = 0 ; i < n - 1 ; i++)
		for (j = i + 1 ; j < n ; j++)
			if (t[i] > t[j]) {
				x = t[i] ;
				t[i] = t[j] ;
				t[j] = x ;
			}
}


int	ok(void)
/*
 * when flstep is true, ask the user to validate choices
 */
{
	if (!flstep)
		return TRUE ;
	printf("<CR>, n(o), f(inish), q(uit) ? ") ;
	switch (getchar()) {
	case 'n' :
		getchar() ;
		return FALSE ;
	case 'f' :
		getchar() ;
		flstep = FALSE ;
		return TRUE ;
	case 'q' :
		getchar() ;
		myexit(-1) ;
		break ;
	default :
		return TRUE ;
	}
}


void	piege(void)
/*
 * used for debugging with edge
 */
{
	printf("****************** coucou\n") ;
}


