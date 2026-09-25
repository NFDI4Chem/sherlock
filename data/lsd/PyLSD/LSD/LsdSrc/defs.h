/*
    file: LsdSrc/defs.h
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

#define MAX_ATOMES 100
#define MAX_LIST 20
#define MAX_NC 12
#define MAX_NLPR 8
#define MAX_COM 58
#define MAX_PRIO 16
#define MAX_CARTES 300
#define MAX_NUCL 17
#define MAX_HEAP 10000
#define MAX_CORR 200
#define MAX_BOND 100
#define MAX_NLIA 4
#define MAX_EQU 30
#define MAX_LN 20 
#define MAX_LNBUF 5
#define MAX_NSSL 60
#define MAX_NASS 60
#define MAX_REELS (3*MAX_ATOMES)
#define DYNSTRING_BLOCKSIZE 80
#define DYNTEXT_BLOCKSIZE 10
#define CALL_TIME_EVERY 10000
#define DEFAULT_PATH "Filters"
#define SOLUTION_COUNT_FILE "solncounter"
#define LSD_STOP_FILE "stoplsd"

#define FALSE 0
#define TRUE 1

struct atom {
/* characterises atoms */
	int	utile ;/* true if used */
	int	element ;/* element index */
	int	valence ;/* valency */
	int	mult ;/* number of bound hydrogen atoms, multiplicity */
	int	hybrid ;/* 0 (sp3), 1 (sp2) or 2 (sp) hybridisation state */
	int charge ;/* charge */
	int ics ;/* index of chemical shift */
	int ihcs[2] ;/* index of eventual chemical shift of hydrogen bound to the atom */
	int	equstat ;/* equivalence class */
	int	iequ ;/* index within its equivalence class */
	int	nlmax ;/* highest possible number of bonds */
	int	nlia ;/* current number of bonds */
	int	lia[MAX_NLIA] ;/* indexes of the neigbours */
	int	ordre[MAX_NLIA] ;/* bond order, 1 or 2 */
	int	corrgrp ;/* true if correlates with a group */
	int	ingrp ;/* true if is group member */
	int	nctot ;/* number of correlations */
	int	nc ;/* number of valid correlations */
	int	ex[MAX_NC] ;/* index of correlations */
	int	other[MAX_NC] ;/* correlating atom or group (if <0) */
	int	asgn ;/* assignment according to the substructure */
	int	nlpr ;/* number of property lists */
	int	lpr[MAX_NLPR] ;/* indexes of the property lists */
	int	occ[MAX_NLPR] ;/* number of occurences within a list */
	int	extocc[MAX_NLPR] ; /* number of occurences outside a list */
	int	flocc[MAX_NLPR] ;/* TRUE of FALSE, use the number of occurences within a list */
	int	flextocc[MAX_NLPR] ; /* TRUE of FALSE, use the number of occurences outside a list */
	int	icnx ; /* index of connexe part */
	int	iidx ; /* inchizer index */
};

struct comtab {
/* informations about commands */
	char	*op ; /* command mnemonic */
	int	fields ; /* number of arguments in the command */
	char	*forme ; /* type of the arguments */
	int	prior ; /* priority during commands analysis */
};

struct crt {
/* digested commands */
	int	done ; /* false only for hmbc and cosy commands with variants */
	int	ln ; /* command number in the input file */
	int	opnum ; /* command index */
	int	pars[5] ; /* arguments coded as integer */
};

struct noy {
/* nuclei table */
	char	*sym ; /* input atomic symbol */
	int	val[4] ; /* list of possible valences (in function of charge) */
	float	mass ; /* atomic mass */
	int	count ; /* how many of such nuclei in the molecule */
	char	*symout ; /* output atomic symbol */
};

struct pilemem {
/* informations for recursivity simulation */
	int	*adr ; /* the address of an integer */
	int	val ; /* its value */
};

struct pilehist {
/* informations for history listing */
	int	h1 ;
	int	h2 ; 
	int corr ;
};

union piledata {
/* informations to pushed onto the heap */
	struct pilemem pmem ;
	struct pilehist phist ; 
};


struct pile {
/* a complete heap item */
	int	niv ; /* bloc number */
	int	ptyp ; /* item type (0 for pilemem, 1-5 for pile hist) */
	union piledata pinfo ; /* what has to be stored */
};

struct co {
/* about one particular equivalence class */
	int	valco ; /* coded status of the atoms within the class */
	int	ico ; /* class population */
};

struct ss_atom {
/* characterises subatoms */
	int	ss_ref ; /* its number in the input file */
	int	ss_element ; /* element index */
	int	ss_mult ; /* multiplicity */
	int	ss_hybrid ; /* hybridisation state */
	int	ss_nlia ; /* number of subbonds */
	int	ss_lia[MAX_NLIA] ; /* index of the subneighbours */
	int	ss_init ; /* 0, 1, 2 during subbonds sorting */
	int	ss_asgn ; /* index of the atom it is assigned to */
};

struct lcorr {
/* informations about HMBC and COSY correlations */
	int at1 ; /* first correlating atom */
	int at2 ; /* second correlating atom */
	int inputat1 ; /* first correlating atom of input file */
	int inputat2 ; /* second correlating atom of input file */
	int inputath1 ; /* first hydrogen number in correlation of input file (in case of COSY) */
	int inputath2 ; /* second hydrogen number in correlation of input file (COSY and HMBC) */
	int inputmin ; /* minimum correlation path length of input file */
	int inputmax ; /* maximum correlation path length of input file */
	int min ; /* minimum correlation path length after treatment */
	int max ; /* maximum correlation path length after treatment */
	int j1ok ; /* indicates if 2J is possible */
	int j2ok ; /* indicates if 3J is possible */
	int jnok ; /* indicates if nJ are possible */
	int origin ; /* origin of the correlation (0 from input file, 1 modified during treatment, 2 created during treatment) */
	int type ; /* type of correlation (0 for HMBC, 1 for COSY) */
};

typedef unsigned char	ucar ;

#include "lvm.h"
#include "dyntext.h"
