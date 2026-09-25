/*
    file: Mol2ab/mol2ab.c
    
    Copyright (c) 2009, Jean-Marc Nuzillard et Bertrand Plainchont, UMR 7312, 
    University of Reims Champagne-Ardenne, CNRS.
    All rights reserved.

    Redistribution and use in source and binary forms, 
    with or without modification, are permitted provided that the following
    conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
    THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
    OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
    WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
    OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
    EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/	

/*
 * Mol2ab converts structures in MOL format into substructures for LSD.
 * Usage: mol2ab Directory pilot-file
 * The pilot file contains lines with two fields:
 * - the first field is the skeleton name that will be associated to the .mol file.
 * - the second field is a .mol file name that must exist in the current directory,
 * The Directory must exist in the current directory.
 * It will be filled with the initial .mol files (renamed file1.mol ... filennn.mol)
 * and their LSD equivalent files (file1 ... filennn) if nnn is the number of lines
 * in the pilot file.
 * A toc file is also created so that LSD can retrieve substructure files from skeleton names.
 *
 * Structure to substructure conversion follows the following rules:
 * - Any hydrogen atom from the structure may be abstracted
 * - Any single bond may be converted into a double bond
 * - Any double bond may be converted into a triple bond
 * - A triple bond remains a triple bond
 * - Atomic elements are preserved
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define LONGMAX 255
#define MAXVAL 4 /* valence maximum for A (any) atom */

void myexit(int v) ;
void mol2ab(char *dir, char *pilot) ;
void eol(char *molname) ;
void readmol(char *molname) ;
void analyze(char *molname) ;
int	z(char *nucname) ;
void writemol(char *dir, int f, char *molname, char *skelname) ;
void copymol(char *dir, int f, char *molname, char *skelname) ;

 
/*
 * Nucleus definition
 */
typedef struct t_nuc {
	char *sym ;				/* pointer to atom symbol */
	int	val ;				/* valency */
} nuc ;
	
/*
 * A short nucleus table
 */
nuc nucl[] = {
	{ "C", 4 },
	{ "O", 2 },
	{ "N", 3 },
	{ "S", 2 },
	{ "F", 1 },
	{ "Cl", 1 },
	{ "Br", 1 },
	{ "I", 1 },
	{ "Si", 4 },
	{ "P", 3 },
	{ "A", MAXVAL }
};
#define MAX_NUCL (sizeof(nucl)/sizeof(nuc))

/*
 * Atom definition
 */
typedef struct t_atom {
	char elt[3] ;			/* element symbol */
	int valency ;			/* valency */
	int hybrid ;			/* hybridization state: 0 for sp3, 1 for sp2, 2 for sp atom */
	int totalbond ;			/* number of neighbors */
	int maxh ;				/* number of H atoms */
	char *multtext ;		/* multiplicity part of the LSD SSTR command */
	char *hybridtext ;		/* hybridization part of the LSD SSTR command */
	int ndb ;				/* number of double bond */
} atom ;

/*
 * Bond defintion
 */
typedef struct t_bond {
	int at1 ;				/* index of atom 1 */
	int at2 ;				/* index of atom 2 */
	int bondtype ;			/* 1 for a single bond, 2 for a double bond, 3 for a triple bond */
} bond ;

/* Global variables, related to the currently proceded substructure */
FILE *mol ;
int natom, nbond ;
atom *tabatom = NULL ;
bond *tabbond = NULL ;


int main(int argc, char *argv[])
{
	if(argc != 3) {
		fprintf(stderr, "%s: usage: mol2ab Directory pilot-file\n", argv[0]) ;
		exit(1) ;
	}	
	mol2ab(argv[1], argv[2]) ;
	myexit(0) ;
}

/*
 * exit without memory leak
 */
void myexit(int v)
{
	if(tabatom != NULL) free(tabatom) ;
	if(tabbond != NULL) free(tabbond) ;
	exit(v) ;
}

/*
 * mol2ab creates files in dir according to pilot
 */
void mol2ab(char *dir, char *pilot)
{
	FILE *toc ;
	FILE *fpilot ;
	int nmol ;
	char molname[LONGMAX] ;
	char skelname[LONGMAX] ;
	char buf[LONGMAX] ;
	char dirtoc[LONGMAX] ;
	
	sprintf(dirtoc, "%s/toc", dir) ;
/* dirtoc is the name of the toc file */
	if((toc = fopen(dirtoc, "w")) == NULL) {
		fprintf(stderr, "Cannot create file %s. Check for the existence of directory %s.\n", dirtoc, dir) ;
		myexit(1) ;
	}
/* ready to write in toc */
	if((fpilot = fopen(pilot, "r")) == NULL) {
		fprintf(stderr, "Cannot open file %s for reading.\n", pilot) ;
		fclose(toc) ;
		myexit(1) ;
	}
/* ready to read the pilot file */
	nmol = 0 ;
	while(!feof(fpilot)) {
		nmol++ ;
/* line count in the pilot file */
		if(fscanf(fpilot, "%s%s\n", skelname, molname) != 2) {
			fprintf(stderr, "Cannot read line %d from %s.\n", nmol, pilot) ;
			fclose(toc) ;
			fclose(fpilot) ;
			myexit(1) ;
		}
#ifdef DEBUG
printf("Reading %s\n", molname) ;
#endif
		if((mol = fopen(molname, "r")) == NULL) {
			fprintf(stderr, "Cannot open file %s for reading. Skipped.\n", molname) ;
			continue ;
		}
		readmol(molname) ;
/* read the current molecule */
		rewind(mol) ;
/* rewind the file so that it will be possible to copy it in dir */
		copymol(dir, nmol, molname, skelname) ;
/* copy the molecule file */
		fclose(mol) ;
/* finished with the current molecule */
		analyze(molname) ;
/* process the moecule so that it can become an LSD substructure */
		writemol(dir, nmol, molname, skelname) ;
/* write the LSD substructure file */
		fprintf(toc, "%s\n", skelname) ;
/* append the skeleton name to toc so that LSD can find the file from the skeleton name */
		free(tabbond) ; tabatom = NULL ;
		free(tabatom) ;	tabbond = NULL ;
	}
/* free the memory that was reserved for the molecule */
	fclose(toc) ;
	fclose(fpilot) ;
/* close pilot and toc files */
}

/*
 * moves the mol file pointer to the next end of line
 */
void eol(char *molname)
{
	static char buf[LONGMAX] ;
	
	if(fgets(buf, LONGMAX, mol) == NULL) {
		fprintf(stderr, "Error while reading %s.\n", molname) ;
		fclose(mol) ;
		myexit(1) ;
	}
}

/*
 * readmol reads a molecule from file pointer mol
 */
void readmol(char *molname)
{
	char buf[LONGMAX] ;
	atom *pcuratom ;
	bond *pcurbond ;
	int i ;
	
	eol(molname) ; eol(molname) ; eol(molname) ;
/* skip the 3 first lines in the current .mol file */
	if(fscanf(mol, "%d%d", &natom, &nbond) != 2)
	{
		fprintf(stderr,"Error in %s while reading the number of atoms/bonds.\n", molname) ;
		fclose(mol) ;
		exit(1) ;
	}
	eol(molname) ;
/* the 4th line contains the number of atoms and of bonds */
	
#ifdef DEBUG
	printf("%d atoms and %d bonds.\n", natom, nbond) ;
#endif
	
	tabatom = (atom *)malloc(natom * sizeof(atom)) ;
	tabbond = (bond *)malloc(nbond * sizeof(bond)) ;
	if(tabatom == NULL || tabbond == NULL) {
		fprintf(stderr, "Cannot allocate memory for molecule in %s.\n", molname) ;
		fclose(mol) ;
		exit(1) ;
	}
/* memory allocation for the current bond */
	
	for(i = 0, pcuratom = tabatom ; i < natom ; i++, pcuratom++) {
/* read the atom block of the .mol file */
		if(fscanf(mol, "%*s%*s%*s%s", pcuratom->elt) != 1) {
			fprintf(stderr, "Cannot read symbol of atom %d in %s.\n", i+1, molname) ;
			fclose(mol) ;
			exit(1) ;
		}
		eol(molname) ;
/* the element symbol is the only useful information in an atom line of a .mol file */
	}
	
#ifdef DEBUG
	printf("Atom table:\n") ;
	for(i = 0 ; i < natom ; i++) {
		printf("%s\n", tabatom[i].elt) ;
	}
#endif
	
	for(i = 0, pcurbond = tabbond ; i < nbond ; i++, pcurbond++) {
/* read the bond block of the current .mol file */
		if(fscanf(mol, "%d%d%d", &pcurbond->at1, &pcurbond->at2, &pcurbond->bondtype) != 3) {
			fprintf(stderr, "Cannot read bond %d in %s.\n", i+1, molname) ;
			fclose(mol) ;
			exit(1) ;
		}
		eol(molname) ;
/* get bonded atom indexes and bond type */
	}
	
#ifdef DEBUG
	printf("Bond table:\n") ;
	for(i = 0 ; i < nbond ; i++) {
		printf("%d  %d  %d\n", tabbond[i].at1, tabbond[i].at2, tabbond[i].bondtype) ;
	}
#endif

}

/*
 * analyze transforms the molecule into a substructure
 */
void analyze(char *molname)
{
	int v, i;
	atom *pa1, *pa2 ;
	atom *pcuratom ;
	bond *pcurbond ;

	for(i = 0, pcuratom = tabatom ; i < natom ; i++, pcuratom++) {
		pcuratom->totalbond = 0 ;
		pcuratom->ndb = 0 ;
		pcuratom->hybrid = 0 ; 
/* default hybridization state is sp3 */
		pcuratom->hybridtext = "(1 2 3)" ;
/* subatom is sp, sp2 or sp3 */
	}
/* initializations */
	
	for(i = 0, pcurbond = tabbond ; i < nbond ; i++, pcurbond++) {
/* loop over bonds */
		pa1 = tabatom + (pcurbond->at1 - 1) ;
		pa2 = tabatom + (pcurbond->at2 - 1) ;
		if(!strcmp(pa1->elt, "H") || !strcmp(pa2->elt, "H")) {
/* skip bonds with an explicit Hydrogen */
			continue ;
		}
		if(pcurbond->bondtype == 3) {
/* triple bonds */
			pa1->hybrid = 2 ;
			pa2->hybrid = 2 ;
			pa1->hybridtext = "1" ;
			pa2->hybridtext = "1" ;
/* atoms that are triply bonded are sp atoms */
		} else if(pcurbond->bondtype == 2) {
/* double bonds */
			pa1->ndb++ ;
			pa2->ndb++ ;
/* double bond count */
		} else {
			if(pcurbond->bondtype != 1) {
/* bond must be either simple, double or triple */
				fprintf(stderr, "%s: invalid bond between atoms "
				"%d and %d (type is %d, only 1, 2 and 3 are allowded.\n", 
				molname, pcurbond->at1, pcurbond->at2, pcurbond->bondtype) ;
				myexit(1) ;
			}
		}
		pa1->totalbond++ ;
		pa2->totalbond++ ;
/* update number of neighbors of each atom */
	}
	
	for(i = 0, pcuratom = tabatom ; i < natom ; i++, pcuratom++) {
/* loop over atoms */
		if(!strcmp(pcuratom->elt, "H")) {
/* skip explicit Hydrogen atoms */
			continue ;
		}
		if(pcuratom->ndb == 2) {
			pcuratom->hybrid = 2 ;
			pcuratom->hybridtext = "1" ;
/* an atom with 2 double bonds is sp */
		} else if(pcuratom->ndb == 1) {
			pcuratom->hybrid = 1 ;
			pcuratom->hybridtext = "(1 2)" ;
/* an atom with 1 double bond is sp2 */
		}
		v = z(pcuratom->elt) ;
/* get atom index in element table */
		if(v < 0) {
/* unknown atom */
			fprintf(stderr, "%s: unknown atom: %s.\n", molname, pcuratom->elt) ;
			myexit(1) ;
		}
		pcuratom->valency = nucl[v].val ;
/* get atom valency */
		pcuratom->maxh = pcuratom->valency - pcuratom->hybrid - pcuratom->totalbond ;
/* calculate number of H in the structure, that is the maximum number of H in the substructure */
		switch (pcuratom->maxh) {
			case 0 :  pcuratom->multtext = "0" ; break ;
			case 1 :  pcuratom->multtext = "(0 1)"; break ;
			case 2 :  pcuratom->multtext = "(0 1 2)"; break ;
			case 3 :  pcuratom->multtext = "(0 1 2 3)"; break ;
			default : fprintf(stderr, "%s: Atom %d has %d bonds.\n", molname, i+1, pcuratom->maxh) ; myexit(1) ;
		}
/* subatom multiplicity */
		if(pcuratom->hybrid == 1 && pcuratom->maxh == 0) {
/* hybridization for sp2 atoms without H atoms: they cannot become sp by H abstraction */
			pcuratom->hybridtext = "2" ;
		}
		if(pcuratom->hybrid == 0 && pcuratom->maxh == 0) {
/* hybridization for sp3 atoms without H atoms: they cannot become sp2 by H abstraction */
			pcuratom->hybridtext = "3" ;
		}
		if(pcuratom->hybrid == 0 && pcuratom->maxh == 1) {
/* hybridization for sp3 atoms with only one H atom: they can become sp2 but not sp by H abstraction */
			pcuratom->hybridtext = "(2 3)" ;
		}
	}

#ifdef DEBUG
	printf("Atom table :\n") ;
	for(i = 0 ; i < natom ; i++) {
		printf("%d %s val%d hyb%d tot%d maxh%d \n", i+1, tabatom[i].elt, tabatom[i].valency, tabatom[i].hybrid, tabatom[i].totalbond, tabatom[i].maxh) ;
	}
#endif
}

/*
 * z returns the index of an atom symbol in global table nucl
 */
int	z(char *nucname)
{
	int	i ;
	for (i = 0 ; (i < MAX_NUCL) && strcmp(nucname, nucl[i].sym) ; i++)  
		;
	return (i == MAX_NUCL) ? (-1) : i ;
}

/*
 * writes a molecule as its corresponding LSD substructure
 */
void writemol(char *dir, int f, char *molname, char *skelname)
{
	FILE *sstr ;
	int i ;
	char dirout[LONGMAX] ;
	atom *pa1, *pa2 ;
	atom *pcuratom ;
	bond *pcurbond ;
	
	sprintf(dirout, "%s/file%d", dir, f) ;
/* substructure file name */
	if((sstr = fopen(dirout, "w")) == NULL) {
		fprintf(stderr, "Cannot open %s for writing.\n", dirout) ;
		myexit(1) ;
	}
/* substructure file opening */
	fprintf(sstr, "; Coded from file: %s\n", molname) ;
	fprintf(sstr, "; Skeleton: %s\n\n", skelname) ;
/* header comments */
	for(i = 0, pcuratom = tabatom ; i < natom ; i++, pcuratom++) {
		if(!strcmp(pcuratom->elt, "H")) {
			continue ;
		}
		fprintf(sstr, "SSTR S%d %s %s %s\n", i+1, pcuratom->elt, pcuratom->hybridtext, pcuratom->multtext) ;
	}
/* SSTR commands for atoms */
	fprintf(sstr, "\n") ;
	for(i = 0, pcurbond = tabbond ; i < nbond ; i++, pcurbond++) {
		pa1 = tabatom + (pcurbond->at1 - 1) ;
		pa2 = tabatom + (pcurbond->at2 - 1) ;
		if(!strcmp(pa1->elt, "H") || !strcmp(pa2->elt, "H")) {
			continue ;
		}		
		fprintf(sstr, "LINK S%d S%d\n", pcurbond->at1, pcurbond->at2) ;
	}
/* LINK commands for bonds */
	fclose(sstr) ;
/* done */
}

/*
 * copymol copies a MOL file in dir, comments are added to show where is comes from
 */
void copymol(char *dir, int f, char *molname, char *skelname) 
{
	FILE *molcopy ;
	char dirmolcopy[LONGMAX] ;
	char buffer[LONGMAX] ;
	int Nblus ;
	
	sprintf(dirmolcopy, "%s/file%d.mol", dir, f) ;
/* new name of MOL file */
	if((molcopy = fopen(dirmolcopy, "w")) == NULL) {
		fprintf(stderr, "Cannot open %s for writing.\n", dirmolcopy) ;
		myexit(1) ;
	}
	fgets(buffer, LONGMAX, mol) ;
/* reads first line in orginal file */
	fprintf(molcopy, "%s", buffer) ;
/* writes a copy of it */
	fgets(buffer, LONGMAX, mol) ;
/* reads the second line */
	fprintf(molcopy, "%s", buffer) ;
/* writes a copy of it */
	fgets(buffer, LONGMAX, mol) ;
/* reads the third line. its content will be ignored */
	fprintf(molcopy, "Copy of file: %s for skeleton %s\n", molname, skelname) ;
/* third line with the original file name and the corresponding skeleton name */
	while((fgets(buffer, LONGMAX, mol)) != 0) {
		fputs(buffer, molcopy) ;
	}
/* copy all following lines */
	fclose(molcopy) ;
/* done */
}
