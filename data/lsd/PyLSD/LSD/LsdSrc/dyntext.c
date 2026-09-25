/*
    file: LsdSrc/dyntext.c
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
#include <stdlib.h>
#include <string.h>
#include "lvm.h"

typedef struct Dynstring
/* is a string of unknown length that may grow
as big as long as there is free memory left. */
{
/* the string itself */
	char *string ;
/* the string length, final \0 included.
length is 1 for an empty chain */
	int length ;
/* size in bytes of dynamically allocated memory blocks */
	int szblock ;
/* number of already allocated memory blocks */
	int nbblock ;
} dynstring ;

typedef struct Line
/* a C string, an attached index, and an attached pointer */
{
/* a fragment filename */
	char *str ;
/* the external fragment index */
	int index ;
/* the serialized fragment data */
	int *serial ;
} line ;

typedef struct Dyntext
/* is an array of lines (C string and int value) of unknown length that may grow
as big as long as there is free memory left. */
{
/* the text itself */
	line **lines ;
/* the text length, final NULL included.
length is 1 for an empty text */
	int length ;
/* size in bytes of dynamically allocated memory blocks */
	int szblock ;
/* number of already allocated memory blocks */
	int nbblock ;
} dyntext ;

dynstring *dynstring_new(int szblock) ;
static void dynstring_lengthen(dynstring *p) ;
void dynstring_reset(dynstring *p) ;
void dynstring_free(dynstring **pp) ;
char *dynstring_getstring(dynstring *p) ;
int dynstring_getlength(dynstring *p) ;
int dynstring_pushc(dynstring *p, char c) ;
int dynstring_pushs(dynstring *p, char *s) ;
int dynstring_pushsn(dynstring *p, char *s, int n) ;
int dynstring_pushd(dynstring *p, dynstring *q) ;

line *line_new(char *s, int n, int *data) ;
void line_free(line **ppl) ;
void line_print(line *pl) ;

dyntext *dyntext_new(int szblock) ;
static void dyntext_lengthen(dyntext *p) ;
void dyntext_reset(dyntext *p) ;
void dyntext_free(dyntext **pp) ;
line *dyntext_getline(dyntext *p, int i) ;
char *dyntext_getstr(dyntext *p, int i) ;
int dyntext_getindex(dyntext *p, int i) ;
int *dyntext_getserial(dyntext *p, int i) ;
void dyntext_setstr(dyntext *p, int i, char *s) ;
void dyntext_setindex(dyntext *p, int i, int n) ;
void dyntext_setserial(dyntext *p, int i, int *n) ;
int dyntext_getlength(dyntext *p) ;
int dyntext_push(dyntext *p, dynstring *q, int n, int *data) ;
void dyntext_print(dyntext *p) ;
static line *dyntext_line_from_index(dyntext *p, int n) ;
char *dyntext_str_from_index(dyntext *p, int n) ;
int *dyntext_serial_from_index(dyntext *p, int n) ;
int dyntext_index_exists(dyntext *p, int n) ;

dynstring *dynstring_new(int szblock)
/* creates a new dynamic string, with memory allocation
when szblock characters are pushed in. */
{
	dynstring *p = NULL ;

/* checks block size validity */
	if(szblock < 1) return p ;
/* allocation of the dynstring object itself */
	p = (dynstring *)malloc(sizeof(dynstring)) ;
	if(p == NULL) return p ;
/* object members initialization */
	p->length = 0 ;
	p->string = NULL ;
	p->szblock = szblock ;
	p->nbblock = 0 ;
/* allocates the first memory block */
	dynstring_lengthen(p) ;
/* sets string to the empty string */
	dynstring_reset(p) ;

	return p ;
}

static void dynstring_lengthen(dynstring *p)
/* allocates one more memory block. */
{
	p->nbblock++ ;
	p->string = (char *)realloc(p->string, p->nbblock * p->szblock * sizeof(char)) ;
}

void dynstring_reset(dynstring *p)
/* sets string to the empty string. */
{
	if(p->string == NULL) return ;
	p->string[0] = '\0' ;
	p->length = 1 ;
}

void dynstring_free(dynstring **pp)
/* frees a dynamic string. */
{
	dynstring *p = *pp ;
	char *s ;

	if(p == NULL) return ;
	s = p->string ;
/* frees string */
	if(s != NULL) free((char *)s) ;
/* frees the dynamic string itself */
	free((char *)p) ;
	*pp = NULL ;
}

char *dynstring_getstring(dynstring *p)
/* returns the string address of a dynamic string */
{
	if(p == NULL) return NULL ;
	return p->string ;
}

int dynstring_getlength(dynstring *p)
/* returns the string length of a dynamic string,
the lentgth of an empty string is zero. */
{
	if(p == NULL) return 0 ;
	return p->length - 1 ;
}

int dynstring_pushc(dynstring *p, char c)
{
/* appends a character to a dynamic string,
returns 1 if OK, 0 otherwise. */
	char *s ;
	int l ;

	if(p == NULL || p->string == NULL) return 0 ;
	l = p->length ;
	if(l / p->szblock == p->nbblock) {
/* the last allocated memory block is full, requires one more */
		dynstring_lengthen(p) ;
	}
	s = p->string ;
	if(s == NULL) return 0 ;
/* strores the new character */
	s[l - 1] = c ;
/* makes a valid C string of string */
	s[l] = '\0' ;
/* updates length */
	p->length++ ;
	return 1 ;
}

int dynstring_pushs(dynstring *p, char *s)
/* appends a C string to a dynamic string,
returns 1 if OK, 0 otherwise. */
{
	char *pc ;
	char c ;
	int r ;

	for(pc = s ; (c = *pc) != '\0' && (r = dynstring_pushc(p, c)) ; pc++) ;
	return r ;
}

int dynstring_pushsn(dynstring *p, char *s, int n)
/* appends a C string to a dynamic string,
returns 1 if OK, 0 otherwise. 
the number of pushed characters cannot be more than n */
{
	char *pc ;
	char c ;
	int r ;
	int i ;

	for(pc = s, i = 0 ; i < n && (c = *pc) != '\0' && (r = dynstring_pushc(p, c)) ; pc++, i++) ;
	return r ;
}

int dynstring_pushd(dynstring *p, dynstring *q)
/* appends a dynamic sting q to dynamic string p,
returns 1 if OK, 0 otherwise. */
{
	if(q == NULL) return 0 ;
	return dynstring_pushs(p, q->string) ;
}

/********************************************************/

line *line_new(char *s, int n, int *data)
/* creates a line from a C string s, and an integer value n */
{
	line *pl ;
	char *news ;

	if(s == NULL) return NULL ;
/* creates a line */
	pl = (line *)malloc(sizeof(line)) ;
	if(pl == NULL) return NULL ;
/* creates the C string */
	news = (char *)malloc((strlen(s)+1) * sizeof(char)) ;
	if(news == NULL) {
		free(pl) ;
		return NULL ;
	}
/* saves a copy of the string */
	strcpy(news, s) ;
/* sets line attributes */
	pl->str = news ;
	pl->index = n ;
	pl->serial = data ;

	return pl ;
}

void line_free(line **ppl)
/* frees a line. ppl is the address of the pointer to the line. */
{
	line *pl = *ppl ;
	char *s ;
	int *data ;

	if(pl == NULL || (s = pl->str) == NULL) return ;
	data = pl->serial ;
	if(data != NULL) {
		if(pl->index == -1) {
/* data is a pointer to a lvm */
			lvm_clean((plvm)data) ;
		} else {
/* data is a pointer to a serialized sub-structure */
			free(data) ;
		}
	}
	free(s) ;
	free((char *)pl) ;
	*ppl = NULL ;
}

void line_print(line *pl)
/* prints a line */
{
	char *s ;

	if (pl == NULL) {
		printf("%s", "null_line") ;
	} else {
		s = pl->str ;
		if(s == NULL) {
			printf("%s in line at %p", "null_string", (void *)pl) ;
		} else {
			printf("(line at %p)%d:%s(string at %p) serial:%p",
				(void *)pl, pl->index, pl->str, 
				(void *)pl->str, (void *)pl->serial) ;
		}
	}
}

/********************************************************/

dyntext *dyntext_new(int szblock)
/* creates a new dynamic text, with memory allocation
when szblock characters are pushed in. */
{
	dyntext *p = NULL ;

/* checks block size validity */
	if(szblock < 1) return p ;
/* allocation of the dynstring object itself */
	p = (dyntext *)malloc(sizeof(dyntext)) ;
	if(p == NULL) return p ;
/* object members initialization */
	p->length = 0 ;
	p->lines = NULL ;
	p->szblock = szblock ;
	p->nbblock = 0 ;
/* allocates the first memory block */
	dyntext_lengthen(p) ;
/* sets string to the empty string */
	dyntext_reset(p) ;
#ifdef DEBUG
dyntext_print(p) ;
#endif
	return p ;
}

static void dyntext_lengthen(dyntext *p)
/* allocates one more memory block. */
{
	p->nbblock++ ;
	p->lines = (line **)
		realloc(p->lines, p->nbblock * p->szblock * sizeof(line *)) ;
}

void dyntext_reset(dyntext *p)
/* sets text to the empty text,
the one with a single null line. */
{
	line **ppl ;

	if(p == NULL || (ppl = p->lines) == NULL) return ;
	ppl[0] = NULL ;
	p->length = 1 ;
}

void dyntext_free(dyntext **pp)
/* frees a dynamic text. */
{
	dyntext *p = *pp ;
	int i ;
	int l ;
	char *s ;
	line **ppl ;
	line *pl ;

	if(p == NULL || (ppl = p->lines) == NULL) return ;
/* frees lines */
	l = p->length - 1 ;
	for(i=0 ; i<l ; i++) {
/* frees a line */
		line_free(&(ppl[i])) ;
	}
/* frees the dynamic text itself */
	free((char *)ppl) ;
	free((char *)p) ;
	*pp = NULL ;
}

line *dyntext_getline(dyntext *p, int i)
/* returns the address of line i of a dynamic text. */
{
	line **ppl ;
	line *pl ;
	int l ;

	if(p == NULL) return NULL ;
	l = p->length - 2 ;
	if(i>l || (ppl = p->lines) == NULL) return NULL ;
	pl = ppl[i] ;
	return pl ;
}

char *dyntext_getstr(dyntext *p, int i)
/* returns the address of line i of a dynamic text. */
{
	line *pl ;

	return ((pl = dyntext_getline(p, i)) == NULL) ? NULL : pl->str ;
}

int dyntext_getindex(dyntext *p, int i)
/* returns the index of string at line i of a dynamic text, -1 if not defined */
{
	line *pl ;

	return ((pl = dyntext_getline(p, i)) == NULL) ? -1 : pl->index ;
}

int *dyntext_getserial(dyntext *p, int i)
/* returns the pointer to serial data at line i of a dynamic text,
-1 if not defined */
{
	line *pl ;

	return ((pl = dyntext_getline(p, i)) == NULL) ? NULL : pl->serial ;
}

void dyntext_setstr(dyntext *p, int i, char *s)
/* sets to s the characters of line i of a dynamic text */
{
	line *pl ;
	char *olds, *news ;

	if((pl = dyntext_getline(p, i)) != NULL) {
		olds = pl->str ;
		if(olds != NULL) free(olds) ;
		news = (char *)malloc(strlen(s) * (sizeof(char) + 1)) ;
		if(news == NULL) return ;
		strcpy(news, s) ;
		pl->str = news ;
	}
}

void dyntext_setindex(dyntext *p, int i, int n)
/* sets to n the index of line i of a dynamic text */
{
	line *pl ;

	if((pl = dyntext_getline(p, i)) != NULL) pl->index = n ;
}

void dyntext_setserial(dyntext *p, int i, int *n)
/* sets to n the pointer to serial data of line i of a dynamic text */
{
	line *pl ;

	if((pl = dyntext_getline(p, i)) != NULL) pl->serial = n ;
}

int dyntext_getlength(dyntext *p)
/* returns the length of a dynamic text,
the lentgth of an empty text is zero,
the length of a null text is -1. */
{
	if(p == NULL) return -1 ;
	return p->length - 1 ;
}

int dyntext_push(dyntext *p, dynstring *q, int n, int *data)
/* appends a line as dynstring q and its index n to a dynamic text p,
returns the line number if OK, -1 otherwise.
Line numbers start at 0. */
{
	int l ;
	char *qs, *news ;
	line **ppl, *newl ;

	if(p == NULL || p->lines == NULL) return -1 ;
	if(q == NULL || (qs = q->string) == NULL) return -1 ;
	l = p->length ;
	if(l / p->szblock == p->nbblock) {
/* the last allocated memory block is full, requires one more */
		dyntext_lengthen(p) ;
	}
/* creates the new string */
	newl = line_new(qs, n, data) ;
#ifdef DEBUG
line_print(newl) ;
#endif
	if(newl == NULL) return -1 ;
/* updates lines */
	ppl = p->lines ;
	ppl[l - 1] = newl ;
/* ensures the list of pointers is NULL terminated */
	ppl[l] = NULL ;
/* updates length */
	p->length++ ;

#ifdef DEBUG
dyntext_print(p) ;
#endif
	return l - 1 ;
}

void dyntext_print (dyntext *p)
/* prints a dynamic text, for debugging */
{
	line **ppl ;
	int i ;
	int l ;

	if(p == NULL) {
		printf("null dyntext\n") ;
	} else {
		printf("(dyntext at %p)", (void *)p) ;
		ppl = p->lines ;
		if(ppl == NULL) {
			printf("%s", "null line array\n") ;
		} else {
			l = p->length ;
			printf("(lines at %p)(%d line(s))\n", (void *)ppl, l) ;
			printf("(%d block(s) of %d line(s))\n", p->nbblock, p->szblock) ;
			for(i=0 ; i<l ; i++) {
				printf("%d:", i) ;
				line_print(ppl[i]) ;
				printf("\n") ;
			}
		}
	}
}

static line *dyntext_line_from_index(dyntext *p, int n)
/* returns the first line whose index is n in dynamic text p. */
{
	line **ppl, *pl ;

	if(p == NULL || (ppl = p->lines) == NULL) return NULL ;
	while((pl = *ppl) != NULL && pl->index != n) ppl++ ;
	return (pl == NULL) ? NULL : pl ;
}

char *dyntext_str_from_index(dyntext *p, int n)
/* returns the first string whose index is n in dynamic text p. */
{
	line *pl = dyntext_line_from_index(p, n) ;
	return (pl == NULL) ? NULL : pl->str ;
}

int *dyntext_serial_from_index(dyntext *p, int n)
/* returns the first serialized data whose index is n in dynamic text p. */
{
	line *pl = dyntext_line_from_index(p, n) ;
	return (pl == NULL) ? NULL : pl->serial ;
}

int dyntext_index_exists(dyntext *p, int n)
{
	return dyntext_line_from_index(p, n) != NULL ;
}
