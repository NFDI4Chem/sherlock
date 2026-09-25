/*
    file: LsdSrc/dyntext.h
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

typedef struct Dynstring
/* is a string of unknown length that may grow
as big as long as ther is free memory left. */
{
	char *string ;
	int length ;
	int szblock ;
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
/* is an array of lines (or text) of unknown length that may grow
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

extern dynstring *dynstring_new(int szblock) ;
extern void dynstring_reset(dynstring *p) ;
extern void dynstring_free(dynstring **pp) ;
extern char *dynstring_getstring(dynstring *p) ;
extern int dynstring_getlength(dynstring *p) ;
extern int dynstring_pushc(dynstring *p, char c) ;
extern int dynstring_pushs(dynstring *p, char *s) ;
extern int dynstring_pushsn(dynstring *p, char *s, int n) ;
extern int dynstring_pushd(dynstring *p, dynstring *q) ;

extern line *line_new(char *s, int n, int *data) ;
extern void line_print(line *p) ;
extern void line_free(line **ppl) ;

extern dyntext *dyntext_new(int szblock) ;
extern void dyntext_reset(dyntext *p) ;
extern void dyntext_free(dyntext **pp) ;
extern line *dyntext_getline(dyntext *p, int i) ;
extern char *dyntext_getstr(dyntext *p, int i) ;
extern int dyntext_getindex(dyntext *p, int i) ;
extern int *dyntext_getserial(dyntext *p, int i) ;
extern void dyntext_setstr(dyntext *p, int i, char *s) ;
extern void dyntext_setindex(dyntext *p, int i, int n) ;
extern void dyntext_setserial(dyntext *p, int i, int *n) ;
extern int dyntext_getlength(dyntext *p) ;
extern int dyntext_push(dyntext *p, dynstring *q, int n, int *data) ;
extern void dyntext_print(dyntext *p) ;
extern char *dyntext_str_from_index(dyntext *p, int n) ;
extern int *dyntext_serial_from_index(dyntext *p, int n) ;
extern int dyntext_index_exists(dyntext *p, int n) ;
