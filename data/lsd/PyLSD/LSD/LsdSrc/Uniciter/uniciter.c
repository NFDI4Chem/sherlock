/*
    file: LsdSrc/Uniciter/uniciter.c
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

#include <stdlib.h>
#include <string.h>
#include "avl.h"

/*
 * see interface description in uniciter.h
 */
void *new_uniciter(void) ;
int is_unique(void *pu, char *s) ;
int uniciter_status(void) ;
void free_uniciter(void *pu) ;

static int comp (const void *pa, const void *pb, void *param) ;
static void free_string(void *s, void *param) ;
static void *my_avl_malloc (struct libavl_allocator *allocator, size_t size) ;
static void my_avl_free (struct libavl_allocator *allocator, void *block) ;

static struct libavl_allocator my_avl_allocator = { my_avl_malloc, my_avl_free } ;
static int status = 0 ;

static int comp (const void *pa, const void *pb, void *param)
/* AVL comparison function for char * strings */
{
	return strcmp (pa, pb) ;
}

static void free_string(void *s, void *param)
{
	free(s) ;
}

static void *my_avl_malloc (struct libavl_allocator *allocator, size_t size)
/* memory allocator that updates the status value */
{
	void *ptr ;
	
	if (size == 0) {
		status = 2 ;
		return NULL ;
	}
	ptr = malloc(size) ;
	if (ptr == NULL) {
		status = 1 ;
	}
	return ptr ;
}

static void my_avl_free (struct libavl_allocator *allocator, void *block)
/* memory liberator that updates the status value */ 
{
	if (block == NULL) {
		status = 3 ;
	} else {
		free (block);
	}
}

void *new_uniciter(void)
/* gets a new AVL tree according to the local comparison and allocation functions */
{
	void *ptr = avl_create(comp, NULL, &my_avl_allocator) ;
	if(ptr == NULL) {
		status = 1 ;
	}
	return ptr ;
}

int is_unique(void *pu, char *s)
/* test for the non-existence of string s in the AVL tree.
appends a new sting to the AVL tree */
{
	struct avl_table *ptab ;
	char *locals ;
	void *exists ;

	ptab = (struct avl_table *)pu ;
	locals = (char *)malloc(strlen(s) + 1) ;
	if(locals == NULL) {
		status = 1 ;
		return 0 ;
	}
	strcpy(locals, s) ;
	exists = avl_insert(ptab, locals) ;
	
	return (exists == NULL) ? 1 : 0 ;
}

int uniciter_status(void)
/* get current global status */
{
	return status ;
}

void free_uniciter(void *pu)
/* free uniciter data */
{
	if(pu == NULL) return ;
	struct avl_table *ptab = (struct avl_table *)pu ;
	avl_destroy(ptab, free_string) ;
}
