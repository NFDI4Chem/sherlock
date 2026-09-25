/*
    file: LsdSrc/Uniciter/uniciter.h
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


/*
 * The uniciter library is a ultra-thin software layer above the GNU AVL library
 * (see avl.h) that was designed to identify redundant (char *) strings.
 * A uniciter address (void *pu) is returned by new_uniciter().
 * A call to is_unique(void *pu, char *s) returns 1 if strings was
 * never submitted to the uniciter at address pu.
 * An integer status code is returned by uniciter_status():
 *  0: OK
 *  1: memory allocation failed
 *  2: allocation of 0 bytes (internal error)
 *  3: liberation of NULL pointer (internal error)
 * The status should be checked after each uniciter operation.
 * The status is common to all uniciters.
 * The uniciter operation *pu allocates memory that is freed by
 * calling free_uniciter(pu).
 */

extern void	*new_uniciter(void) ;
extern int	is_unique(void *pu, char *s) ;
extern int	uniciter_status(void) ;
extern void	free_uniciter(void *pu) ;
