/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.01
 * July 21, 2006
 * Developed at NIST
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "e_mode.h"
#include "inchi_api.h"
#include "e_inchi_atom.h"
#include "e_ichisize.h"
#include "e_util.h"
/******************************************************************************************************/
void e_FreeInchi_Atom( inchi_Atom **at )
{
    if ( at && *at ) {
        e_inchi_free( *at );
        *at = NULL;
    }
}
/******************************************************************************************************/
void e_FreeInchi_Stereo0D( inchi_Stereo0D **stereo0D )
{
    if ( stereo0D && *stereo0D ) {
        e_inchi_free( *stereo0D );
        *stereo0D = NULL;
    }
}
/******************************************************************************************************/
inchi_Atom *e_CreateInchi_Atom( int num_atoms )
{
   inchi_Atom *p = (inchi_Atom* ) e_inchi_calloc(num_atoms, sizeof(inchi_Atom) );
   return p;
}
/******************************************************************************************************/
inchi_Stereo0D *e_CreateInchi_Stereo0D( int num_stereo0D )
{
   return (inchi_Stereo0D* ) e_inchi_calloc(num_stereo0D, sizeof(inchi_Stereo0D) );
}
/******************************************************************************************************/
void e_FreeInchi_Input( inchi_Input *inp_at_data )
{
    e_FreeInchi_Atom( &inp_at_data->atom );
    e_FreeInchi_Stereo0D( &inp_at_data->stereo0D );
    memset( inp_at_data, 0, sizeof(*inp_at_data) ); 
}
/*********************************************************/
int e_RemoveRedundantNeighbors( inchi_Input *inp_at_data )
{
    int i, j, k;
    inchi_Atom *at        = inp_at_data->atom;
    int         num_atoms = inp_at_data->num_atoms;
    for ( i = 0; i < num_atoms; i ++ ) {
        for ( j = k = 0; j < at[i].num_bonds; j ++ ) {
            if (  at[i].neighbor[j] < i ) {
                at[i].neighbor[k]    = at[i].neighbor[j];
                at[i].bond_type[k]   = at[i].bond_type[j];
                at[i].bond_stereo[k] = at[i].bond_stereo[j];
                k ++;
            }
        }
        for ( j = k; j < at[i].num_bonds; j ++ ) {
            at[i].neighbor[j]    = 0;
            at[i].bond_type[j]   = 0;
            at[i].bond_stereo[j] = 0;
        }
        at[i].num_bonds = k;
    }
    return 0;
}
