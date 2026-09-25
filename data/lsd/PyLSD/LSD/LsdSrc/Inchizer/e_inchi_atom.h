/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.01
 * July 21, 2006
 * Developed at NIST
 */

#ifndef __INCHI_ATOM_H__
#define __INCHI_ATOM_H__

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

void           e_FreeInchi_Atom( inchi_Atom **at );
void           e_FreeInchi_Stereo0D( inchi_Stereo0D **stereo0D );
inchi_Atom     *e_CreateInchi_Atom( int num_atoms );
inchi_Stereo0D *e_CreateInchi_Stereo0D( int num_stereo0D );
void           e_FreeInchi_Input( inchi_Input *inp_at_data );
int            e_RemoveRedundantNeighbors( inchi_Input *inp_at_data );

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /* __INCHI_ATOM_H__ */
