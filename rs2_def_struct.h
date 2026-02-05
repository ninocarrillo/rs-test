/* 
 * File:   rs2_def_struct.h
 * Author: nino
 *
 * Created on November 11, 2017, 7:54 AM
 */



#ifndef RS2_DEF_STRUCT_H
#define	RS2_DEF_STRUCT_H

#include "gf2_def_struct.h"

// Change to suit your application. Bigger = more memory allocated.
#define MAX_GENPOLY_ROOTS 32

typedef struct {
    GF2_def_struct *GF;
    int Genpoly[MAX_GENPOLY_ROOTS + 1];
    int ErrorLocs[MAX_GENPOLY_ROOTS];
    int ErrorMags[MAX_GENPOLY_ROOTS];
    int FirstRoot;
    int NumRoots;
    int FieldOrder;
    int ErrorCount;
} RS2_def_struct;

#endif	/* RS2_DEF_STRUCT_H */

