/* 
 * File:   rs2_def_struct.h
 * Author: nino
 *
 * Created on November 11, 2017, 7:54 AM
 */



#ifndef RS2_DEF_STRUCT_H
#define	RS2_DEF_STRUCT_H

#include "stdint.h"
#include "gf2_def_struct.h"

// Change to suit your application. Bigger = more memory allocated.
#define MAX_GENPOLY_ROOTS 32

typedef struct {
    uint16_t Genpoly[MAX_GENPOLY_ROOTS + 1];
    int16_t FirstRoot;
    int16_t NumRoots;
    int16_t FieldOrder;
    GF2_def_struct *GF;
    int16_t ErrorCount;
    uint16_t ErrorLocs[MAX_GENPOLY_ROOTS];
    uint16_t ErrorMags[MAX_GENPOLY_ROOTS];
    int16_t MinimumErrorDistance;
} RS2_def_struct;

#endif	/* RS2_DEF_STRUCT_H */

