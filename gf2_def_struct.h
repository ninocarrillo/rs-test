/* 
 * File:   gf2_def_struct.h
 * Author: nino
 *
 * Created on November 11, 2017, 7:42 AM
 */

#ifndef GF2_DEF_STRUCT_H
#define	GF2_DEF_STRUCT_H

#include "stdint.h"

// Change to suit your application. Bigger = more memory used!
#define MAX_GF_BITS 10

// Don't change below this line.
// MAX_FIELD_SIZE controls how much memory is allocated for tables.
#define MAX_FIELD_SIZE (1 << MAX_GF_BITS)

typedef struct {
    int Power;
    int GenPoly;
    int Order;
    int LFSR;
    int Mask;
    int Table[MAX_FIELD_SIZE - 1];
    int Index[MAX_FIELD_SIZE];
    int Inverse[MAX_FIELD_SIZE];
} GF2_def_struct;

#endif	/* GF2_DEF_STRUCT_H */

