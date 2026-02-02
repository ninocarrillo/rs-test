/* 
 * File:   rs2.h
 * Author: nino
 *
 * Created on May 17, 2017, 7:34 AM
 */

#ifndef RS2_H
#define	RS2_H

#include "stdint.h"
#include "rs2_def_struct.h"

// init_rs2
// Initializes Reed Solomon parameters in memory for library functions.
// Arg1: starting root for generator polynomial
// Arg2: number of roots in generator polynomial
void InitRS2(int16_t, int16_t, RS2_def_struct*);

// rs2_enc
// Computes Reed Solomon parity symbols for input array and appends them to the
// end of the original array.
// Arg1: pointer to first word of input array
// Arg2: word count of input array
void RSEncode16(uint16_t *, int16_t, RS2_def_struct*);

// rs2_enc
// Computes Reed Solomon parity symbols for input array and appends them to the
// end of the original array.
// Arg1: pointer to first word of input array
// Arg2: word count of input array
void RSEncode8(uint8_t *, int16_t, RS2_def_struct*);

// rs2_dec
// Computes error locations and values based on Reed Solomon parity symbols
// appended to end of input array. Corrects input array if possible.
// Arg1: pointer to first word of input array
// Arg2: word count of input array
// Returns number of errors corrected. Returns negative if correction failed.
// Clobbers input array.
int16_t RSDecode16(uint16_t *, int16_t, RS2_def_struct*);


// rs2_dec
// Computes error locations and values based on Reed Solomon parity symbols
// appended to end of input array. Corrects input array if possible.
// Arg1: pointer to first word of input array
// Arg2: word count of input array
// Returns number of errors corrected. Returns negative if correction failed.
// Clobbers input array.
int16_t RSDecode8(uint8_t *, int16_t, RS2_def_struct*);

// get_genpoly_coef
// Returns generator polynomial coefficient for term x^Arg1.
int16_t rs2_get_genpoly_coef(int16_t, RS2_def_struct*);


#endif	/* RS2_H */

