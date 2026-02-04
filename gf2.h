/*
 * File:   GF2.h
 * Author: nino
 *
 * Created on May 14, 2017, 6:00pm
 */

#ifndef GF2_H
#define	GF2_H

#include "stdint.h"
#include "gf2_def_struct.h"

// init_gf2
// Initializes Galois Field tables in memory to support library functions.
// Arg1: Power of base-2 field, ie 2^Power. 2^8 for 8-bit field.
// Arg2: Generator (reducing) polynomial, in binary. 285 = x^8+x^4+x^3+x^2+1
// Returns 0 for maximal field, number of repetition cycles otherwise.
int InitGF2(int, GF2_def_struct*);

// gf2_get_order
// Returns order (size) of specified field.
int gf2_get_order(GF2_def_struct*);

// gf2_mul
// Performs Galois Field multiplication by addition of exponents.
// Arg1: multiplicand
// Arg2: multiplier
// Returns product.
int gf2_mul(int, int, GF2_def_struct*);

// gf2_div
// Performs Galois Field division by subtraction of exponents.
// Arg1: dividend
// Arg2: divisor
// Returns quotient.
int gf2_div(int, int, GF2_def_struct*);

// gf2_conv
// Convolves two Galois Field polynomials.
// Arg1: pointer to least significant coefficient of polynomial 1
// Arg2: number of coefficients in polynomial 1
// Arg3: pointer to least significant coefficient of polynomial 2
// Arg4: number of coefficients in polynomial 2
// Returns convolved polynomial of length (p1n + p2n - 1).
int gf2_conv(int*, int, int*, int, GF2_def_struct*);

// gf2_pow
// Returns field primitive (2) raised to Arg.
int gf2_pow(int, GF2_def_struct*);

// gf2_log
// Returns log base <field primitive> (2) of Arg.
int gf2_log(int, GF2_def_struct*);

// gf2_inv
// Returns 1/Arg.
int gf2_inv(int, GF2_def_struct*);

#endif	/* GF2_H */
