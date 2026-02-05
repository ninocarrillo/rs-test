#include "rs2.h"
#include "gf2.h"
#include <stdio.h>

void InitRS2(int first_root, int num_roots, RS2_def_struct *rs) {
    rs->FirstRoot = first_root;
    rs->NumRoots = num_roots;
    rs->FieldOrder = GF2GetOrder(rs->GF);
    // Generate Reed Solomon generator polynomial through convolution of polynomials.
    // rs->GenPoly = (x + a^b)(x + a^b+1)...(x + a^b+r-1)
    // start with rs->GenPoly = x + a^b
    // lowest order coefficient in lowest index of array
	// b represents the "first consecutive root" of generator polynomial.
    rs->Genpoly[0] = GF2Pow(rs->FirstRoot, rs->GF);
    rs->Genpoly[1] = 1;
    int factorpoly[2];
    // preload the x^1 coefficient in the factor polynomial
    factorpoly[1] = 1;
    for (int i = 1; i < num_roots; i++) {
        factorpoly[0] = GF2Pow(i + rs->FirstRoot, rs->GF);
        GF2Conv(&rs->Genpoly[0], i + 1, &factorpoly[0], 2, rs->GF);
    }
}

void RSEncode(int *message, int message_size, RS2_def_struct *rs) {
	for (int i = 0; i < rs->NumRoots; i++) {
		message[i + message_size] = 0;
	}
	int quotient[MAX_GENPOLY_ROOTS + 1];
	for (int i = 0; i < rs->NumRoots + 1; i++) {
		quotient[i] = message[i];
	}
	for (int i = 0; i < message_size; i++) {
		int x = quotient[0];
		for (int j = 1; j < rs->NumRoots + 1; j++) {
			quotient[j - 1] = GF2Mul(x, rs->Genpoly[rs->NumRoots - j], rs->GF) ^ quotient[j];
		}
		quotient[rs->NumRoots] = message[i + rs->NumRoots + 1];
	}
	for (int i = 0; i < rs->NumRoots; i++) {
		message[i + message_size] = quotient[i];
	}
}

int calc_syndromes(RS2_def_struct *rs, int block_size, int *data_block, int *syndromes) {
	// Calculate one syndrome for each root of rs->GenPoly.
	// Each syndrome is the evaluation of the message polynomial at a root of rs->GenPoly
	int nonzero = 0; // Count how many non-zero syndromes are calculated.
	for (int i = 0; i < rs->NumRoots; i++) {
		syndromes[i] = 0;
        int x = GF2Pow(rs->FirstRoot + i, rs->GF);
		for (int j = 0; j < block_size - 1; j++) {
			syndromes[i] = GF2Mul(syndromes[i] ^ data_block[j], x, rs->GF);
		}
		syndromes[i] = syndromes[i] ^ data_block[block_size - 1];
		if (syndromes[i]) {
			nonzero++;// Count how many non-zero syndromes are calculated.
		}
	}
	return nonzero;
}

void calc_berlekamp(RS2_def_struct *rs, int *syndromes, int *error_locator_poly) {
    int next_error_locator_poly[MAX_GENPOLY_ROOTS];
    int correction_poly[MAX_GENPOLY_ROOTS];
	error_locator_poly[0] = 1;
	int order_tracker = 0;
	for (int i = 0; i < rs->NumRoots; i++) {
		next_error_locator_poly[i] = 0;
		correction_poly[i] = 0;
	}
	correction_poly[1] = 1;
	for (int step_factor = 1; step_factor <= rs->NumRoots; step_factor++) {
        // first calculate error value, e
		int y = step_factor - 1;
		int e = syndromes[y];
		for (int i = 1; i <= order_tracker; i++) {
			int x = y - i;
			e = e ^ GF2Mul(error_locator_poly[i], syndromes[x], rs->GF);
		}
        // now update the estimate of the error locator polynomial
		if (e != 0) {
			for (int i = 0; i <= order_tracker; i++) {
				next_error_locator_poly[i] = error_locator_poly[i] ^ GF2Mul(e, correction_poly[i], rs->GF);
			}
            // and update the value of the correction polynomial
			e = GF2Inv(e, rs->GF);
			for (int i = 0; i <= rs->NumRoots; i++) {		// refine loop limit?
				correction_poly[i] = GF2Mul(error_locator_poly[i], e, rs->GF);
			}
			for (int i = 0; i <= rs->NumRoots; i++) {	// refine loop limit?
				error_locator_poly[i] = next_error_locator_poly[i];         // this should only happen if e != 0!
			}
		}
		if ((2 * order_tracker) < step_factor) {
			order_tracker = step_factor - order_tracker;
		}
        // multiply error_value_poly by x (increase power by one)
		for (int i = rs->NumRoots; i > 0; i--) {
			correction_poly[i] = correction_poly[i - 1];
		}
		correction_poly[0] = 0;
	}
}

int calc_chien(RS2_def_struct *rs, int block_size, int *error_locator_poly) {
	// Calculate error locations and error count from error locator polynomial.
	// Brute force search for roots of error locator polynomial. Solutions
	// found when polynomial evaluates to zero.
	rs->ErrorCount = 0;
	for (int j = 0; j < block_size; j++) {
		int x = 0;
		int y = j + rs->FieldOrder - block_size;  // account for code shortening by modifying index variable
		for (int i = 1; i <= rs->NumRoots/2; i++) {
			if (error_locator_poly[i]) {
				int z = y * i;  // calculate power by multiplying exponents
				z = z + GF2Log(error_locator_poly[i], rs->GF); //multiply by adding exponents
				while (z > (rs->FieldOrder - 2)) {
					z = z - rs->FieldOrder + 1;
				}
				x = x ^ GF2Pow(z, rs->GF);
			}
		}
		x = x ^ error_locator_poly[0];
		if (x == 0) {
			rs->ErrorLocs[rs->ErrorCount] = j;
			rs->ErrorCount++;
            // Todo: check for an ambiguous solution
            
		}
	}
	return rs->ErrorCount;
}

void calc_error_value_poly(RS2_def_struct *rs, int *syndromes, int *error_locator_poly, int *error_value_poly) {
	for (int i = 0; i < rs->ErrorCount; i++) {
		error_value_poly[i] = syndromes[i];
		for (int j = 1; j <= i; j++) {
			error_value_poly[i] = error_value_poly[i] ^ GF2Mul(syndromes[i - j], error_locator_poly[j], rs->GF);
		}
	}
}

void calc_forney(RS2_def_struct *rs, int block_size, int *error_locator_poly, int *error_value_poly, int *syndromes) {
	// Forney algorithm to determine error values

	int e, x, y, z;
	
	for (int i = 0; i < rs->ErrorCount; i++) { // compute an error value for each error location
		e = block_size - (rs->ErrorLocs[i] + 1);
		z = error_value_poly[0];
		for (int j = 1; j < rs->ErrorCount; j++) { // calculate numerator
			x = GF2Clamp((rs->FieldOrder - 1) - (e * j), rs->GF);
			z ^= GF2Mul(error_value_poly[j], GF2Pow(x, rs->GF), rs->GF);
		}
		z = GF2Mul(z, GF2Pow(e, rs->GF), rs->GF);

		y = error_locator_poly[1];
		for (int j = 3; j <= rs->NumRoots / 2; j = j + 2) { 
			x = GF2Mod(e * (j - 1), rs->GF);
			x = rs->FieldOrder - (x + 1);
			y = y ^ GF2Mul(error_locator_poly[j], GF2Pow(x, rs->GF), rs->GF);
		}
		y = GF2Pow(GF2Clamp(rs->FieldOrder - (GF2Log(y, rs->GF) + 1), rs->GF), rs->GF);
		rs->ErrorMags[i] = GF2Mul(y, z, rs->GF);
	}
}

int RSDecode(int *data_block, int block_size, RS2_def_struct *rs) {
	int syndromes[MAX_GENPOLY_ROOTS];
    int error_locator_poly[MAX_GENPOLY_ROOTS];
    int error_value_poly[MAX_GENPOLY_ROOTS + 1];
    
    // clear out arrays
	for (int i = 0; i <  rs->NumRoots; i++) {
		error_locator_poly[i] = 0;
		error_value_poly[i] = 0;
	}
	error_value_poly[rs->NumRoots] = 0;
	
	calc_syndromes(rs, block_size, data_block, syndromes);
	
	calc_berlekamp(rs, syndromes, error_locator_poly);
	
	calc_chien(rs, block_size, error_locator_poly);

	calc_error_value_poly(rs, syndromes, error_locator_poly, error_value_poly);

	calc_forney(rs, block_size, error_locator_poly, error_value_poly, syndromes);

	for (int i = 0; i < rs->ErrorCount; i++) {
		// Correct each detected error
		data_block[rs->ErrorLocs[i]] = data_block[rs->ErrorLocs[i]] ^ rs->ErrorMags[i];
	}

	// check for success by calculating syndromes (should be zero if no errors)
	int nonzero = calc_syndromes(rs, block_size, data_block, syndromes);
	
	if (nonzero) {
		return -nonzero;
	}
	return rs->ErrorCount; // return number of errors corrected    
}