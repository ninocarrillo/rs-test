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

int calc_syndromes(RS2_def_struct *rs) {
	// Calculate one syndrome for each root of rs->GenPoly.
	// Each syndrome is the evaluation of the message polynomial at a root of rs->GenPoly
	int nonzero = 0; // Count how many non-zero syndromes are calculated.
	for (int i = 0; i < rs->NumRoots; i++) {
		rs->Syndromes[i] = 0;
        int x = GF2Pow(rs->FirstRoot + i, rs->GF);
		for (int j = 0; j < rs->BlockSize - 1; j++) {
			rs->Syndromes[i] = GF2Mul(rs->Syndromes[i] ^ rs->DataBlock[j], x, rs->GF);
		}
		rs->Syndromes[i] = rs->Syndromes[i] ^ rs->DataBlock[rs->BlockSize - 1];
		if (rs->Syndromes[i]) {
			nonzero++;// Count how many non-zero syndromes are calculated.
		}
	}
	return nonzero;
}

void save_syndromes(RS2_def_struct *rs) {
	for (int i = 0; i < rs->NumRoots; i++) {
		rs->SavedSyndromes[i] = rs->Syndromes[i];
	}
}

void calc_berlekamp(RS2_def_struct *rs) {
    int next_poly[MAX_GENPOLY_ROOTS];
    int correction_poly[MAX_GENPOLY_ROOTS];
	for (int i = 0; i <  rs->NumRoots; i++) {
		rs->ErrorLocatorPoly[i] = 0;
	}	
	rs->ErrorLocatorPoly[0] = 1;
	int order_tracker = 0;
	for (int i = 0; i < rs->NumRoots; i++) {
		next_poly[i] = 0;
		correction_poly[i] = 0;
	}
	correction_poly[1] = 1;
	for (int step_factor = 1; step_factor <= rs->NumRoots; step_factor++) {
        // first calculate error value, e
		int y = step_factor - 1;
		int e = rs->Syndromes[y];
		for (int i = 1; i <= order_tracker; i++) {
			int x = y - i;
			e = e ^ GF2Mul(rs->ErrorLocatorPoly[i], rs->Syndromes[x], rs->GF);
		}
        // now update the estimate of the error locator polynomial
		if (e != 0) {
			for (int i = 0; i <= order_tracker; i++) {
				next_poly[i] = rs->ErrorLocatorPoly[i] ^ GF2Mul(e, correction_poly[i], rs->GF);
			}
            // and update the value of the correction polynomial
			e = GF2Inv(e, rs->GF);
			for (int i = 0; i <= rs->NumRoots; i++) {		// refine loop limit?
				correction_poly[i] = GF2Mul(rs->ErrorLocatorPoly[i], e, rs->GF);
			}
			for (int i = 0; i <= rs->NumRoots; i++) {	// refine loop limit?
				rs->ErrorLocatorPoly[i] = next_poly[i];         // this should only happen if e != 0!
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

int calc_chien(RS2_def_struct *rs) {
	// Calculate error locations and error count from error locator polynomial.
	// Brute force search for roots of error locator polynomial. Solutions
	// found when polynomial evaluates to zero.
	rs->ErrorCount = 0;
	// Step through each index position in the code block
	for (int candidate_location = 0; candidate_location < rs->BlockSize; candidate_location++) {
		int evaluation = rs->ErrorLocatorPoly[0];
		// account for code shortening by modifying candidate_root based on block size and field order:
		int candidate_root = (candidate_location + rs->FieldOrder) - rs->BlockSize;  
		for (int i = 1; i <= rs->NumRoots/2; i++) {
			if (rs->ErrorLocatorPoly[i]) {
				// Calculate power by multiplying exponents, then multiply by adding exponents:
				int x = (candidate_root * i) + GF2Log(rs->ErrorLocatorPoly[i], rs->GF);
				// Sum the evaluation, xor is addition in GF
				evaluation ^= GF2Pow(GF2Mod(x, rs->GF), rs->GF);
			}
		}
		// If evaluation is zero, we have found a root of the error locator polynomial.
		if (evaluation == 0) {
			rs->ErrorIndices[rs->ErrorCount] = candidate_location;
			rs->ErrorLocatorRoots[rs->ErrorCount] = candidate_root;
			rs->ErrorCount++;
            // Todo: check for an ambiguous solutions
		}
	}
	return rs->ErrorCount;
}

void calc_error_value_poly(RS2_def_struct *rs) {
	for (int i = 0; i < rs->ErrorCount; i++) {
		rs->ErrorMagPoly[i] = rs->Syndromes[i];
		for (int j = 1; j <= i; j++) {
			rs->ErrorMagPoly[i] ^= GF2Mul(rs->Syndromes[i - j], rs->ErrorLocatorPoly[j], rs->GF);
		}
	}
}

void calc_forney(RS2_def_struct *rs) {
	// Forney algorithm to determine error values
	int denominator, numerator;
	for (int i = 0; i < rs->ErrorCount; i++) {
		// compute an error value for each error location
		// Divide the error value polynomial by the derivitave of the error locator polynomial,
		// both evaluated at the root of the error locator polynomial corresponding to the error location.
		printf("\r\n                     --------- e: %i", rs->ErrorLocatorRoots[i]);
		numerator = rs->ErrorMagPoly[0];
		for (int j = 1; j < rs->ErrorCount; j++) { // calculate numerator
			numerator ^= GF2Mul(rs->ErrorMagPoly[j], GF2Pow(GF2Mod(rs->ErrorLocatorRoots[i] * j, rs->GF), rs->GF), rs->GF);
		}
		// Apply adjustment for first consecutive root:
		numerator = GF2Mul(numerator, GF2Pow(GF2Mod(-rs->ErrorLocatorRoots[i], rs->GF), rs->GF), rs->GF);
		
		denominator = rs->ErrorLocatorPoly[1];
		for (int j = 3; j <= rs->NumRoots / 2; j += 2) {
			denominator ^= GF2Mul(rs->ErrorLocatorPoly[j], GF2Pow(GF2Mod(rs->ErrorLocatorRoots[i] * (j - 1), rs->GF), rs->GF), rs->GF);
		}
		
		// Take inverse of denominator term so division becomes multiplication.
		rs->ErrorMags[i] = GF2Mul(GF2Inv(denominator, rs->GF), numerator, rs->GF);
	}
}

int RSDecode(int *data_block, int block_size, RS2_def_struct *rs) {
	rs->BlockSize = block_size;
	rs->DataBlock = data_block;
    
	calc_syndromes(rs);
	
	calc_berlekamp(rs);
	
	calc_chien(rs);

	calc_error_value_poly(rs);

	calc_forney(rs);

	for (int i = 0; i < rs->ErrorCount; i++) {
		// Correct each detected error
		data_block[rs->ErrorIndices[i]] = data_block[rs->ErrorIndices[i]] ^ rs->ErrorMags[i];
	}

	save_syndromes(rs);

	// check for success by calculating syndromes (should be zero if no errors)
	int nonzero = calc_syndromes(rs);
	
	if (nonzero) {
		return -nonzero;
	}
	return rs->ErrorCount; // return number of errors corrected    
}