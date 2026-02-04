#include "rs2.h"
#include "gf2.h"
#include <stdio.h>

void InitRS2(int first_root, int num_roots, RS2_def_struct *rs) {
    rs->FirstRoot = first_root;
    rs->NumRoots = num_roots;
    rs->FieldOrder = GF2GetOrder(rs->GF);
    // Generate Reed Solomon generator polynomial through convolution of polynomials.
    // rs->genpoly = (x + a^b)(x + a^b+1)...(x + a^b+r-1)
    // start with rs->genpoly = x + a^b
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
	int quotient[MAX_GENPOLY_ROOTS + 1];
	for (int i = 0; i < rs->NumRoots; i++) {
		message[i + message_size] = 0;
	}
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

int calc_syndromes(int *data_block, int block_size, int *syndromes, RS2_def_struct *rs) {
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

int RSDecode(int *data_block, int block_size, RS2_def_struct *rs) {
	int syndromes[MAX_GENPOLY_ROOTS];
    int error_locator[MAX_GENPOLY_ROOTS];
    int correction_poly[MAX_GENPOLY_ROOTS + 1];
    int next_error_locator[MAX_GENPOLY_ROOTS];
	int step_factor, order_tracker;
	int x, y, z, e;
    
	calc_syndromes(data_block, block_size, syndromes, rs);

    // Berlekamp's Algorothm
    // calculate the error locator
    // uses:
    // correction polynomial C
    // step parameter K
    // order tracker L
    // error value x
    
    // clear out arrays
	for (int i = 0; i <  rs->NumRoots; i++) {
		error_locator[i] = 0;
		correction_poly[i] = 0;
		next_error_locator[i] = 0;
	}
	correction_poly[rs->NumRoots] = 0;
	
	
	error_locator[0] = 1;
	correction_poly[1] = 1;	// C = x
	order_tracker = 0;
	for (step_factor = 1; step_factor <= rs->NumRoots; step_factor++) {
        // first calculate error value, e
		y = step_factor - 1;
		e = syndromes[y];
		for (int i = 1; i <= order_tracker; i++) {
			x = y - i;
			e = e ^ GF2Mul(error_locator[i], syndromes[x], rs->GF);
		}
        // now update the estimate of V[x]
		if (e != 0) {
			for (int i = 0; i <= order_tracker; i++) {
				next_error_locator[i] = error_locator[i] ^ GF2Mul(e, correction_poly[i], rs->GF);
			}
            // and update the value of C
			e = GF2Inv(e, rs->GF);
			for (int i = 0; i <= rs->NumRoots / 2; i++) {		// need to refine loop limit t, could be based on L?
				correction_poly[i] = GF2Mul(error_locator[i], e, rs->GF);
			}
			for (int i = 0; i <= rs->NumRoots / 2; i++) {	// refine loop limit
				error_locator[i] = next_error_locator[i];         // this should only happen if e != 0!
			}
		}
		if ((2 * order_tracker) < step_factor) {
			order_tracker = step_factor - order_tracker;
		}
        // multiply C(x) by x (increase power by one)
		for (int i = rs->NumRoots; i > 0; i--) {
			correction_poly[i] = correction_poly[i - 1];
		}
		correction_poly[0] = 0;
	}
	
	
    // now solve the error locator polynomial to find the error positions
    // by using the Chien Search
    // Uses V[i] calculated above
	rs->ErrorCount = 0;  // rs->error_count is the number of errors found
	for (int j = 0; j < block_size; j++) {
		x = 0;
		y = j + rs->FieldOrder - block_size;  // account for code shortening by modifying index variable
		for (int i = 1; i <= rs->NumRoots/2; i++) {
			if (error_locator[i]) {
				z = y * i;  // calculate power by multiplying exponents
				z = z + GF2Log(error_locator[i], rs->GF); //multiply by adding exponents
				while (z > (rs->FieldOrder - 2)) {
					z = z - rs->FieldOrder + 1;
				}
				x = x ^ GF2Pow(z, rs->GF);
			}
		}
		x = x ^ error_locator[0];
		if (x == 0) {
			printf("\r\nRS Decode Error Location %li: %li", rs->ErrorCount, j);
			rs->ErrorLocs[rs->ErrorCount] = j;
			rs->ErrorCount++;
            // Todo: check for an ambiguous solution
            
		}
	}

	// Forney algorithm to determine error values
	// first calculate omega, the error magnitude polynomial
	for (int i = 0; i < rs->ErrorCount; i++) {
		correction_poly[i] = syndromes[i];
		for (int j = 1; j <= i; j++) {
			correction_poly[i] = correction_poly[i] ^ GF2Mul(syndromes[i - j], error_locator[j], rs->GF);
		}
	}
	for (int i = 0; i < rs->ErrorCount; i++) { // compute an error value for each error location
		e = block_size - rs->ErrorLocs[i] - 1;
		z = correction_poly[0];
		for (int j = 1; j < rs->ErrorCount; j++) { // calculate numerator
			x = e * j;
			while (x > (rs->FieldOrder - 2)) {
				x = x + 1 - rs->FieldOrder;
			}
			x = rs->FieldOrder - x - 1;
			while (x > (rs->FieldOrder - 2)) {
				x = x + 1 - rs->FieldOrder;
			}
			z = z ^ GF2Mul(correction_poly[j], GF2Pow(x, rs->GF), rs->GF);
		}
		z = GF2Mul(z, GF2Pow(e, rs->GF), rs->GF);
		y = error_locator[1];
		for (int j = 3; j <= rs->NumRoots / 2; j = j + 2) { // operate on each odd coefficient of V
			x = e * (j - 1);
			while (x > (rs->FieldOrder - 2)) {
				x = x + 1 - rs->FieldOrder;
			}
			x = rs->FieldOrder - x - 1;
			while (x > (rs->FieldOrder - 2)) {
				x = x + 1 - rs->FieldOrder;
			}
			y = y ^ GF2Mul(error_locator[j], GF2Pow(x, rs->GF), rs->GF);
		}
		y = GF2Log(y, rs->GF);
		y = rs->FieldOrder - y - 1;
		if (y == (rs->FieldOrder - 1)) {
			y = 0;
		}
		y = GF2Pow(y, rs->GF);
		rs->ErrorMags[i] = GF2Mul(y, z, rs->GF);
		data_block[rs->ErrorLocs[i]] = data_block[rs->ErrorLocs[i]] ^ rs->ErrorMags[i];
	}
	// check for success by calculating syndromes (should be zero if no errors)
	// for (int i = 0; i < rs->NumRoots; i++) {
		// syndromes[i] = 0;
        // x = GF2Pow(rs->FirstRoot + i, rs->GF);
		// for (int j = 0; j < block_size - 1; j++) {
			// syndromes[i] = GF2Mul(syndromes[i] ^ data_block[j], x, rs->GF);
		// }
		// syndromes[i] = syndromes[i] ^ data_block[block_size - 1];
		// if (syndromes[i] != 0) {
			// return -1; // decode fail
		// }
	// }
	int nonzero = calc_syndromes(data_block, block_size, syndromes, rs);
	if (nonzero) {
		return nonzero;
	}
	return rs->ErrorCount; // return number of errors corrected    
}