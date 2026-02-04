#include "rs2.h"
#include "gf2.h"

void InitRS2(int first_root, int num_roots, RS2_def_struct *rs) {
    rs->FirstRoot = first_root;
    rs->NumRoots = num_roots;
    rs->FieldOrder = gf2_get_order(rs->GF);
    // Generate Reed Solomon generator polynomial through convolution of polynomials.
    // rs->genpoly = (x + a^b)(x + a^b+1)...(x + a^b+r-1)
    // start with rs->genpoly = x + a^b
    // lowest order coefficient in lowest index of array
    rs->Genpoly[0] = gf2_pow(first_root, rs->GF);
    rs->Genpoly[1] = 1;
    
    int factorpoly[2];
    // preload the x^1 coefficient in the factor polynomial
    factorpoly[1] = 1;
    int i;
    for (i = first_root + 1; i < first_root + num_roots; i++) {
        factorpoly[0] = gf2_pow(i, rs->GF);
        gf2_conv(&rs->Genpoly[0], i + 1 - first_root, &factorpoly[0], 2, rs->GF);
    }
    rs->MinimumErrorDistance = 0;
}

void RSEncode(int *message, int block_size, RS2_def_struct *rs) {
	int quotient[MAX_GENPOLY_ROOTS + 1];
	for (int i = 0; i < rs->NumRoots; i++) {
		message[i + block_size] = 0;
	}
	for (int i = 0; i < rs->NumRoots + 1; i++) {
		quotient[i] = message[i];
	}
	for (int i = 0; i < block_size; i++) {
		int x = quotient[0];
		for (int j = 1; j < rs->NumRoots + 1; j++) {
			quotient[j - 1] = gf2_mul(x, rs->Genpoly[rs->NumRoots - j], rs->GF) ^ quotient[j];
		}
		quotient[rs->NumRoots] = message[i + rs->NumRoots + 1];
	}
	int i = block_size;
	for (int j = 0; j < rs->NumRoots; j++) {
		message[i] = quotient[j];
		i++;
	}
}

int RSDecode(int *data_block, int block_size, RS2_def_struct *rs) {
	int syndromes[MAX_GENPOLY_ROOTS];
    int error_locator[MAX_GENPOLY_ROOTS];
    int correction_poly[MAX_GENPOLY_ROOTS + 1];
    int next_error_locator[MAX_GENPOLY_ROOTS];
	int step_factor, order_tracker;
    int i, j;
	int x, y, z, e;
    
// find the syndromes, S
// calculate one syndrome for each root of rs->genpoly
	for (i = 0; i < rs->NumRoots; i++) {
		syndromes[i] = 0;
        x = gf2_pow(rs->FirstRoot + i, rs->GF);
		y = block_size - 1;
		for (j = 0; j < y; j++) {
			syndromes[i] = syndromes[i] ^ data_block[j];
			syndromes[i] = gf2_mul(syndromes[i], x, rs->GF);
		}
		syndromes[i] = syndromes[i] ^ data_block[j];
	}

    // Berlekamp's Algorothm
    // calculate the error locator
    // uses:
    // correction polynomial C
    // step parameter K
    // order tracker L
    // error value x
    
    // clear out arrays
	for (i = 0; i <  rs->NumRoots; i++) {
		error_locator[i] = 0;
		correction_poly[i] = 0;
		next_error_locator[i] = 0;
	}
	error_locator[0] = 1;
	correction_poly[1] = 1;	// C = x
	order_tracker = 0;
	for (step_factor = 1; step_factor <= rs->NumRoots; step_factor++) {
        // first calculate error value, e
		y = step_factor - 1;
		e = syndromes[y];
		for (i = 1; i <= order_tracker; i++) {
			x = y - i;
			e = e ^ gf2_mul(error_locator[i], syndromes[x], rs->GF);
		}
        // now update the estimate of V[x]
		if (e != 0) {
			for (i = 0; i <= order_tracker; i++) {
				next_error_locator[i] = error_locator[i] ^ gf2_mul(e, correction_poly[i], rs->GF);
			}
            // and update the value of C
			e = gf2_inv(e, rs->GF);
			for (i = 0; i <= rs->NumRoots / 2; i++) {		// need to refine loop limit t, could be based on L?
				correction_poly[i] = gf2_mul(error_locator[i], e, rs->GF);
			}
			for (i = 0; i <= rs->NumRoots / 2; i++) {	// refine loop limit
				error_locator[i] = next_error_locator[i];         // this should only happen if e != 0!
			}
		}
		if ((2 * order_tracker) < step_factor) {
			order_tracker = step_factor - order_tracker;
		}
        // multiply C(x) by x (increase power by one)
		for (i = rs->NumRoots; i > 0; i--) {
			correction_poly[i] = correction_poly[i - 1];
		}
		correction_poly[0] = 0;
	}
    // now solve the error locator polynomial to find the error positions
    // by using the Chien Search
    // Uses V[i] calculated above
	rs->ErrorCount = 0;  // rs->error_count is the number of errors found
	for (j = 0; j < block_size; j++) {
		x = 0;
		y = j + rs->FieldOrder - block_size;  // account for code shortening by modifying index variable
		for (i = 1; i <= rs->NumRoots / 2; i++) {
			if (error_locator[i]) {
				z = y * i;  // calculate power by multiplying exponents
				z = z + gf2_log(error_locator[i], rs->GF); //multiply by adding exponents
				while (z > (rs->FieldOrder - 2)) {
					z = z - rs->FieldOrder + 1;
				}
				x = x ^ gf2_pow(z, rs->GF);
			}
		}
		x = x ^ error_locator[0];
		if (x == 0) {
			rs->ErrorLocs[rs->ErrorCount] = j;
			rs->ErrorCount++;
            // Now check for an ambiguous solution
            
		}
	}
		if(1) {
        //if (rs->ErrorCount <= ((rs->NumRoots / 2) - rs->MinimumErrorDistance)) {
        // Forney algorithm to determine error values
        // first calculate omega, the error magnitude polynomial
        for (i = 0; i < rs->ErrorCount; i++) {
            correction_poly[i] = syndromes[/*rs->FirstRoot + */i];
            for (j = 1; j <= i; j++) {
                correction_poly[i] = correction_poly[i] ^ gf2_mul(syndromes[/*rs->FirstRoot + */i - j], error_locator[j], rs->GF);
            }
        }
        for (i = 0; i < rs->ErrorCount; i++) { // compute an error value for each error location
            e = block_size - rs->ErrorLocs[i] - 1;
            z = correction_poly[0];
            for (j = 1; j < rs->ErrorCount; j++) { // calculate numerator
                x = e * j;
                while (x > (rs->FieldOrder - 2)) {
                    x = x + 1 - rs->FieldOrder;
                }
                x = rs->FieldOrder - x - 1;
                while (x > (rs->FieldOrder - 2)) {
                    x = x + 1 - rs->FieldOrder;
                }
                z = z ^ gf2_mul(correction_poly[j], gf2_pow(x, rs->GF), rs->GF);
            }
            z = gf2_mul(z, gf2_pow(e, rs->GF), rs->GF);
            y = error_locator[1];
            for (j = 3; j <= rs->NumRoots / 2; j = j + 2) { // operate on each odd coefficient of V
                x = e * (j - 1);
                while (x > (rs->FieldOrder - 2)) {
                    x = x + 1 - rs->FieldOrder;
                }
                x = rs->FieldOrder - x - 1;
                while (x > (rs->FieldOrder - 2)) {
                    x = x + 1 - rs->FieldOrder;
                }
                y = y ^ gf2_mul(error_locator[j], gf2_pow(x, rs->GF), rs->GF);
            }
            y = gf2_log(y, rs->GF);
            y = rs->FieldOrder - y - 1;
            if (y == (rs->FieldOrder - 1)) {
                y = 0;
            }
            y = gf2_pow(y, rs->GF);
            rs->ErrorMags[i] = gf2_mul(y, z, rs->GF);
            data_block[rs->ErrorLocs[i]] = data_block[rs->ErrorLocs[i]] ^ rs->ErrorMags[i];
        }
        // check for success by calculating syndromes (should be zero if no errors)
        for (i = 0; i < rs->NumRoots; i++) {
            syndromes[i] = 0;
            x = gf2_pow(rs->FirstRoot + i, rs->GF);
            y = block_size - 1;
            for (j = 0; j < y; j++) {
                syndromes[i] = syndromes[i] ^ data_block[j];
                syndromes[i] = gf2_mul(syndromes[i], x, rs->GF);
            }
            syndromes[i] = syndromes[i] ^ data_block[j];
            if (syndromes[i] != 0) {
                return -1; // decode fail
            }
        }
    } else {
        rs->ErrorCount = -2;
    }
	return rs->ErrorCount; // return number of errors corrected    
}