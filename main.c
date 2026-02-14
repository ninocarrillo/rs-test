#include <stdio.h>
#include <stdlib.h>
#include "gf2.h"
#include "rs2.h"

#define MAX_BUFFER MAX_FIELD_SIZE

void GenRandomMessage(int *buffer, int mask, int size) {
	for (int i = 0; i < size; i++) {
		buffer[i] = rand() & mask;
	}
}
void CopyMessage(int *in, int *out, int size) {
	for (int i = 0; i < size; i++) {
		out[i] = in[i];
	}
}

void GenErrorVector(int *buffer, int mask, int size, int count) {
	int error_locs[MAX_BUFFER];
	// Clear error buffers
	for (int i = 0; i < size; i++) {
		buffer[i] = 0;
		error_locs[i] = 0;
	}
	// Generate count unique error locations in range 0:(size-1)
	int error_index = 0;
	while(error_index < count) {
		int candidate_location = rand() % size;
		int is_unique = 1;
		int i = 0;
		while (is_unique && (i < error_index)) {
			if (error_locs[i++] == candidate_location) {
				is_unique = 0;
			}
		}
		if (is_unique) {
			error_locs[error_index++] = candidate_location;
		}
	}
	for (int i = 0; i < count; i++) {
		int x = 0;
		while (x == 0) {
			x = rand() & mask;
		}
		buffer[error_locs[i]] = x;
	}
}

void CombineVectors(int *in1, int *in2, int *out, int count) {
	for (int i = 0; i < count; i++) {
		out[i] = in1[i] ^ in2[i];
	}
}

int CompareVectors(int *a, int *b, int size) {
	int errors = 0;
	for (int i = 0; i < size; i++) {
		if (a[i] ^ b[i]) {
			errors++;
		}
	}
	return errors;
}

int main(int arg_count, char* arg_values[]) {
	
	if (arg_count < 8) {
		printf("Not enough arguments.\r\n");
		printf("Usage:\r\nrs-test <gf poly> <rs first root> <block size> <message size> <max error count> <runs> <seed>\r\n");
		printf("\r\nExample: rs-test 285 0 15 13 7 100000 0");
		printf("\r\n\n     gf poly:");
		printf("\r\n              Integer number representing the Galois Field reducing polynomial, in GF(2).");
		printf("\r\n              The specified polynomial also defines the Galois Field element size. ");
		printf("\r\n              For example, x^8+x^4+x^3+x^2+1 is represented as 285, and valid for GF(2^8).");
		printf("\r\n\n     rs first root:");
		printf("\r\n              Integer number, typically 0 or 1, sometimes called First Consecutive Root");
		printf("\r\n              or FCR in literature.");
		printf("\r\n\n     block size:");
		printf("\r\n              Integer number of elements in the data block, including parity. Sometimes called");
		printf("\r\n              'n' in literature.");
		printf("\r\n\n     message size:");
		printf("\r\n              Integer number of payload message elements in the data block. Sometimes called");
		printf("\r\n              'k' in literature. The number of parity symbols computed per block is n-k.");
		printf("\r\n       max error count:");
		printf("\r\n               Maximum number of errors per block.");
		printf("\r\n\n     runs:");
		printf("\r\n              Integer number of random test cases to perform at each error count. The program");
		printf("\r\n              will generate a random message of specified length for each run, and corrupt");
		printf("\r\n              the message with a precise number of random errors in random locations.");
		printf("\r\n              Error count will span from zero to (n-k).");
		printf("\r\n\n     seed:");
		printf("\r\n              Integer number used to seed random number generator, for test repeatability.");
		printf("\r\n");

		return(-1);
	}

	int gf_poly = atoi(arg_values[1]);
	int rs_first_root = atoi(arg_values[2]);
	int block_size = atoi(arg_values[3]);
	int message_size = atoi(arg_values[4]);
	int max_errors = atoi(arg_values[5]);
	int run_count = atoi(arg_values[6]);
	int seed = atoi(arg_values[7]);
	int parity_size = block_size - message_size;
	srand(seed);
	
	// Initialize Galois Field.
	GF2_def_struct gf;
	int gf_status = InitGF2(gf_poly, &gf);
	if (gf_status > 0) {
		printf("\r\nGalois Field generator polynomial %i is not irreducible, field repeated %i times.\r\n", gf_poly, gf_status);
		return(-1);
	} else if (gf_status < 0) {
		printf("\r\nGalois Field generator polynomial %i is even, must be odd.\r\n", gf_poly);
		return(-1);
	} else {
	}
	
	if (block_size > (gf.Order - 1)) {
		printf("\r\nBlock size %i is too large. Must be less than field order %i.\r\n", block_size, gf.Order);
		return(-1);
	}
	
	if (block_size <= message_size) {
		printf("\r\nMessage size %i is too large. Must be less than block size %i.\r\n", message_size, block_size);
		return(-1);
	}

	if (max_errors >= block_size) {
		printf("\r\nMax error count %i is too large. Must be less than block size %i.\r\n", max_errors, block_size);
		return(-1);
	}
	if (max_errors < 1) {
		printf("\r\nMax error count %i is too small. Must be greater than zero.\r\n", max_errors);
		return(-1);
	}
	

	printf("\r\nGalois Field generator polynomial %i is irreducible.", gf_poly);
	printf("\r\nGalois Field contains %i elements.", gf.Order);
	printf("\r\nGalois Field element size is %i bits.", gf.Power);

	printf("\r\nGalois Field Table:");
	for (int i = 0; i < gf.Order; i++) {
		if ((i % 16) == 0) {
			printf("\r\n");
		}
		if (i == 0) {
			printf("%5i", 0);
		} else {
			printf("%5i", gf.Table[i-1]);
		}
	}


	printf("\r\nSize of int variable is %li bits.", sizeof(int)*8);
	
	RS2_def_struct rs;
	rs.GF = &gf;
	InitRS2(rs_first_root, parity_size, &rs);

	printf("\r\nReed Solomon Generator Polynomial, highest coefficient first:\r\n");
	for(int i = 0; i < rs.NumRoots + 1; i++){
		printf("%i ", rs.Genpoly[i]);
	}
	printf("\r\n");

	int original_message[MAX_BUFFER];
	int error_vector[MAX_BUFFER];
	int corrupt_message[MAX_BUFFER];
	int decoded_message[MAX_BUFFER];
	int reencoded_message[MAX_BUFFER];

	int decoder_indicated_failures[MAX_FIELD_SIZE];
	int failures[MAX_FIELD_SIZE];
	int undetected_failures[MAX_FIELD_SIZE];
	int successes[MAX_FIELD_SIZE];
	int artificial_codewords[MAX_FIELD_SIZE];
	for (int i = 0; i <= MAX_FIELD_SIZE; i++) {
		failures[i] = 0;
		undetected_failures[i] = 0;
		successes[i] = 0;
		decoder_indicated_failures[i] = 0;
		artificial_codewords[i] = 0;
	}

	printf("\r\nStarting %i runs.\r\n", (max_errors + 1) * run_count);
	int master_count = 1;
	
	for (int error_count = 0; error_count <= max_errors; error_count++) {
		for (int run_number = 1; run_number <= run_count; run_number++) {
			// printf("\r\n\nError Count %i, Run %i, ", error_count, run_number);
			printf("\r%i", master_count++);
			// Generate a random message to encode.
			GenRandomMessage(original_message, gf.Order - 1, message_size);
			// printf("\r\nMessage:");
			// for (int i = 0; i < message_size; i++) {
				// printf(" %X", original_message[i]);
			// }
			// Encode message in Reed Solomon block.
			RSEncode(original_message, message_size, &rs);
			// printf("\r\nEncodedMessage:");
			// for (int i = 0; i < block_size; i++) {
			// 	printf(" %X", original_message[i]);
			// }

			GenErrorVector(error_vector, gf.Order - 1, block_size, error_count);
			// printf("\r\n             Error Vector:");
			// for (int i = 0; i < block_size; i++) {
				// printf(" %i", error_vector[i]);
			// }

			CombineVectors(original_message, error_vector, corrupt_message, block_size);
			// printf("\r\nCorrupt Message:");
			// for (int i = 0; i < block_size; i++) {
				// printf(" %X", corrupt_message[i]);
			// }

			CopyMessage(corrupt_message, reencoded_message, message_size);
			RSEncode(reencoded_message, message_size, &rs);
			// Check if the randomly corrupted message is also a valid codeword
			if ((CompareVectors(corrupt_message, reencoded_message, block_size) == 0) && (error_count > 1)) {
				artificial_codewords[error_count]++;
			}

			int corrected_count = RSDecode(corrupt_message, block_size, &rs);
			if (corrected_count < 0) {
				decoder_indicated_failures[error_count]++;
			}
			// printf("\r\nCorrected %i errors in message:", corrected_count);
			// for (int i = 0; i < block_size; i++) {
				// printf(" %X", corrupt_message[i]);
			// }

			int errors = CompareVectors(corrupt_message, original_message, block_size);
			// printf("\r\nBlock size: %i, Errors: %i", block_size, errors);
			if (errors > 0) {
				failures[error_count]++;
				if (corrected_count >= 0) {
					undetected_failures[error_count]++;
				}
			} else {
				successes[error_count]++;
				// printf("\r\nSuccessful message:", corrected_count);
				// for (int i = 0; i < block_size; i++) {
					// printf(" %X,%X", corrupt_message[i], original_message[i]);
				// }				
			}
			if ((errors > 0) && (error_count <= parity_size/2)) {
			//if (corrected_count > 0) {
				
				printf("\r\n          Original Message, Encoded:");
				for (int i = 0; i < block_size; i++) {
					printf(" %i", original_message[i]);
				}
				printf("\r\n          Actual Error Vector:");
				for (int i = 0; i < block_size; i++) {
					printf(" %i", error_vector[i]);
				}
				printf("\r\n          Corrupt Message:");
				for (int i = 0; i < block_size; i++) {
					printf(" %i", corrupt_message[i]);
				}
				printf("\r\n          Syndromes:");
				for (int i = 0; i < rs.NumRoots; i++) {
					printf(" %i", rs.SavedSyndromes[i]);
				}
				printf("\r\n          Detected error indices: ");
				for (int i = 0; i < rs.ErrorCount; i++) {
					printf(" %i", rs.ErrorIndices[i]);
				}
				printf("\r\n          Detected error roots: ");
				for (int i = 0; i < rs.ErrorCount; i++) {
					printf(" %i", rs.ErrorLocatorRoots[i]);
				}
				printf("\r\n          Detected error magnitudes: ");
				for (int i = 0; i < rs.ErrorCount; i++) {
					printf(" %i", rs.ErrorMags[i]);
				}
				printf("\r\n          Error Locator Poly:");
				for (int i = 0; i <= rs.NumRoots/2; i++) {
					printf(" %i", rs.ErrorLocatorPoly[i]);
				}
				printf("\r\n          Error Magnitude Poly:");
				for (int i = 0; i <= rs.NumRoots/2; i++) {
					printf(" %i", rs.ErrorMagPoly[i]);
				}
				// printf("\r\n          RS Gen Poly:");
				// for (int i = 0; i < rs.NumRoots+1; i++) {
				// 	printf(" %i", rs.Genpoly[i]);
				// }
			}
		}
	}

	printf("\r\nDecode Success by Error Count:");
	for (int i = 0; i <= max_errors; i++) {
		printf("\r\n%i, %i", i, successes[i]);
	}
	printf("\r\nDecoder Indicated Failures by Error Count:");
	for (int i = 0; i <= max_errors; i++) {
		printf("\r\n%i, %i", i, decoder_indicated_failures[i]);
	}
	printf("\r\nActual Decode Failures by Error Count:");
	for (int i = 0; i <= max_errors; i++) {
		printf("\r\n%i, %i", i, failures[i]);
	}
	printf("\r\nUndetected Decode Failures by Error Count:");
	for (int i = 0; i <= max_errors; i++) {
		printf("\r\n%i, %i", i, undetected_failures[i]);
	}
	printf("\r\nArtificial Codewords Generated by Error Count:");
	for (int i = 0; i <= max_errors; i++) {
		printf("\r\n%i, %i", i, artificial_codewords[i]);
	}
	printf("\r\nDone.\r\n");
}