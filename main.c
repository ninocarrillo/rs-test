#include <stdio.h>
#include <stdlib.h>
#include "gf2.h"
#include "rs2.h"

#define MAX_BUFFER 2000

void GenRandomMessage(int *buffer, int mask, int size) {
	for (int i = 0; i < size; i++) {
		buffer[i] = rand() & mask;
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

	printf("\n");
	for(int i = 0; i < arg_count; i++) {
		printf("%s ", arg_values[i]);
	}
	printf("\n");
	
	if (arg_count < 6) {
		printf("Not enough arguments.\r\n");
		printf("Usage:\r\nrs-test <gf poly> <message size> <block size> <rs first root> <runs>\r\n");
		return(-1);
	}

	int gf_poly = atoi(arg_values[1]);
	int message_size = atoi(arg_values[2]);
	int block_size = atoi(arg_values[3]);
	int parity_size = block_size - message_size;
	int rs_first_root = atoi(arg_values[4]);
	int run_count = atoi(arg_values[5]);

	printf("\r\nSize of int variable is %li bits.", sizeof(int)*8);
	
	// Initialize Galois Field.
	GF2_def_struct gf;
	int gf_status = InitGF2(gf_poly, &gf);
	printf("\r\nGalois Field Table, order %i:\r\n", gf.Order);
	for (int i = 0; i < (gf.Order-1); i++) {
		printf("%i ", gf.Table[i]);
	}
	if (gf_status > 0) {
		printf("\r\nGalois Field generator polynomial is not irreducible, field repeated %i times.", gf_status);
		return(-1);
	} else if (gf_status < 0) {
		printf("\r\nGalois Field generator polynomial is even, must be odd.");
		return(-1);
	} else {
		printf("\r\nGalois Field generator polynomial is irreducible.");
	}
	printf("\r\nGalois Field order: %i", gf.Order);
	printf("\r\nGalois Field power: %i", gf.Power);
	
	if (block_size > (gf.Order - 1)) {
		printf("\r\nBlock size %i is too large. Must be less than field order %i.\r\n", block_size, gf.Order);
		return(-1);
	}
	
	if (block_size <= message_size) {
		printf("\r\nMessage size %i is too large. Must be less than block size %i.\r\n", message_size, block_size);
		return(-1);
	}


	
	RS2_def_struct rs;
	rs.GF = &gf;
	InitRS2(rs_first_root, parity_size, &rs);

	printf("\r\nRS Genpoly: ");
	for(int i = 0; i < rs.NumRoots + 1; i++){
		printf("%i ", rs.Genpoly[i]);
	}
	printf("\r\n");

	int original_message[MAX_BUFFER];
	int error_vector[MAX_BUFFER];
	int corrupt_message[MAX_BUFFER];
	int decoded_message[MAX_BUFFER];

	int decoder_indicated_failures[parity_size+1];
	int failures[parity_size+1];
	int undetected_failures[parity_size+1];
	int successes[parity_size+1];
	for (int i = 0; i <= parity_size; i++) {
		failures[i] = 0;
		undetected_failures[i] = 0;
		successes[i] = 0;
		decoder_indicated_failures[i] = 0;
	}

	printf("\r\nStarting %i runs.\r\n", (parity_size + 1) * run_count);
	int master_count = 1;
	
	for (int error_count = 0; error_count <= parity_size; error_count++) {
		for (int run_number = 1; run_number <= run_count; run_number++) {
			// printf("\r\n\nError Count %i, Run %i, ", error_count, run_number);
			printf("\r%i", master_count++);
			// Generate a random message to encode.
			GenRandomMessage(original_message, gf.Order - 1, message_size);
			// printf("\r\nMessage:");
			// for (int i = 0; i < message_size; i++) {
				// printf(" %lX", original_message[i]);
			// }
			// Encode message in Reed Solomon block.
			RSEncode(original_message, message_size, &rs);
			// printf("\r\nEncodedMessage:");
			// for (int i = 0; i < block_size; i++) {
			// 	printf(" %X", original_message[i]);
			// }

			GenErrorVector(error_vector, gf.Order - 1, block_size, error_count);
			printf("\r\n             Error Vector:");
			for (int i = 0; i < block_size; i++) {
				printf(" %lX", error_vector[i]);
			}

			CombineVectors(original_message, error_vector, corrupt_message, block_size);
			// printf("\r\nCorrupt Message:");
			// for (int i = 0; i < block_size; i++) {
				// printf(" %lX", corrupt_message[i]);
			// }

			int corrected_count = RSDecode(corrupt_message, block_size, &rs);
			if (corrected_count < 0) {
				decoder_indicated_failures[error_count]++;
			}
			// printf("\r\nCorrected %i errors in message:", corrected_count);
			// for (int i = 0; i < block_size; i++) {
				// printf(" %lX", corrupt_message[i]);
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
					// printf(" %lX,%lX", corrupt_message[i], original_message[i]);
				// }				
			}
			//if (error_count == (block_size-message_size)/2) {
				if (rs.ErrorCount > 0) {
					printf("\r\n          Calculated Error Vector:");
					for (int i = 0; i < block_size; i++) {
						printf(" %i", error_vector[i]);
					}
					printf("\r\n          Detected error locations: ");
					for (int i = 0; i < rs.ErrorCount; i++) {
						printf(" %i", rs.ErrorLocs[i]);
					}
					printf("\r\n          Detected error magnitudes: ");
					for (int i = 0; i < rs.ErrorCount; i++) {
						printf(" %i", rs.ErrorMags[i]);
					}

				}
			//}
		}
	}

	printf("\r\nDecode Success by Error Count:");
	for (int i = 0; i <= parity_size; i++) {
		printf("\r\n%i, %i", i, successes[i]);
	}
	printf("\r\nDecoder Indicated Failures by Error Count:");
	for (int i = 0; i <= parity_size; i++) {
		printf("\r\n%i, %i", i, decoder_indicated_failures[i]);
	}
	printf("\r\nActual Decode Failures by Error Count:");
	for (int i = 0; i <= parity_size; i++) {
		printf("\r\n%i, %i", i, failures[i]);
	}
	printf("\r\nUndetected Decode Failures by Error Count:");
	for (int i = 0; i <= parity_size; i++) {
		printf("\r\n%i, %i", i, undetected_failures[i]);
	}
	printf("\r\nDone.\r\n");
}