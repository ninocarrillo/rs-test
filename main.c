#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "gf2.h"
#include "rs2.h"

#define MAX_BUFFER 2000

void GenRandomMessage(uint16_t *buffer, int mask, int size) {
	for (int i = 0; i < size; i++) {
		buffer[i] = rand() & mask;
	}
}

void GenErrorVector(uint16_t *buffer, int mask, int size, int count) {
	uint16_t error_locs[MAX_BUFFER];
	// Clear error buffers
	for (int i = 0; i < size; i++) {
		buffer[i] = 0;
		error_locs[i] = 0;
	}
	// Generate count unique error locations in range 0:(size-1)
	int error_index = 0;
	while(error_index < count) {
		int candidate_location = rand() % (size + 1);
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
		buffer[error_locs[i]] = rand() & mask;
	}
}

void CombineVectors(uint16_t *in1, uint16_t *in2, uint16_t *out, int count) {
	for (int i = 0; i < count; i++) {
		out[i] = in1[i] ^ in2[i];
	}
}

int CompareVectors(uint16_t *a, uint16_t *b, int size) {
	int errors = 0;
	for (int i = 0; i < size; i++) {
		if (a[i] != b[i]) {
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
		printf("Usage:\r\nrs-test <message size> <block size> <gf poly> <rs first root> <runs>\r\n");
		return(-1);
	}

	int message_size = atoi(arg_values[1]);
	int block_size = atoi(arg_values[2]);
	int parity_size = block_size - message_size;
	int gf_poly = atoi(arg_values[3]);
	int rs_first_root = atoi(arg_values[4]);
	int run_count = atoi(arg_values[5]);
	
	// Determine the field order based on generator polynomial
	int gf_order = 1;
	int gf_power = 0;
	while ((gf_order<<1) < gf_poly) {
		gf_order <<= 1;
		gf_power++;
	}
	
	if (block_size > (gf_order - 1)) {
		printf("\r\nBlock size %i is too large. Must be less than field order %i.\r\n", block_size, gf_order);
		return(-1);
	}
	
	GF2_def_struct gf;
	InitGF2(gf_power, gf_poly, &gf);
	// printf("\r\nGalois Field Table, order %i:\r\n", gf.Order);
	// for (int i = 0; i < (gf.Order-1); i++) {
	// 	printf("%i ", gf.Table[i]);
	// }
	
	RS2_def_struct rs;
	rs.GF = &gf;
	InitRS2(rs_first_root, parity_size, &rs);

	printf("\r\nRS Genpoly: ");
	for(int i = 0; i < rs.NumRoots + 1; i++){
		printf("%i ", rs.Genpoly[i]);
	}
	printf("\r\n");

	uint16_t original_message[MAX_BUFFER];
	uint16_t error_vector[MAX_BUFFER];
	uint16_t corrupt_message[MAX_BUFFER];
	uint16_t decoded_message[MAX_BUFFER];

	int failures[parity_size];
	int undetected_failures[parity_size];
	for (int i = 0; i < parity_size; i++) {
		failures[i] = 0;
		undetected_failures[i] = 0;
	}


	int master_count = 1;
	
	for (int error_count = 1; error_count <= parity_size; error_count++) {
		for (int run_number = 1; run_number <= run_count; run_number++) {
			printf("\r%i", master_count++);
			// printf("\r\nError Count %i, Run %i, ", error_count, run_number);

			// Generate a random message to encode.
			GenRandomMessage(original_message, 0xFF, message_size);
			// printf("\r\nMessage:");
			// for (int i = 0; i < message_size; i++) {
			// 	printf(" %02X", original_message[i]);
			// }
			// Encode message in Reed Solomon block.
			RSEncode16(original_message, message_size, &rs);
			// printf("\r\nEncodedMessage:");
			// for (int i = 0; i < block_size; i++) {
			// 	printf(" %02X", original_message[i]);
			// }

			GenErrorVector(error_vector, 0xFF, block_size, error_count);
			// printf("\r\nError Vector:");
			// for (int i = 0; i < block_size; i++) {
			// 	printf(" %02X", error_vector[i]);
			// }

			CombineVectors(original_message, error_vector, corrupt_message, block_size);
			// printf("\r\nCorrupt Message:");
			// for (int i = 0; i < block_size; i++) {
			// 	printf(" %02X", corrupt_message[i]);
			// }

			int corrected_count = RSDecode16(corrupt_message, block_size, &rs);
			// printf("\r\nCorrected %i errors in message:", corrected_count);
			// for (int i = 0; i < block_size; i++) {
			// 	printf(" %02X", corrupt_message[i]);
			// }

			int errors = CompareVectors(corrupt_message, original_message, message_size);
			if (errors > 0) {
				failures[error_count - 1]++;
				if (corrected_count >= 0) {
					undetected_failures[error_count - 1]++;
				}
			}
		}
	}

	printf("\r\nDecode Failures by Error Count:");
	for (int i = 0; i < parity_size; i++) {
		printf("\r\n%i, %i", i+1, failures[i]);
	}
	printf("\r\nUndetected Decode Failures by Error Count:");
	for (int i = 0; i < parity_size; i++) {
		printf("\r\n%i, %i", i+1, undetected_failures[i]);
	}
	printf("\r\nDone.\r\n");
}