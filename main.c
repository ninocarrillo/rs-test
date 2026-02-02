#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "gf2.h"
#include "rs2.h"

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
	printf("\r\nGalois Field Table, order %i:\r\n", gf.Order);
	for (int i = 0; i < (gf.Order-1); i++) {
		printf("%i ", gf.Table[i]);
	}
	
	RS2_def_struct rs;
	rs.GF = &gf;
	InitRS2(rs_first_root, parity_size, &rs);
	
	
	for (int error_count = 1; error_count <= parity_size; error_count++) {
		printf("\r\nError count = %i\r\n", error_count);
		for (int run_number = 1; run_number <= run_count; run_number++) {
			printf(" %i", run_number);
		}
	}
}