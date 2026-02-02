#include "gf2.h"

void lfsr_step(GF2_def_struct *gf) {
// utilize Galois configuration to implement LFSR
   unsigned bit0;
   bit0 = gf->LFSR & 1; // save the feedback bit
   gf->LFSR >>= 1; // shift memory register right one bit
   if (bit0) { // if the feedback bit is one
      gf->LFSR ^= (-1) & (gf->GenPoly >> 1); // then XOR the tapped bits
   }
   return;
}

int16_t gf2_get_order(GF2_def_struct *gf) {
    return gf->Order;
}

uint16_t gf2_pow(uint16_t i, GF2_def_struct *gf) {
    return gf->Table[i];
}

uint16_t gf2_log(uint16_t i, GF2_def_struct *gf) {
    return gf->Index[i];
}

uint16_t gf2_inv(uint16_t i, GF2_def_struct *gf) {
    return gf->Inverse[i];
}

int16_t InitGF2(int16_t power, uint32_t genpoly, GF2_def_struct *gf) {
    int32_t i, j;
    gf->Power = power;
    gf->Order = 1;
    gf->GenPoly = genpoly;
    for (i = 0; i < power; i++) {
        gf->Order *= 2;
    }
    // generate the field table and index
    int16_t status = 0;
    gf->LFSR = 1; // start with GF element a^0
    for (i = gf->Order - 2; i >= 0; i--) {
        lfsr_step(gf);
        gf->Table[i] = gf->LFSR;
        gf->Index[gf->LFSR] = i;
        if ((gf->LFSR == 1) && (i > 0)) {
            status++; // number of times sequence repeated during generation
        }
    }  
    gf->Index[0] = 0;
    // generate the inverse table
    gf->Inverse[0] = 0;
    for (i = 1; i < gf->Order; i++) {
        j = 1;
        while (gf2_mul(i, j, gf) != 1) {
            j++;
        }
        gf->Inverse[i] = j;
    }
    return status;
}

uint16_t gf2_mul(uint16_t a, uint16_t b, GF2_def_struct *gf) {
	if ((a == 0) | (b == 0)) {
		return 0;
	}
	a = gf->Index[a];
	b = gf->Index[b];
	a = a + b;
	while (a > (gf->Order - 2)) {
		a = a - (gf->Order - 1);
	}
	return gf->Table[a];
}

uint16_t gf2_div(uint16_t a_arg, uint16_t b_arg, GF2_def_struct *gf) {
	volatile int16_t a, b;
	a = a_arg;
	b = b_arg;
	if (b == 0) {
		return 0xFFFF;
	}
	if (a == 0) {
		return 0;
	}
	a = gf->Index[a];
	b = gf->Index[b];
	a = a - b;
	while (a < 0) {
		a = a + (gf->Order - 1);
	}
	return gf->Table[a];
}

uint16_t gf2_conv(uint16_t *p1, int16_t p1n, uint16_t *p2, int16_t p2n, GF2_def_struct *gf) {
// convolves two gf polynomials
// p1 points to polynomial1 containing p1n uint16_t elements
// p2 points to polynomial2 containing p2n uint16_t elements
// p1 and p2 are stored with lowest power at lowest address
// returns p1n + p2n - 1
// places result in p1
// p1 must have length p1n + p2n - 1
	uint16_t i, j, k;
	k = p1n + p2n - 1;
	uint16_t pr[k];
	for (i = 0; i < k; i++) {
		pr[i] = 0;
	}
	for (i = 0; i < p1n; i++) {
		for (j = 0; j < p2n; j++) {
			pr[i + j] = pr[i + j] ^ gf2_mul(p1[i], p2[j], gf);
		}
	}
	for (i = 0; i < k; i++) {
		p1[i] = pr[i];
	}
	return k;
}
