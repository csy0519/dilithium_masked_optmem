#ifndef MASKED_SMALLPOLY_H
#define MASKED_SMALLPOLY_H

#include <stdint.h>
#include "params.h"
#include "reduce.h"
#include "masked_representations.h"
#include "masked.h"

typedef struct{
    int16_t coeffs[N];
} halfpoly;


void PQCLEAN_DILITHIUM2_MASK_OPTMEM_fnt_257(int16_t a[N]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_invfnt_257(int16_t a[N]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_basemul(int16_t *r, int16_t *a, int16_t *b);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_caddqp(int16_t *a);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_csi(int16_t* c, uint32_t* p, uint16_t buff[N + 2*BSSIZE]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_sub_eta(halfpoly *c, const halfpoly *b);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_decompress_chall(halfpoly *c, chall_poly *cp);

#endif