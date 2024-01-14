#ifndef MASKEDPOLY_H
#define MASKEDPOLY_H

#include "poly.h"
#include "polyvec.h"
#include "gadgets.h"
#include "masked.h"
#include "masked_representations.h"
#include "params.h"

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyveck_unmask(polyveck v[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecl_unmask(polyvecl v[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_unmask(poly a[NSHARES]);*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_bpoly_unmask(uint32_t a[NSHARES * COEF_NBITS * N/32]);

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_negate_unmask(poly a[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecl_ntt(polyvecl v[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyveck_ntt(polyveck v[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_ntt(poly a[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyveck_invntt_tomont(polyveck v[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecl_invntt_tomont(polyvecl v[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_invntt_tomont(poly a[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyveck_pointwise_poly_montgomery(polyveck r[NSHARES], const poly *a, const polyveck v[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecl_pointwise_poly_montgomery(polyvecl r[NSHARES], const poly *a, const polyvecl v[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_pointwise_poly_montgomery(poly r[NSHARES], const poly *a, const poly v[NSHARES]);*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_pack(uint8_t *b, const poly *a, const size_t mask_stride, uint32_t buf[32]); 

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_unpack(uint32_t *r, const uint8_t *p);

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecleta_unmask(polyvecl *r, const BitslicePolyveclEta *p, uint32_t buf[32 + POLYETA_PACKEDBYTES/4]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecketa_unmask(polyveck *r, const BitslicePolyveckEta *p, uint32_t buf[32 + POLYETA_PACKEDBYTES/4]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_unmask(poly *r, const BitslicePolyEta *p, const size_t mask_stride, uint32_t buf[32 + POLYETA_PACKEDBYTES/4]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecketa_expand(polyveck r[NSHARES], const BitslicePolyveckEta v[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecleta_expand(polyvecl r[NSHARES], const BitslicePolyveclEta v[NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_expand(poly *r, const BitslicePolyEta *a);*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecketa_secb2a_modp(polyveck r[NSHARES],
        uint32_t buf[2 * BSSIZE * NSHARES + 2 * COEF_NBITS * NSHARES + 2 * (COEF_NBITS + 1) * NSHARES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecleta_secb2a_modp(polyvecl r[NSHARES],
        uint32_t buf[2 * BSSIZE * NSHARES + 2 * COEF_NBITS * NSHARES + 2 * (COEF_NBITS + 1) * NSHARES]);

/*
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_bpoly_secb2a_modp(uint32_t *a, size_t mask_stride, 
        uint32_t buf[2 * BSSIZE * NSHARES + 2 * COEF_NBITS * NSHARES + 2 * (COEF_NBITS + 1) * NSHARES]);
*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_bpoly_secb2a_mod257(uint32_t *a, size_t mask_stride, 
        uint32_t buf[COEF_257_NBITS * (NSHARES - 1) + (COEF_257_NBITS + 1) * NSHARES + 4 * NSHARES + 4]);

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_bitslice2dense(poly *v, uint32_t *u);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyveck_negate(polyveck v[NSHARES]);*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_stream_add_uniform_gamma1(poly *r,
        uint32_t a[NSHARES * COEF_257_NBITS * N/32],
        const uint8_t seed[CRHBYTES],
        uint16_t nonce,
        uint8_t buf[sizeof(stream256_state) + 16 + 4*((2*COEF_NBITS + 5)*NSHARES - COEF_NBITS)]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_stream_sub_from_w0(poly *r,
        uint32_t a[NSHARES * COEF_257_NBITS * N/32], 
        comp_w0_poly *w0,
        uint8_t buf[4*((2*COEF_NBITS + 5)*NSHARES - COEF_NBITS)]);

#endif