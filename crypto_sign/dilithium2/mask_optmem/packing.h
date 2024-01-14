#ifndef PQCLEAN_DILITHIUM2_MASK_OPTMEM_PACKING_H
#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_PACKING_H
#include "params.h"
#include "masked.h"
#include "polyvec.h"
#include "poly.h"
#include <stdint.h>

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_rho_key(uint8_t *sk, uint8_t *pk, const uint8_t rho[SEEDBYTES], const uint8_t key[SEEDBYTES]);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_s1(uint8_t *sk, const poly *s1, const unsigned int i, uint32_t buf[32]);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_s2(uint8_t *sk, const poly *s2, const unsigned int i, uint32_t buf[32]);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_t1(uint8_t *pk, const poly *t1, const unsigned int i);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_t0(uint8_t *sk, const poly *t0, const unsigned int i);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_tr(uint8_t *sk, const uint8_t tr[SEEDBYTES]);

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_sig(uint8_t sig[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES], const uint8_t c[SEEDBYTES], const polyvecl *z, const polyveck *h);*/

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_hint(uint8_t* sig, const uint8_t *h, size_t n, size_t index);*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_finish_hint(uint8_t *h, uint32_t n);

//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_pk(uint8_t rho[SEEDBYTES], polyveck *t1, const uint8_t pk[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_PUBLICKEYBYTES]);


void PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_sk(uint8_t rho[SEEDBYTES],
                                        uint8_t tr[SEEDBYTES],
                                        uint8_t key[SEEDBYTES],
                                        polyveck *t0,
                                        uint32_t BMasks1[NSHARES * COEF_NBITS * N/32 * L],
                                        uint32_t BMasks2[NSHARES * COEF_NBITS * N/32 * K],
                                        const uint8_t sk[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_SECRETKEYBYTES]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_s1(const uint8_t *sk, poly s1[NSHARES], size_t index);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_s2(const uint8_t *sk, poly s2[NSHARES], size_t index);

/*
int PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_sig(uint8_t c[SEEDBYTES], polyvecl *z, polyveck *h, const uint8_t sig[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES]);
*/

int PQCLEAN_DILITHIUM2_MASK_OPTMEM_verify_hint(const uint8_t sig[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES]);

int PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_hint(uint8_t *h, const uint8_t sig[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES], const size_t i);

int PQCLEAN_DILITHIUM2_MASK_OPTMEM_unmask_sk(uint8_t *umsk_sk, const uint8_t *sk);

#endif
