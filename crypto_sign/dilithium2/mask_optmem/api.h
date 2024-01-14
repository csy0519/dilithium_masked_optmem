#ifndef PQCLEAN_DILITHIUM2_MASK_OPTMEM_API_H
#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_API_H

#include <stddef.h>
#include <stdint.h>

#define NSHARES 2

#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_PUBLICKEYBYTES 1312
#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_SECRETKEYBYTES 2528 + (NSHARES - 1) * (4*96 + 4*96)
#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES 2420
#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_GEN_BUF_BYTES 4376
#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_SIGN_BUF_BYTES 4516
#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_VER_BUF_BYTES 2276
#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_ALGNAME "Dilithium2"


int PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_keypair(uint8_t *pk, uint8_t *sk, uint8_t *buf);

int PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_signature(
    uint8_t *sig, size_t *siglen,
    const uint8_t *m, size_t mlen, const uint8_t *sk, 
    uint8_t *buf);

int PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_verify(
    const uint8_t *sig, size_t siglen,
    const uint8_t *m, size_t mlen, const uint8_t *pk, 
    uint8_t *buf);

int PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign(
    uint8_t *sm, size_t *smlen,
    const uint8_t *m, size_t mlen, const uint8_t *sk, 
    uint8_t *buf);

int PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_open(
    uint8_t *m, size_t *mlen,
    const uint8_t *sm, size_t smlen, const uint8_t *pk, 
    uint8_t *buf);

#endif
