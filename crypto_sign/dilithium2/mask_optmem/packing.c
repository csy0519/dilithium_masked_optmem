#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "masked.h"
#include "masked_poly.h"

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_rho_key(uint8_t *sk, uint8_t *pk, const uint8_t rho[SEEDBYTES], const uint8_t key[SEEDBYTES]) {
    unsigned int i;

    for (i = 0; i < SEEDBYTES; ++i) {
        sk[i] = rho[i];
    }
    sk += SEEDBYTES;

    for (i = 0; i < SEEDBYTES; ++i) {
        sk[i] = key[i];
    }

    for (i = 0; i < SEEDBYTES; ++i) {
        pk[i] = rho[i];
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_s1(uint8_t *sk, const poly *s1, const unsigned int i, uint32_t buf[32]){
    sk += 3 * SEEDBYTES + i * NSHARES * POLYETA_PACKEDBYTES;

    PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_pack(sk, s1, POLYETA_PACKEDBYTES, buf);
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_s2(uint8_t *sk, const poly *s2, const unsigned int i, uint32_t buf[32]){
    sk += 3 * SEEDBYTES + NSHARES * L * POLYETA_PACKEDBYTES + i * NSHARES * POLYETA_PACKEDBYTES;

    PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_pack(sk, s2, POLYETA_PACKEDBYTES, buf);
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_t1(uint8_t *pk, const poly *t1, const unsigned int i){
    pk += SEEDBYTES + i * POLYT1_PACKEDBYTES;

    PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt1_pack(pk, t1);
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_t0(uint8_t *sk, const poly *t0, const unsigned int i){
    sk += 3 * SEEDBYTES + NSHARES * (L + K) * POLYETA_PACKEDBYTES + i * POLYT0_PACKEDBYTES;

    PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt0_pack(sk, t0);
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_tr(uint8_t *sk, const uint8_t tr[SEEDBYTES]){
    unsigned int i;

    sk += 2 * SEEDBYTES;
    
    for (i = 0; i < SEEDBYTES; ++i) {
        sk[i] = tr[i];
    }
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_pk
*
* Description: Unpack public key pk = (rho, t1).
*
* Arguments:   - const uint8_t rho[]: output byte array for rho
*              - const polyveck *t1: pointer to output vector t1
*              - uint8_t pk[]: byte array containing bit-packed pk
**************************************************/
/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_pk(uint8_t rho[SEEDBYTES],
                                        polyveck *t1,
                                        const uint8_t pk[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_PUBLICKEYBYTES]) {
    unsigned int i;

    for (i = 0; i < SEEDBYTES; ++i) {
        rho[i] = pk[i];
    }
    pk += SEEDBYTES;

    for (i = 0; i < K; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt1_unpack(&t1->vec[i], pk + i * POLYT1_PACKEDBYTES);
    }
}*/

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_sk
*
* Description: Unpack secret key sk = (rho, tr, key, t0, s1, s2).
*
* Arguments:   - const uint8_t rho[]: output byte array for rho
*              - const uint8_t tr[]: output byte array for tr
*              - const uint8_t key[]: output byte array for key
*              - const polyveck *t0: pointer to output vector t0
*              - const polyvecl *s1: pointer to output vector s1
*              - const polyveck *s2: pointer to output vector s2
*              - uint8_t sk[]: byte array containing bit-packed sk
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_sk(uint8_t rho[SEEDBYTES],
                                        uint8_t tr[SEEDBYTES],
                                        uint8_t key[SEEDBYTES],
                                        polyveck *t0,
                                        uint32_t BMasks1[NSHARES * COEF_NBITS * N/32 * L],
                                        uint32_t BMasks2[NSHARES * COEF_NBITS * N/32 * K],
                                        const uint8_t sk[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_SECRETKEYBYTES]) {
    unsigned int i;

    for (i = 0; i < SEEDBYTES; ++i) {
        rho[i] = sk[i];
    }
    sk += SEEDBYTES;

    for (i = 0; i < SEEDBYTES; ++i) {
        key[i] = sk[i];
    }
    sk += SEEDBYTES;

    for (i = 0; i < SEEDBYTES; ++i) {
        tr[i] = sk[i];
    }
    sk += SEEDBYTES;

    for (i = 0; i < L; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_unpack(&BMasks1[i * NSHARES * COEF_NBITS * N/32], sk + i * NSHARES * POLYETA_PACKEDBYTES);
    }

    sk += NSHARES * L * POLYETA_PACKEDBYTES;

    for (i = 0; i < K; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_unpack(&BMasks2[i * NSHARES * COEF_NBITS * N/32], sk + i * NSHARES * POLYETA_PACKEDBYTES);
    }

    sk += NSHARES * K * POLYETA_PACKEDBYTES;

    for (i = 0; i < K; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt0_unpack(&t0->vec[i], sk + i * POLYT0_PACKEDBYTES);
    }
}

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_s1(const uint8_t *sk, poly s1[NSHARES], size_t index){
    unsigned int d;
    BitslicePolyEta *BMasks1 = (BitslicePolyEta *) ((uint8_t *) s1 + NSHARES * POLY_BYTES - NSHARES * sizeof(BitslicePolyEta));

    sk += 3 * SEEDBYTES;

    for (d = 0; d < NSHARES; ++d){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_unpack(&BMasks1[d], sk + (d * K + index) * POLYETA_PACKEDBYTES);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_expand(&s1[d], &BMasks1[d]);
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_s2(const uint8_t *sk, poly s2[NSHARES], size_t index){
    unsigned int d;
    BitslicePolyEta *BMasks2 = (BitslicePolyEta *) ((uint8_t *) s2 + NSHARES * POLY_BYTES - NSHARES * sizeof(BitslicePolyEta));

    sk += 3 * SEEDBYTES + (NSHARES * L) * POLYETA_PACKEDBYTES;

    for (d = 0; d < NSHARES; ++d){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_unpack(&BMasks2[d], sk + (d * K + index) * POLYETA_PACKEDBYTES);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_expand(&s2[d], &BMasks2[d]);
    }
}*/

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_sig
*
* Description: Bit-pack signature sig = (c, z, h).
*
* Arguments:   - uint8_t sig[]: output byte array
*              - const uint8_t *c: pointer to PQCLEAN_DILITHIUM2_MASK_OPTMEM_challenge hash length SEEDBYTES
*              - const polyvecl *z: pointer to vector z
*              - const polyveck *h: pointer to hint vector h
**************************************************/
/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_sig(uint8_t sig[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES],
                                       const uint8_t c[SEEDBYTES],
                                       const polyvecl *z,
                                       const polyveck *h) {
    unsigned int i, j, k;

    for (i = 0; i < SEEDBYTES; ++i) {
        sig[i] = c[i];
    }
    sig += SEEDBYTES;

    for (i = 0; i < L; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyz_pack(sig + i * POLYZ_PACKEDBYTES, &z->vec[i]);
    }
    sig += L * POLYZ_PACKEDBYTES;

    Encode h
    for (i = 0; i < OMEGA + K; ++i) {
        sig[i] = 0;
    }

    k = 0;
    for (i = 0; i < K; ++i) {
        for (j = 0; j < N; ++j) {
            if (h->vec[i].coeffs[j] != 0) {
                sig[k++] = (uint8_t) j;
            }
        }

        sig[OMEGA + i] = (uint8_t) k;
    }
}*/

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_hint(uint8_t* sig, const uint8_t *h, size_t n, size_t index){
    unsigned int i, k = 0;

    sig += SEEDBYTES + L * POLYZ_PACKEDBYTES;

    sig[OMEGA + index] = n;
    
    k = (index == 0) ? 0 : sig[OMEGA + index - 1];

    for (i = k; i < n; ++i) {
        sig[i] = h[i - k];
    }

    /For strong unforgeability
    if (index == K - 1) {
        for (i = n; i < OMEGA; ++i){
            sig[i] = 0;
        }
    }
}*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_finish_hint(uint8_t *h, uint32_t n) {
    unsigned int i;
    for (i = n; i < OMEGA; ++i){
        h[i] = 0;
    }
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_sig
*
* Description: Unpack signature sig = (c, z, h).
*
* Arguments:   - uint8_t *c: pointer to output PQCLEAN_DILITHIUM2_MASK_OPTMEM_challenge hash
*              - polyvecl *z: pointer to output vector z
*              - polyveck *h: pointer to output hint vector h
*              - const uint8_t sig[]: byte array containing
*                bit-packed signature
*
* Returns 1 in case of malformed signature; otherwise 0.
**************************************************/
/*int PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_sig(uint8_t c[SEEDBYTES],
                                        polyvecl *z,
                                        polyveck *h,
                                        const uint8_t sig[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES]) {
    unsigned int i, j, k;

    for (i = 0; i < SEEDBYTES; ++i) {
        c[i] = sig[i];
    }
    sig += SEEDBYTES;

    for (i = 0; i < L; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyz_unpack(&z->vec[i], sig + i * POLYZ_PACKEDBYTES);
    }
    sig += L * POLYZ_PACKEDBYTES;

    /Decode h /
    k = 0;
    for (i = 0; i < K; ++i) {
        for (j = 0; j < N; ++j) {
            h->vec[i].coeffs[j] = 0;
        }

        if (sig[OMEGA + i] < k || sig[OMEGA + i] > OMEGA) {
            return 1;
        }

        for (j = k; j < sig[OMEGA + i]; ++j) {
            / Coefficients are ordered for strong unforgeability /
            if (j > k && sig[j] <= sig[j - 1]) {
                return 1;
            }
            h->vec[i].coeffs[sig[j]] = 1;
        }

        k = sig[OMEGA + i];
    }

    / Extra indices are zero for strong unforgeability /
    for (j = k; j < OMEGA; ++j) {
        if (sig[j]) {
            return 1;
        }
    }

    return 0;
}*/

int PQCLEAN_DILITHIUM2_MASK_OPTMEM_verify_hint(const uint8_t sig[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES]) {
    unsigned int i, j, k;

    sig += SEEDBYTES + L * POLYZ_PACKEDBYTES;

    k = 0;
    for (i = 0; i < K; ++i) {
        if (sig[OMEGA + i] < k || sig[OMEGA + i] > OMEGA) {
            return 1;
        }

        for (j = k; j < sig[OMEGA + i]; ++j) {
            /* Coefficients are ordered for strong unforgeability */
            if (j > k && sig[j] <= sig[j - 1]) {
                return 1;
            }
        }

        k = sig[OMEGA + i];
    }

    /* Extra indices are zero for strong unforgeability */
    for (j = k; j < OMEGA; ++j) {
        if (sig[j]) {
            return 1;
        }
    }

    return 0;
}

int PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_hint(uint8_t *h, const uint8_t sig[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES], const size_t i) {
    unsigned int j, k;

    sig += SEEDBYTES + L*POLYZ_PACKEDBYTES;

    k = (i == 0) ? 0 : sig[OMEGA + i - 1];

    for (j = k; j < sig[OMEGA + i]; ++j) {
        h[j - k] = sig[j];
    }

    return sig[OMEGA + i] - k;
}

int PQCLEAN_DILITHIUM2_MASK_OPTMEM_unmask_sk(uint8_t *umsk_sk, const uint8_t *sk){
    unsigned int i;
    uint32_t BMasks1[NSHARES * COEF_NBITS * N/32 * L];
    uint32_t BMasks2[NSHARES * COEF_NBITS * N/32 * K];
    polyvecl s1;
    polyveck s2, t0;
    uint8_t rho[SEEDBYTES], tr[SEEDBYTES], key[SEEDBYTES]; 


    PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_sk(rho, tr, key, &t0, BMasks1, BMasks2, sk);

    for (i = 0; i < L; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_bpoly_unmask(&BMasks1[i * NSHARES * COEF_NBITS * N/32]);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_bitslice2dense(&s1.vec[i], &BMasks1[i * NSHARES * COEF_NBITS * N/32]);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_sub_eta(&s1.vec[i], &s1.vec[i]);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_reduce(&s1.vec[i]);
    }

    for (i = 0; i < K; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_bpoly_unmask(&BMasks2[i * NSHARES * COEF_NBITS * N/32]);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_bitslice2dense(&s2.vec[i], &BMasks2[i * NSHARES * COEF_NBITS * N/32]);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_sub_eta(&s2.vec[i], &s2.vec[i]);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_reduce(&s2.vec[i]);
    }

    for (i = 0; i < SEEDBYTES; ++i) {
        umsk_sk[i] = rho[i];
    }
    umsk_sk += SEEDBYTES;

    for (i = 0; i < SEEDBYTES; ++i) {
        umsk_sk[i] = key[i];
    }
    umsk_sk += SEEDBYTES;

    for (i = 0; i < SEEDBYTES; ++i) {
        umsk_sk[i] = tr[i];
    }
    umsk_sk += SEEDBYTES;

    for (i = 0; i < L; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyeta_pack(umsk_sk + i * POLYETA_PACKEDBYTES, &s1.vec[i]);
    }

    umsk_sk += L * POLYETA_PACKEDBYTES;

    for (i = 0; i < K; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyeta_pack(umsk_sk + i * POLYETA_PACKEDBYTES, &s2.vec[i]);
    }
    umsk_sk += K * POLYETA_PACKEDBYTES;

    for (i = 0; i < K; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt0_pack(umsk_sk + i * POLYT0_PACKEDBYTES, &t0.vec[i]);
    }

    return PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_SECRETKEYBYTES - (NSHARES - 1) * (K + L) * POLYETA_PACKEDBYTES;
}
