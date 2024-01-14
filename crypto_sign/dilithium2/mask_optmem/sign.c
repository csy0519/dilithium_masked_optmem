#include "shake.h"
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"
#include "sign.h"
#include "symmetric.h"
#include "masked.h"
#include "masked_poly.h"
#include "masked_representations.h"
#include "masked_smallpoly.h"
#include <stdint.h>

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_keypair
*
* Description: Generates public and private key.
*
* Arguments:   - uint8_t *pk: pointer to output public key (allocated
*                             array of PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key (allocated
*                             array of PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_SECRETKEYBYTES bytes)
*              - uint8_t *buf address of the working memory. 
*                           
*
* Returns 0 (success)
**************************************************/
int PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_keypair(uint8_t *pk, uint8_t *sk, uint8_t *buf) {
    unsigned int i;
    uint8_t *seedbuf, *genbuf;
    uint8_t *tr, *rho, *rhoprime, *key;
    poly *s1, *s2, *that, *t1, *t0;
    comp_polyveck *t;

    /* Memory repartition */
    rhoprime = buf;
    t = (comp_polyveck *) (rhoprime + CRHBYTES);
    s1 = (poly *) ((uint8_t *) t + K * sizeof(comp_poly));
    genbuf = (uint8_t *) s1 + sizeof(poly);

    /* Aliases */
    seedbuf = rhoprime;
    rho = seedbuf;
    key = seedbuf + SEEDBYTES + CRHBYTES;
    tr = rhoprime;

    s2 = s1;
    that = s1;
    t1 = s1;
    t0 = s1;

    /* Get randomness for rho, rhoprime and key */
    randombytes(seedbuf, SEEDBYTES);
    shake256(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES, genbuf);
    PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_rho_key(sk, pk, rho, key);

    /* Rho is now stored in pk */
    rho = pk;

    /* Copy rhoprime to another location*/
    for (i = 0; i < CRHBYTES; ++i){
        rhoprime[i] = seedbuf[i + SEEDBYTES];
    }

    PQCLEAN_DILITHIUM2_MASK_OPTMEM_cpolyveck_init(t, 0);

    /* Computes Ahat * s1hat and saves s1 */
    for (i = 0; i < L; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_uniform_eta(s1, rhoprime, i, genbuf);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_s1(sk, s1, i, (uint32_t*) genbuf);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_ntt(s1);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyvec_stream_acc_matrix_pointwise(t, s1, rho, i, genbuf);
    }

    /* Computes t = A * s1 + s2, t1, t0 and generates s2*/
    for (i = 0; i < K; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompress(that, &t->vec[i]);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_invntt_tomont(that);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress(&t->vec[i], that);

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_uniform_eta(s2, rhoprime, L + i, genbuf);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_s2(sk, s2, i, (uint32_t*) genbuf);

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_add_cpoly(t1, s2, &t->vec[i]);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_caddq(t1);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress(&t->vec[i], t1);
        
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_power2round_a1(t1, t1);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_t1(pk, t1, i);
        
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_power2round_a0(t0, t1, &t->vec[i]);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_t0(sk, t0, i);
    }

    /* Compute H(rho, t1) and write secret key */
    shake256(tr, SEEDBYTES, pk, PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_PUBLICKEYBYTES, genbuf);
    PQCLEAN_DILITHIUM2_MASK_OPTMEM_pack_tr(sk, tr);

    return 0;
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_signature
*
* Description: Computes signature.
*
* Arguments:   - uint8_t *sig:   pointer to output signature (of length PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES)
*              - size_t *siglen: pointer to output length of signature
*              - uint8_t *m:     pointer to message to be signed
*              - size_t mlen:    length of message
*              - uint8_t *sk:    pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_signature(uint8_t *sig,
        size_t *siglen,
        const uint8_t *m,
        size_t mlen,
        const uint8_t *sk,
        uint8_t *buf) {
    unsigned int i;
    uint32_t n;
    uint8_t *hash_buf;
    uint32_t *mult_buf;
    const uint8_t *rho, *tr, *key, *sk_s1, *sk_s2, *sk_t0;
    uint8_t *mu, *rhoprime, *h, *cbuff;
    uint16_t nonce = 0;
    uint32_t *b_s1, *b_s2;
    poly *y, *z, *w1, *t2; 
    halfpoly *c;
    chall_poly *cp;
    shake256incctx *state;
    comp_polyveck *w;
    comp_w0_polyveck *w0;
    comp_poly *dc;

    /* Memory repartition */
    mu = buf;
    rhoprime = mu + CRHBYTES;
    w = (comp_polyveck *) (rhoprime + CRHBYTES);
    c = (halfpoly *) ((uint8_t *) w + sizeof(comp_w0_polyveck));
    b_s1 = (uint32_t *) ((uint8_t *) c + sizeof(halfpoly));
    mult_buf = (uint32_t *) ((uint8_t *)b_s1 + NSHARES * N/32 * COEF_257_NBITS * 4);
    cp = (chall_poly *) ((uint8_t *) mult_buf + sizeof(shake256incctx) + 16 + 4*51*NSHARES - 23*4 + 4 + 128);
    hash_buf = (uint8_t *) w + sizeof(poly) + sizeof(comp_polyveck);

    /* Aliases */
    z = (poly *) c;
    b_s2 = b_s1;
    w1 = (poly *)((uint8_t *) w + sizeof(comp_polyveck));
    t2 = z;
    y = w1;
    dc = (comp_poly *) ((uint8_t *) t2 + sizeof(poly));
    state = (shake256incctx *) hash_buf;
    cbuff = (uint8_t *) w1;
    w0 = (comp_w0_polyveck *) w;

    /* "Unpacking" of sk */ 
    rho = sk;
    key = sk + SEEDBYTES;
    tr = sk + 2*SEEDBYTES;
    sk_s1 = sk + 3*SEEDBYTES;
    sk_s2 = sk + 3*SEEDBYTES + NSHARES * L * POLYETA_PACKEDBYTES;
    sk_t0 = sk + 3*SEEDBYTES + NSHARES * (L + K) * POLYETA_PACKEDBYTES;

    h = sig + SEEDBYTES + L * POLYZ_PACKEDBYTES;

    /* Compute CRH(tr, msg) */
    shake256_inc_init(state);
    shake256_inc_absorb(state, tr, SEEDBYTES);
    shake256_inc_absorb(state, m, mlen);
    shake256_inc_finalize(state);
    shake256_inc_squeeze(mu, CRHBYTES, state);

    shake256_inc_init(state);
    shake256_inc_absorb(state, key, SEEDBYTES);
    shake256_inc_absorb(state, mu, CRHBYTES);
    shake256_inc_finalize(state);
    shake256_inc_squeeze(rhoprime, CRHBYTES, state);

    nonce = (uint16_t) -L;
rej:
    nonce += L;
    PQCLEAN_DILITHIUM2_MASK_OPTMEM_cpolyveck_init(w, 0);

    /* Compute what and store it in compressed form */
    for (i = 0; i < L; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_uniform_gamma1(y, rhoprime, (uint16_t) (nonce + i), hash_buf);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_ntt(y);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_reduce(y);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyvec_stream_acc_matrix_pointwise(w, y, rho, i, hash_buf);
    }

    /* Compute w = NTT(what) and w1 */
    for (i = 0; i < K; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompress(w1, &w->vec[i]);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_invntt_tomont(w1);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_caddq(w1);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress_w_pack_w1(&w0->vec[i], sig + i * POLYW1_PACKEDBYTES, w1);
        /*PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress(&w->vec[i], w1);

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompose_a1(w1, w1);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyw1_pack(sig + i * POLYW1_PACKEDBYTES, w1);*/
    }

    /* Compute CRH(mu, w1) and store it in the signature */
    shake256_inc_init(state);
    shake256_inc_absorb(state, mu, CRHBYTES);
    shake256_inc_absorb(state, sig, K * POLYW1_PACKEDBYTES);
    shake256_inc_finalize(state);
    shake256_inc_squeeze(sig, SEEDBYTES, state);

    /* Compute the challenge */
    PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compact_challenge(cp, sig, cbuff, hash_buf);

    /* Compute z = y + cs1*/
    for (i = 0; i < L; ++i) {
        /* Unpack s1_i and convert from binary masking to arithmetic masking */
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_unpack(b_s1, sk_s1 + i * NSHARES * POLYETA_PACKEDBYTES);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_bpoly_secb2a_mod257(b_s1, COEF_257_NBITS, (uint32_t *) mult_buf);

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_decompress_chall(c, cp);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_fnt_257(c->coeffs);

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_csi(c->coeffs, b_s1, (uint16_t *) mult_buf);

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_stream_add_uniform_gamma1(z, b_s1, rhoprime, (uint16_t) (nonce + i), (uint8_t *) mult_buf);
 
        if (PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_chknorm(z, GAMMA1 - BETA)) {
            goto rej;
        }

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyz_pack(sig + SEEDBYTES + i * POLYZ_PACKEDBYTES, z);
    }

    /* Compute w - cs2 */
    for (i = 0; i < K; ++i) {
        /* Unpack s2_i and convert from binary masking to arithmetic masking */
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_unpack(b_s2, sk_s2 + i * NSHARES * POLYETA_PACKEDBYTES);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_bpoly_secb2a_mod257(b_s2, COEF_257_NBITS, (uint32_t *) mult_buf);

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_decompress_chall(c, cp);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_fnt_257(c->coeffs);

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_csi(c->coeffs, b_s2, (uint16_t *) mult_buf);

        /* Compute w_i - s2_i */
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_stream_sub_from_w0(z, b_s2, &w0->vec[i], (uint8_t *) mult_buf);
        //PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_caddq(y);

        if (PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_chknorm(z, GAMMA2 - BETA)) {
            goto rej;
        }

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress_w0(&w0->vec[i], z);
    }

    /* Compute w - cs2 + ct0 and find the hint */
    n = 0;
    PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompress_chall(t2, cp);
    PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_ntt(t2);
    PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_reduce(t2);
    PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress(dc, t2);

    for (i = 0; i < K; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt0_unpack(t2, sk_t0 + i * POLYT0_PACKEDBYTES);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_ntt(t2);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compressed_pointwise_montgomery(t2, t2, dc);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_invntt_tomont(t2);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_reduce(t2);

        if (PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_chknorm(t2, GAMMA2)) {
            goto rej;
        }

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_cpoly_make_hint_and_pack(h, t2, &w0->vec[i], &n, i);

        if ((n + 1) > OMEGA) {
            goto rej;
        }
    }

    PQCLEAN_DILITHIUM2_MASK_OPTMEM_finish_hint(h, n);

    *siglen = PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES;
    return 0;
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign
*
* Description: Compute signed message.
*
* Arguments:   - uint8_t *sm: pointer to output signed message (allocated
*                             array with PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES + mlen bytes),
*                             can be equal to m
*              - size_t *smlen: pointer to output length of signed
*                               message
*              - const uint8_t *m: pointer to message to be signed
*              - size_t mlen: length of message
*              - const uint8_t *sk: pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign(uint8_t *sm,
        size_t *smlen,
        const uint8_t *m,
        size_t mlen,
        const uint8_t *sk,
        uint8_t *buf) {
    size_t i;

    for (i = 0; i < mlen; ++i) {
        sm[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
    }
    PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_signature(sm, smlen, sm + PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES, mlen, sk, buf);
    *smlen += mlen;
    return 0;
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_verify
*
* Description: Verifies signature.
*
* Arguments:   - uint8_t *m: pointer to input signature
*              - size_t siglen: length of signature
*              - const uint8_t *m: pointer to message
*              - size_t mlen: length of message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signature could be verified correctly and -1 otherwise
**************************************************/
int PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_verify(const uint8_t *sig,
        size_t siglen,
        const uint8_t *m,
        size_t mlen,
        const uint8_t *pk,
        uint8_t *buf) {
    unsigned int i, j;
    int hlen;
    const uint8_t *rho, *c, *t1;
    uint8_t *mu, *c2, *genbuf, *pw1, *h, *cbuff;
    poly *z, *w1; //*ct;
    chall_poly *cp;
    comp_poly *cw1;
    shake256incctx *state;

    /* Memory repartition */
    cw1 = (comp_poly *) buf;
    z = (poly *) ((uint8_t *) cw1 + sizeof(comp_poly));
    state = (shake256incctx *) ((uint8_t *) z + sizeof(poly));
    genbuf = (uint8_t *) state + sizeof(shake256incctx);
    cp = (chall_poly *) ((uint8_t *) genbuf + sizeof(shake256incctx) + 8);

    /* Aliases */
    //ct = z;
    w1 = z;
    mu = (uint8_t *) z;
    pw1 = (uint8_t *) z;
    c2 = (uint8_t *) z;
    cbuff = (uint8_t *) z;
    h = genbuf;

    /* Variables taken from input */
    rho = pk;
    t1 = pk + SEEDBYTES;
    c = sig;

    if (siglen != PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES) {
        return -1;
    }

    if (PQCLEAN_DILITHIUM2_MASK_OPTMEM_verify_hint(sig)) {
        return -1;
    }

    /* Compute CRH(H(rho, t1), msg) */
    shake256(mu, SEEDBYTES, pk, PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_PUBLICKEYBYTES, genbuf);

    shake256_inc_init(state);
    shake256_inc_absorb(state, mu, SEEDBYTES);
    shake256_inc_absorb(state, m, mlen);
    shake256_inc_finalize(state);
    shake256_inc_squeeze(mu, CRHBYTES, state);

    shake256_inc_init(state);
    shake256_inc_absorb(state, mu, CRHBYTES);

    PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compact_challenge(cp, c, cbuff, genbuf);

    /* Negates challenge to compute -ct0 = (-c)t0*/
    PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_negate_chall_c(cp);

    /* Streamlined call to the random oracle */
    for (i = 0; i < K; ++i){
        /* Compute w1_i */
        //PQCLEAN_DILITHIUM2_MASK_OPTMEM_cpoly_init(cw1, 0);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompress_chall(w1, cp);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_ntt(w1);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_reduce(w1);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress(cw1, w1);

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt1_unpack(w1, t1 + i * POLYT1_PACKEDBYTES);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_shiftl(w1);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_ntt(w1);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compressed_pointwise_montgomery(w1, w1, cw1);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress(cw1, w1);

        for (j = 0; j < L; ++j){
            PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyz_unpack(z, sig + SEEDBYTES + j * POLYZ_PACKEDBYTES);

            if ((i == 0) && PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_chknorm(z, GAMMA1 - BETA)) {
                return -1;
            }
            PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_ntt(z);
            PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_stream_acc_matrix_pointwise(cw1, z, rho, i, j, genbuf);
        }

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompress(w1, cw1);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_invntt_tomont(w1);
        //PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress(cw1, w1);

        /* Compute ct1 using textbook multiplication and then w1 - ct1*/
        //PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compute_ct1(ct, cp, t1 + i * POLYT1_PACKEDBYTES);
        //PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_reduce(ct);
        //PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_sub_from_cpoly(w1, cw1, ct);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_caddq(w1);
        
        hlen = PQCLEAN_DILITHIUM2_MASK_OPTMEM_unpack_hint(h, sig, i);

        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_use_hint(w1, w1, h, hlen);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyw1_pack(pw1, w1);

        /* Absorb packed representation of w1_i */
        shake256_inc_absorb(state, pw1, POLYW1_PACKEDBYTES);
    }

    shake256_inc_finalize(state);
    shake256_inc_squeeze(c2, SEEDBYTES, state);
    for (i = 0; i < SEEDBYTES; ++i) {
        if (c[i] != c2[i]) {
            return -1;
        }
    }

    return 0;
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_open
*
* Description: Verify signed message.
*
* Arguments:   - uint8_t *m: pointer to output message (allocated
*                            array with smlen bytes), can be equal to sm
*              - size_t *mlen: pointer to output length of message
*              - const uint8_t *sm: pointer to signed message
*              - size_t smlen: length of signed message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signed message could be verified correctly and -1 otherwise
**************************************************/
int PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_open(uint8_t *m,
        size_t *mlen,
        const uint8_t *sm,
        size_t smlen,
        const uint8_t *pk,
        uint8_t *buf) {
    size_t i;

    if (smlen < PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES) {
        goto badsig;
    }

    *mlen = smlen - PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES;
    if (PQCLEAN_DILITHIUM2_MASK_OPTMEM_crypto_sign_verify(sm, PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES, sm + PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES, *mlen, pk, buf)) {
        goto badsig;
    } else {
        /* All good, copy msg, return 0 */
        for (i = 0; i < *mlen; ++i) {
            m[i] = sm[PQCLEAN_DILITHIUM2_MASK_OPTMEM_CRYPTO_BYTES + i];
        }
        return 0;
    }

badsig:
    /* Signature verification failed */
    *mlen = (size_t) -1;
    for (i = 0; i < smlen; ++i) {
        m[i] = 0;
    }

    return -1;
}
