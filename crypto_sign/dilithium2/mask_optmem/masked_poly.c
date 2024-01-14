#include "masked_poly.h"
#include "randombytes.h"
#include "masked_representations.h"
#include "masked_utils.h"
#include "reduce.h"

static void pack32(uint8_t *out, const uint32_t in){
    unsigned int i = 0;

    /*Little-endian packing*/
    out[i++] = (uint8_t) ((in & 0x000000FF) >> 0);
    out[i++] = (uint8_t) ((in & 0x0000FF00) >> 8);
    out[i++] = (uint8_t) ((in & 0x00FF0000) >> 16);
    out[i++] = (uint8_t) ((in & 0xFF000000) >> 24);
}

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyveck_unmask(polyveck v[NSHARES]){
    unsigned int i;

    for (i = 1; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyveck_add(&v[0], &v[0], &v[i]); 
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecl_unmask(polyvecl v[NSHARES]){
    unsigned int i;

    for (i = 1; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyvecl_add(&v[0], &v[0], &v[i]); 
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_unmask(poly a[NSHARES]){
    unsigned int i;

    for (i = 1; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_add(&a[0], &a[0], &a[i]); 
    }
}*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_bpoly_unmask(uint32_t a[NSHARES * COEF_NBITS * N/32]){
    unsigned int i, j, d;

    for (d = 1; d < NSHARES; ++d){
        for (i = 0; i < N/32; ++i) {
            for (j = 0; j < COEF_257_NBITS; ++j) {
                a[i*COEF_257_NBITS*NSHARES + j] ^= a[d * COEF_257_NBITS + i*COEF_257_NBITS*NSHARES + j];
            }
        } 
    }
}

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_negate_unmask(poly a[NSHARES]){
    unsigned int i;

    for (i = 1; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_sub(&a[0], &a[0], &a[i]); 
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecl_ntt(polyvecl v[NSHARES]){
    unsigned int i;

    for (i = 0; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyvecl_ntt(&v[i]); 
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyveck_ntt(polyveck v[NSHARES]){
    unsigned int i;

    for (i = 0; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyveck_ntt(&v[i]); 
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_ntt(poly a[NSHARES]){
    unsigned int i;

    for (i = 0; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_ntt(&a[i]); 
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyveck_invntt_tomont(polyveck v[NSHARES]){
    unsigned int i;

    for (i = 0; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyveck_invntt_tomont(&v[i]); 
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecl_invntt_tomont(polyvecl v[NSHARES]){
    unsigned int i;

    for (i = 0; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyvecl_invntt_tomont(&v[i]); 
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_invntt_tomont(poly a[NSHARES]){
    unsigned int i;

    for (i = 0; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_invntt_tomont(&a[i]); 
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyveck_pointwise_poly_montgomery(polyveck r[NSHARES], const poly *a, const polyveck v[NSHARES]){
    unsigned int i;

    for (i = 0; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyveck_pointwise_poly_montgomery(&r[i], a, &v[i]); 
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecl_pointwise_poly_montgomery(polyvecl r[NSHARES], const poly *a, const polyvecl v[NSHARES]){
    unsigned int i;

    for (i = 0; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyvecl_pointwise_poly_montgomery(&r[i], a, &v[i]); 
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_pointwise_poly_montgomery(poly r[NSHARES], const poly *a, const poly v[NSHARES]){
    unsigned int i;

    for (i = 0; i < NSHARES; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_pointwise_montgomery(&r[i], a, &v[i]); 
    }
}*/

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecleta_mask(BitslicePolyveclEta r[NSHARES], polyvecl *p, uint32_t buf[POLYETA_PACKEDBYTES/4]){
    unsigned int i;

    for (i = 0; i < L; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_mask(&r[0].BPoly[i], &p->vec[i], L, buf);
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecketa_mask(BitslicePolyveckEta r[NSHARES], polyveck *p, uint32_t buf[POLYETA_PACKEDBYTES/4]){
    unsigned int i;

    for (i = 0; i < K; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_mask(&r[0].BPoly[i], &p->vec[i], K, buf);
    }
}*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_unpack(uint32_t *r, const uint8_t *p){
    unsigned d, i, j;

    /*Little-endian unpacking*/
    for (d = 0; d < NSHARES; ++d){
        for (i = 0; i < N/32; ++i){
            for (j = 0; j < ETA_NBITS; ++j){
                r[d * COEF_257_NBITS + i * COEF_257_NBITS * NSHARES + j] = (uint32_t)(p[4*(d * ETA_NBITS * N/32 + i * ETA_NBITS + j) + 0]);
                r[d * COEF_257_NBITS + i * COEF_257_NBITS * NSHARES + j] |= (uint32_t)(p[4*(d * ETA_NBITS * N/32 + i * ETA_NBITS + j) + 1]) << 8;
                r[d * COEF_257_NBITS + i * COEF_257_NBITS * NSHARES + j] |= (uint32_t)(p[4*(d * ETA_NBITS * N/32 + i * ETA_NBITS + j) + 2]) << 16;
                r[d * COEF_257_NBITS + i * COEF_257_NBITS * NSHARES + j] |= (uint32_t)(p[4*(d * ETA_NBITS * N/32 + i * ETA_NBITS + j) + 3]) << 24;
            }
            for (j = ETA_NBITS; j < COEF_257_NBITS; ++j) {
                r[d * COEF_257_NBITS + i * COEF_257_NBITS * NSHARES + j] = 0;
            }
        }
    }
}

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_mask(BitslicePolyEta *r, const poly *p, const uint32_t mask_stride, uint32_t buf[32]){
    unsigned int d, i, j;

    polyeta_canonical2bitslice_opt(r, p, buf);

    for(d = 1; d < NSHARES; ++d){
        randombytes_masking((uint8_t*) buf, POLYETA_PACKEDBYTES);
        for (i = 0; i < ETA_NBITS; ++i){
            for (j = 0; j < (N/32); ++j){
                r->coeffs[i][j] ^= buf[i*(N/32) + j];
                r[d * mask_stride/4].coeffs[i][j] = buf[i*(N/32) + j];
            }
        }
    }
}*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_pack(uint8_t *b, const poly *a, const size_t mask_stride, uint32_t buf[32]){
    unsigned int d, i;
    uint32_t t, r;

    polyeta_canonical2bitslice_opt(buf, a);

    for (i = 0; i < ETA_NBITS * (N/32); ++i){
        t = buf[i];
        for(d = 1; d < NSHARES; ++d){
            r = rand32();
            t ^= r;
            pack32(b + (4 * i + d * mask_stride), r);
        }
        
        pack32(b + (4 * i), t);
    }
}

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecleta_unmask(polyvecl *r, const BitslicePolyveclEta *p, uint32_t buf[32 + POLYETA_PACKEDBYTES/4]){
    unsigned int i;

    for (i = 0; i < L; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_unmask(&r->vec[i], &p->BPoly[i], L, buf);
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecketa_unmask(polyveck *r, const BitslicePolyveckEta *p, uint32_t buf[32 + POLYETA_PACKEDBYTES/4]){
    unsigned int i;

    for (i = 0; i < K; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_unmask(&r->vec[i], &p->BPoly[i], K, buf);
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_unmask(poly *r, const BitslicePolyEta *p, const size_t mask_stride, uint32_t buf[32 + POLYETA_PACKEDBYTES/4]){
    unsigned int i, j, d;
    BitslicePolyEta *t;

    t = (BitslicePolyEta *) (buf + 32);

    for (i = 0; i < ETA_NBITS; ++i){
        for (j = 0; j < (N/32); ++j)
        {
            t->coeffs[i][j] = p->coeffs[i][j];

            for (d = 1; d < NSHARES; ++d){
                t->coeffs[i][j] ^= p[d * mask_stride/4].coeffs[i][j];
            }
        }
    }

    polyeta_bitslice2canonical_opt(r, t, buf);    
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecketa_expand(polyveck r[NSHARES], const BitslicePolyveckEta v[NSHARES]) {
    unsigned int i, d;

    for (d = 0; d < NSHARES; ++d){
        for (i = 0; i < K; ++i){
            PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_expand(&r[d].vec[i], &v[d].BPoly[i]);
        }
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecleta_expand(polyvecl r[NSHARES], const BitslicePolyveclEta v[NSHARES]) {
    unsigned int i, d;

    for (d = 0; d < NSHARES; ++d){
        for (i = 0; i < L; ++i){
            PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_expand(&r[d].vec[i], &v[d].BPoly[i]);
        }
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_expand(poly *r, const BitslicePolyEta *a){
    unsigned int i, j;

    for (i = 0; i < N/32; ++i){
        for (j = 0; j < ETA_NBITS; ++j){
            r->coeffs[i*32 + j] = a->coeffs[j][i];
        }
    }

    //Set remaining bits to 0
    for (i = 0; i < N/32; ++i){
        for (j = ETA_NBITS; j < 32; ++j){
            r->coeffs[i*32 + j] = 0;
        }
    }
}*/

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecketa_secb2a_modp(polyveck v[NSHARES], 
        uint32_t buf[2 * BSSIZE * NSHARES + 2 * COEF_NBITS * NSHARES + 2 * (COEF_NBITS + 1) * NSHARES]){
    unsigned int i;

    for (i = 0; i < K; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_secb2a_modp(&v->vec[i], POLYVECK_BYTES/4, buf);
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyvecleta_secb2a_modp(polyvecl v[NSHARES],
        uint32_t buf[2 * BSSIZE * NSHARES + 2 * COEF_NBITS * NSHARES + 2 * (COEF_NBITS + 1) * NSHARES]){
    unsigned int i;

    for (i = 0; i < L; ++i){
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyeta_secb2a_modp(&v->vec[i], POLYVECL_BYTES/4, buf);
    }
}*/

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_bpoly_secb2a_modp(uint32_t *a, size_t mask_stride, 
        uint32_t buf[2 * BSSIZE * NSHARES + 2 * COEF_NBITS * NSHARES + 2 * (COEF_NBITS + 1) * NSHARES]){
    unsigned int i;

    for(i = 0; i < N/32; ++i){
        secb2a_modp(NSHARES, 23, Q, &a[i * COEF_NBITS], mask_stride/4, 1, &a[i * COEF_NBITS], mask_stride/4, 1, buf);
    }
}*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_bpoly_secb2a_mod257(uint32_t *a, size_t mask_stride, 
        uint32_t buf[COEF_257_NBITS * (NSHARES - 1) + (COEF_257_NBITS + 1) * NSHARES + 4 * NSHARES + 4]){
    unsigned int i;

    for(i = 0; i < N/32; ++i){
        secb2a_mod257(NSHARES, QP, &a[i * COEF_257_NBITS * NSHARES], mask_stride, 1, buf);
    }
}

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_bitslice2dense(poly *v, uint32_t *u) {
    unsigned d;
    
    for (d = 0; d < NSHARES; ++d) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_bitslice2dense(&v[d], &u[d * COEF_NBITS * N/32]); 
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_polyveck_negate(polyveck v[NSHARES]){
    unsigned int i, d;

    for (d = 0; d < NSHARES; ++d){
        for (i = 0; i < K; ++i){
            PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_negate(&v[d].vec[i]);
        }
    }
}*/

static void load_coefficients(uint32_t *a, uint32_t *coeff_buf) {
    seca2a_centered_modp(NSHARES, 257, a, COEF_257_NBITS, 1, Q, coeff_buf, COEF_NBITS, 1, coeff_buf);
    masked_inplace_bitslice2dense_opt(NSHARES, COEF_NBITS, coeff_buf);
}

static int32_t add_and_unmask(int32_t r, uint32_t *a) {
    for (size_t i = 0; i < NSHARES; ++i){
        r += (int32_t) a[32*i];
        r = PQCLEAN_DILITHIUM2_MASK_OPTMEM_reduce32(r);
    }
    return r;
}

static void unpack_gamma1_and_add(poly *r, uint32_t *a, uint32_t *coeff_buf, const uint8_t *buf, unsigned int buflen, unsigned int *ctr){
    unsigned int i, bound;
    uint32_t t;

    bound = N-(*ctr);
    for (i = 0; (i < bound/4) && (i < buflen/9); ++i){
        if (((*ctr) % 32) == 0){
            load_coefficients(&a[((*ctr)/32) * COEF_257_NBITS * NSHARES], coeff_buf);
        }

        t = buf[9 * i + 0];
        t |= (uint32_t)buf[9 * i + 1] << 8;
        t |= (uint32_t)buf[9 * i + 2] << 16;
        t &= 0x3FFFF;
        r->coeffs[*ctr] = add_and_unmask(GAMMA1 - t, &coeff_buf[(*ctr) % 32]);

        (*ctr)++;

        t  = buf[9 * i + 2] >> 2;
        t |= (uint32_t)buf[9 * i + 3] << 6;
        t |= (uint32_t)buf[9 * i + 4] << 14;
        t &= 0x3FFFF;
        r->coeffs[*ctr] = add_and_unmask(GAMMA1 - t, &coeff_buf[(*ctr) % 32]);

        (*ctr)++;

        t  = buf[9 * i + 4] >> 4;
        t |= (uint32_t)buf[9 * i + 5] << 4;
        t |= (uint32_t)buf[9 * i + 6] << 12;
        t &= 0x3FFFF;
        r->coeffs[*ctr] = add_and_unmask(GAMMA1 - t, &coeff_buf[(*ctr) % 32]);

        (*ctr)++;

        t  = buf[9 * i + 6] >> 6;
        t |= (uint32_t)buf[9 * i + 7] << 2;
        t |= (uint32_t)buf[9 * i + 8] << 10;
        t &= 0x3FFFF;
        r->coeffs[*ctr] = add_and_unmask(GAMMA1 - t, &coeff_buf[(*ctr) % 32]);

        (*ctr)++;
    }
}


void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_stream_add_uniform_gamma1(poly *r,
        uint32_t a[NSHARES * COEF_257_NBITS * N/32],
        const uint8_t seed[CRHBYTES],
        uint16_t nonce,
        uint8_t buf[sizeof(stream256_state) + 16 + 4*((2*COEF_NBITS + 5)*NSHARES - COEF_NBITS)]) {
    
    unsigned int ctr, off, k;
    unsigned int buflen = STREAM256_BLOCKBYTES;
    uint32_t *coeff_buf;
    stream256_state *state = (stream256_state *) (buf + 16);
    coeff_buf = (uint32_t *) (buf + sizeof(stream256_state) + 16);
    buf = buf + 7;

    ctr = 0;

    stream256_init(state, seed, nonce);
    stream256_squeezeblocks(buf + 9, 1, state);

    unpack_gamma1_and_add(r, a, coeff_buf, buf + 9, buflen, &ctr);

    while (ctr < N) {
        off = buflen % 9;
        for (k = 9 - off; k < 9; ++k) {
            buf[k] = buf[STREAM256_BLOCKBYTES + k];
        }

        stream256_squeezeblocks(buf + 9, 1, state);
        buflen = STREAM256_BLOCKBYTES + off;
        unpack_gamma1_and_add(r, a, coeff_buf, buf + 9 - off, buflen, &ctr);
    }

}

static int32_t sub_and_unmask(int32_t r, uint32_t *a) {
    for (size_t i = 0; i < NSHARES; ++i){
        r -= (int32_t) a[32*i];
        r = PQCLEAN_DILITHIUM2_MASK_OPTMEM_reduce32(r);
    }
    return r;
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_stream_sub_from_w0(poly *r,
        uint32_t a[NSHARES * COEF_257_NBITS * N/32], 
        comp_w0_poly *w0,
        uint8_t buf[4*((2*COEF_NBITS + 5)*NSHARES - COEF_NBITS)]) {
    
    unsigned int i,j;
    int32_t t;
    uint32_t *coeff_buf;

    coeff_buf = (uint32_t *) buf;

    for (i = 0; i < N/32; ++i){
        load_coefficients(&a[i * COEF_257_NBITS * NSHARES], coeff_buf);

        for (j = 0; j < 32; ++j){
            t = decompress_w0(&w0->coeffs[5*(16*i + j/2)], (j & 0x1));
            t = sub_and_unmask(t, &coeff_buf[j]);

            r->coeffs[32*i + j] = t; 
        }
    }
}