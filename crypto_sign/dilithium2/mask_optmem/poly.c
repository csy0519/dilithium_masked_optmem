#include "ntt.h"
#include "params.h"
#include "poly.h"
#include "reduce.h"
#include "rounding.h"
#include "masked_representations.h"
#include <stdint.h>

#define DBENCH_START()
#define DBENCH_STOP(t)

void compress_q(uint8_t *b, const uint32_t a){
    b[0] = (uint8_t)((a >>  0) & 0xFF);
    b[1] = (uint8_t)((a >>  8) & 0xFF);
    b[2] = (uint8_t)((a >> 16) & 0xFF);
}

int32_t decompress_q(const uint8_t *a){
    int32_t r;

    r = (int32_t) a[0] << 8;
    r |= (uint32_t) a[1] << 16;
    r |= (uint32_t) a[2] << 24;

    return (r >> 8);
}

static void compress_w0(uint8_t *b, const int32_t w0, const size_t odd) {
    b[0 + 3*odd] = (uint8_t)((w0 >> 0) & 0xFF);
    b[1 + 3*odd] = (uint8_t)((w0 >> 8) & 0xFF);

    if (odd == 1) {
        b[2] = (b[2] & 0x1F) | ((uint8_t)((w0 >> 11) & 0xE0));
    } else {
        b[2] = (b[2] & 0xE3) | ((uint8_t)((w0 >> 14) & 0x1C));
    }
}

void compress_w(uint8_t *b, const int32_t w1, const int32_t w0, const size_t odd){
    compress_w0(b, w0, odd);

    if (odd == 1) {
        b[2] = (b[2] & 0xFD) | ((w1 == 0) << 1);
    } else {
        b[2] = (b[2] & 0xFE) | (w1 == 0);
    }
}

static int32_t decompress_w1(const uint8_t *a, const size_t odd){
    return ((a[2] & (1 << odd)) >> odd);
}

int32_t decompress_w0(const uint8_t *a, const size_t odd){
    int32_t r;

    r = (int32_t) a[0 + 3*odd];
    r |= (int32_t) a[1 + 3*odd] << 8;

    if (odd == 1) {
        r |= (int32_t)(((uint32_t) a[2] & 0xE0) << 24) >> 13;
    } else {
        r |= (int32_t)(((uint32_t) a[2] & 0x1C) << 27) >> 13;
    }

    return r;
}

static void add_24bit_reduce(uint8_t *c, const uint8_t *a, const uint32_t b){
    int32_t t;

    t = decompress_q(a);
    t = PQCLEAN_DILITHIUM2_MASK_OPTMEM_reduce32(t + b);
    compress_q(c, t); 
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_reduce
*
* Description: Inplace reduction of all coefficients of polynomial to
*              representative in [-6283009,6283007].
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_reduce(poly *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        a->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_reduce32(a->coeffs[i]);
    }

    DBENCH_STOP(*tred);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_caddq
*
* Description: For all coefficients of in/out polynomial add Q if
*              coefficient is negative.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_caddq(poly *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        a->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_caddq(a->coeffs[i]);
    }

    DBENCH_STOP(*tred);
}

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_freeze(poly *a) {
    unsigned int i;

    for (i = 0; i < N; ++i) {
        a->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_freeze(a->coeffs[i]);
    }
}*/

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_add
*
* Description: Add polynomials. No modular reduction is performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first summand
*              - const poly *b: pointer to second summand
**************************************************/
/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_add(poly *c, const poly *a, const poly *b)  {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        c->coeffs[i] = a->coeffs[i] + b->coeffs[i];
    }

    DBENCH_STOP(*tadd);
}*/

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_sub
*
* Description: Subtract polynomials. No modular reduction is
*              performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial to be
*                               subtraced from first input polynomial
**************************************************/
/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_sub(poly *c, const poly *a, const poly *b) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        c->coeffs[i] = a->coeffs[i] - b->coeffs[i];
    }

    DBENCH_STOP(*tadd);
}*/

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_shiftl
*
* Description: Multiply polynomial by 2^D without modular reduction. Assumes
*              input coefficients to be less than 2^{31-D} in absolute value.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_shiftl(poly *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        a->coeffs[i] <<= D;
    }

    DBENCH_STOP(*tmul);

}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_sub_eta(poly *c, const poly *b){
    unsigned int i;
    for (i = 0; i < N; ++i) {
        c->coeffs[i] = b->coeffs[i] - ETA;
    }
}

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_negate(poly *a){
    unsigned int i;

    for (i = 0; i < N; ++i){
        a->coeffs[i] = -a->coeffs[i];
    }
}*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_negate_chall_c(chall_poly *c){
    
    c->signs = ~c->signs;
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_ntt
*
* Description: Inplace forward NTT. Coefficients can grow by
*              8*Q in absolute value.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_ntt(poly *a) {
    DBENCH_START();

    PQCLEAN_DILITHIUM2_MASK_OPTMEM_ntt(a->coeffs);

    DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_invntt_tomont
*
* Description: Inplace inverse NTT and multiplication by 2^{32}.
*              Input coefficients need to be less than Q in absolute
*              value and output coefficients are again bounded by Q.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_invntt_tomont(poly *a) {
    DBENCH_START();

    PQCLEAN_DILITHIUM2_MASK_OPTMEM_invntt_tomont(a->coeffs);

    DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_pointwise_montgomery
*
* Description: Pointwise multiplication of polynomials in NTT domain
*              representation and multiplication of resulting polynomial
*              by 2^{-32}.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_pointwise_montgomery(poly *c, const poly *a, const poly *b) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        c->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce((int64_t)a->coeffs[i] * b->coeffs[i]);
    }

    DBENCH_STOP(*tmul);
}*/


void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compressed_pointwise_montgomery(poly *c, const poly *a, const comp_poly *b) {
    unsigned int i;
    int32_t t;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        t = decompress_q(&b->coeffs[3*i]);
        c->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce((int64_t)a->coeffs[i] * t);
    }

    DBENCH_STOP(*tmul);
}

/*static int32_t polyt1_unpack_coeff(const uint8_t *a, size_t i) {
    size_t i1, i0;

    i1 = i / 4;
    i0 = i - 4*i1;
    return ((a[5 * i1 + i0] >> (2 * i0)) | ((uint32_t)a[5 * i1 + i0 + 1] << (8 - 2*i0))) & 0x3FF;
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compute_ct1(poly *r, const chall_poly *c, const uint8_t *a) {
    unsigned int i, j;
    int32_t t;
    uint8_t sign;

    for (j = 0; j < N; ++j){
        t = polyt1_unpack_coeff(a, j) << D;
        sign = (c->signs & 0x1) ^ ((j + c->index[0])/N);
        r->coeffs[(j + c->index[0]) % N] = (1 - 2 * sign) * t;
    }

    for (i = 1; i < TAU; ++i){
        for (j = 0; j < N; ++j){
            t = polyt1_unpack_coeff(a, j) << D;
            sign = ((c->signs >> i) & 0x1) ^ ((j + c->index[i])/N);
            r->coeffs[(j + c->index[i]) % N] += (1 - 2 * sign) * t;
        }
    }
}

static int32_t polyt0_unpack_coeff(const uint8_t *a, size_t i) {
    size_t i1, i0;
    uint32_t t = 0;

    i1 = i / 8;
    i0 = i - 8*i1;

    switch (i0)
    {
        case 0:
            t  = a[13 * i1 + 0];
            t |= (uint32_t)a[13 * i1 + 1] << 8;
            break;

        case 1:
            t  = a[13 * i1 + 1] >> 5;
            t |= (uint32_t)a[13 * i1 + 2] << 3;
            t |= (uint32_t)a[13 * i1 + 3] << 11;
            break;

        case 2:
            t  = a[13 * i1 + 3] >> 2;
            t |= (uint32_t)a[13 * i1 + 4] << 6;
            break;

        case 3:
            t  = a[13 * i1 + 4] >> 7;
            t |= (uint32_t)a[13 * i1 + 5] << 1;
            t |= (uint32_t)a[13 * i1 + 6] << 9;
            break;

        case 4:
            t  = a[13 * i1 + 6] >> 4;
            t |= (uint32_t)a[13 * i1 + 7] << 4;
            t |= (uint32_t)a[13 * i1 + 8] << 12;
            break;;

        case 5:
            t  = a[13 * i1 + 8] >> 1;
            t |= (uint32_t)a[13 * i1 + 9] << 7;
            break;

        case 6:
            t  = a[13 * i1 + 9] >> 6;
            t |= (uint32_t)a[13 * i1 + 10] << 2;
            t |= (uint32_t)a[13 * i1 + 11] << 10;
            break;
        
        case 7:
            t  = a[13 * i1 + 11] >> 3;
            t |= (uint32_t)a[13 * i1 + 12] << 5;
            break;

    }

    t &= 0x1FFF;
    return ((1 << (D - 1)) - t);
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compute_ct0(poly *r, const chall_poly *c, const uint8_t *a) {
    unsigned int i, j;
    int32_t t;
    uint8_t sign;

    for (j = 0; j < N; ++j){
        t = polyt0_unpack_coeff(a, j);
        sign = (c->signs & 0x1) ^ ((j + c->index[0])/N);
        r->coeffs[(j + c->index[0]) % N] = (1 - 2 * sign) * t;
    }

    for (i = 1; i < TAU; ++i){
        for (j = 0; j < N; ++j){
            t = polyt0_unpack_coeff(a, j);
            sign = ((c->signs >> i) & 0x1) ^ ((j + c->index[i])/N);
            r->coeffs[(j + c->index[i]) % N] += (1 - 2 * sign) * t;
        }
    }
}*/

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_power2round
*
* Description: For all coefficients c of the input polynomial,
*              compute c0, c1 such that c mod Q = c1*2^D + c0
*              with -2^{D-1} < c0 <= 2^{D-1}. Assumes coefficients to be
*              standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients c0
*              - const poly *a: pointer to input polynomial
**************************************************/
/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_power2round(poly *a1, poly *a0, const poly *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        a1->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_power2round(&a0->coeffs[i], a->coeffs[i]);
    }

    DBENCH_STOP(*tround);
}*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_power2round_a1(poly *a1, const poly *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        a1->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_power2round_a1(a->coeffs[i]);
    }

    DBENCH_STOP(*tround);
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_power2round_a0(poly *a0, const poly *a1, const comp_poly *a) {
    unsigned int i;
    int32_t t;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        t = decompress_q(&a->coeffs[3*i]);
        a0->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_power2round_a0(a1->coeffs[i], t);
    }

    DBENCH_STOP(*tround);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompose
*
* Description: For all coefficients c of the input polynomial,
*              compute high and low bits c0, c1 such c mod Q = c1*ALPHA + c0
*              with -ALPHA/2 < c0 <= ALPHA/2 except c1 = (Q-1)/ALPHA where we
*              set c1 = 0 and -ALPHA/2 <= c0 = c mod Q - Q < 0.
*              Assumes coefficients to be standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients c0
*              - const poly *a: pointer to input polynomial
**************************************************/
/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompose(poly *a1, poly *a0, const poly *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        a1->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_decompose(&a0->coeffs[i], a->coeffs[i]);
    }

    DBENCH_STOP(*tround);
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompose_a1(poly *a1, const poly *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        a1->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_decompose_a1(a->coeffs[i]);
    }

    DBENCH_STOP(*tround);
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompose_a0(poly *a0, const poly *a1, const comp_poly *a) {
    unsigned int i;
    int32_t t;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        t = decompress_q(&a->coeffs[3*i]);
        a0->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_decompose_a0(a1->coeffs[i], t);
    }

    DBENCH_STOP(*tround);
}*/

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_make_hint
*
* Description: Compute hint polynomial. The coefficients of which indicate
*              whether the low bits of the corresponding coefficient of
*              the input polynomial overflow into the high bits.
*
* Arguments:   - poly *h: pointer to output hint polynomial
*              - const poly *a0: pointer to low part of input polynomial
*              - const poly *a1: pointer to high part of input polynomial
*
* Returns number of 1 bits.
**************************************************/
/*unsigned int PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_make_hint(poly *h, const poly *a0, const poly *a1) {
    unsigned int i, s = 0;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        h->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_make_hint(a0->coeffs[i], a1->coeffs[i]);
        s += h->coeffs[i];
    }

    DBENCH_STOP(*tround);
    return s;
}*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_cpoly_make_hint_and_pack(uint8_t *h, const poly *a, const comp_w0_poly *b, uint32_t *n, size_t index) {
    unsigned int i;
    int32_t t, a0, a1;
    DBENCH_START();

    for (i = 0; i < N && *n < OMEGA; ++i) {
        a0 = decompress_w0(&b->coeffs[5*(i/2)], (i & 0x1));
        a1 = decompress_w1(&b->coeffs[5*(i/2)], (i & 0x1));
        t = PQCLEAN_DILITHIUM2_MASK_OPTMEM_make_hint(a0 + a->coeffs[i], a1);
        h[*n] = i;
        *n += t;
    }

    h[OMEGA + index] = *n;

    DBENCH_STOP(*tround);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_use_hint
*
* Description: Use hint polynomial to correct the high bits of a polynomial.
*
* Arguments:   - poly *b: pointer to output polynomial with corrected high bits
*              - const poly *a: pointer to input polynomial
*              - const poly *h: pointer to input hint polynomial
**************************************************/
/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_use_hint(poly *b, const poly *a, const poly *h) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N; ++i) {
        b->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_use_hint(a->coeffs[i], h->coeffs[i]);
    }

    DBENCH_STOP(*tround);
}*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_use_hint(poly *b, const poly *a, const uint8_t *h, const int hlen) {
    unsigned int i;
    int k;
    uint32_t hint;
    DBENCH_START();

    k = 0;
    for (i = 0; i < N; ++i) {
        hint = (k < hlen && i == h[k]) ? 1 : 0;
        b->coeffs[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_use_hint(a->coeffs[i], hint);
        k += hint;
    }

    DBENCH_STOP(*tround);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_chknorm
*
* Description: Check infinity norm of polynomial against given bound.
*              Assumes input coefficients were reduced by PQCLEAN_DILITHIUM2_MASK_OPTMEM_reduce32().
*
* Arguments:   - const poly *a: pointer to polynomial
*              - int32_t B: norm bound
*
* Returns 0 if norm is strictly smaller than B <= (Q-1)/8 and 1 otherwise.
**************************************************/
int PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_chknorm(const poly *a, int32_t B) {
    unsigned int i;
    int32_t t;
    DBENCH_START();

    if (B > (Q - 1) / 8) {
        return 1;
    }

    /* It is ok to leak which coefficient violates the bound since
       the probability for each coefficient is independent of secret
       data but we must not leak the sign of the centralized representative. */
    for (i = 0; i < N; ++i) {
        /* Absolute value */
        t = a->coeffs[i] >> 31;
        t = a->coeffs[i] - (t & 2 * a->coeffs[i]);

        if (t >= B) {
            DBENCH_STOP(*tsample);
            return 1;
        }
    }

    DBENCH_STOP(*tsample);
    return 0;
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_cpoly_init(comp_poly *a, const uint32_t val){
    unsigned int i;
    
    for (i = 0; i < N; i++) {
        compress_q(&a->coeffs[3*i], val);
    } 
}

/*int PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_lowbits_chknorm(const poly *a, int32_t B) {
    unsigned int i;
    int32_t t, a0;
    DBENCH_START();

    if (B > (Q - 1) / 8) {
        return 1;
    }

    / It is ok to leak which coefficient violates the bound since
       the probability for each coefficient is independent of secret
       data but we must not leak the sign of the centralized representative. /
    for (i = 0; i < N; ++i) {
        /Absolute value /
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_decompose(&a0, a->coeffs[i]);

        t = a0 >> 31;
        t = a0 - (t & 2 * a0);

        if (t >= B) {
            return 1;
        }
    }

    return 0;
}
*/
/*************************************************
* Name:        rej_uniform
*
* Description: Sample uniformly random coefficients in [0, Q-1] by
*              performing rejection sampling on array of random bytes.
*
* Arguments:   - int32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
/*static unsigned int rej_uniform(int32_t *a,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen) {
    unsigned int ctr, pos;
    uint32_t t;
    DBENCH_START();

    ctr = pos = 0;
    while (ctr < len && pos + 3 <= buflen) {
        t  = buf[pos++];
        t |= (uint32_t)buf[pos++] << 8;
        t |= (uint32_t)buf[pos++] << 16;
        t &= 0x7FFFFF;

        if (t < Q) {
            a[ctr++] = t;
        }
    }

    DBENCH_STOP(*tsample);
    return ctr;
}
*/
static unsigned int stream_acc_pointwise_montgomery(uint8_t *b,
                                const int32_t *a,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen) {

    unsigned int ctr, pos;
    uint32_t t;

    ctr = pos = 0;
    while (ctr < len && pos + 3 <= buflen) {
        t  = buf[pos++];
        t |= (uint32_t)buf[pos++] << 8;
        t |= (uint32_t)buf[pos++] << 16;
        t &= 0x7FFFFF;

        if (t < Q) {
            t = PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce((int64_t)a[ctr] * t);
            add_24bit_reduce(b + 3 * ctr, b + 3 * ctr, t);
            ctr++;
        }
    }

    return ctr;
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_uniform
*
* Description: Sample polynomial with uniformly random coefficients
*              in [0,Q-1] by performing rejection sampling on the
*              output stream of SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length SEEDBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_uniform(poly *a,
        const uint8_t seed[SEEDBYTES],
        uint16_t nonce,
        uint8_t buf[sizeof(stream128_state) + POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES + 2]) {
    unsigned int i, ctr, off;
    unsigned int buflen = POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES;
    stream128_state *state = (stream128_state*) buf;
    buf = buf + sizeof(stream128_state);

    stream128_init(state, seed, nonce);
    stream128_squeezeblocks(buf, POLY_UNIFORM_NBLOCKS, state);

    ctr = rej_uniform(a->coeffs, N, buf, buflen);

    while (ctr < N) {
        off = buflen % 3;
        for (i = 0; i < off; ++i) {
            buf[i] = buf[buflen - off + i];
        }

        stream128_squeezeblocks(buf + off, 1, state);
        buflen = STREAM128_BLOCKBYTES + off;
        ctr += rej_uniform(a->coeffs + ctr, N - ctr, buf, buflen);
    }
}*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_stream_acc_matrix_pointwise(comp_poly *b, 
        const poly *a, const uint8_t rho[SEEDBYTES], const size_t i, const size_t j,
        uint8_t buf[sizeof(stream128_state) + 8]) {
    unsigned int k, ctr, off;
    unsigned int buflen = STREAM128_BLOCKBYTES;
    uint16_t nonce = (uint16_t)((i << 8) + j);
    stream128_state *state = (stream128_state*) (buf + 8);
    buf = buf + 6;

    stream128_init(state, rho, nonce);
    stream128_squeezeblocks(buf + 2, 1, state);

    ctr = stream_acc_pointwise_montgomery(b->coeffs, a->coeffs, N, buf + 2, buflen);

    while (ctr < N) {
        off = buflen % 3;
        for (k = 3 - off; k < 3; ++k) {
            buf[k] = buf[STREAM128_BLOCKBYTES + k];
        }

        stream128_squeezeblocks(buf + 2, 1, state);
        buflen = STREAM128_BLOCKBYTES + off;
        ctr += stream_acc_pointwise_montgomery(b->coeffs + 3 * ctr, a->coeffs + ctr, N - ctr, buf + 2 - off, buflen);
    }
}

/*************************************************
* Name:        rej_eta
*
* Description: Sample uniformly random coefficients in [-ETA, ETA] by
*              performing rejection sampling on array of random bytes.
*
* Arguments:   - int32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_eta(int32_t *a,
                            unsigned int len,
                            const uint8_t *buf,
                            unsigned int buflen) {
    unsigned int ctr, pos;
    uint32_t t0, t1;
    DBENCH_START();

    ctr = pos = 0;
    while (ctr < len && pos < buflen) {
        t0 = buf[pos] & 0x0F;
        t1 = buf[pos++] >> 4;

        if (t0 < 15) {
            t0 = t0 - (205 * t0 >> 10) * 5;
            a[ctr++] = 2 - t0;
        }
        if (t1 < 15 && ctr < len) {
            t1 = t1 - (205 * t1 >> 10) * 5;
            a[ctr++] = 2 - t1;
        }
    }

    DBENCH_STOP(*tsample);
    return ctr;
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_uniform_eta
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-ETA,ETA] by performing rejection sampling on the
*              output stream from SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length CRHBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_uniform_eta(poly *a,
        const uint8_t seed[CRHBYTES],
        uint16_t nonce,
        uint8_t buf[sizeof(stream256_state) + POLY_UNIFORM_ETA_NBLOCKS * STREAM256_BLOCKBYTES]) {
    unsigned int ctr;
    unsigned int buflen = POLY_UNIFORM_ETA_NBLOCKS * STREAM256_BLOCKBYTES;
    stream256_state *state = (stream256_state *) buf;

    stream256_init(state, seed, nonce);
    stream256_squeezeblocks(buf, POLY_UNIFORM_ETA_NBLOCKS, state);

    ctr = rej_eta(a->coeffs, N, buf, buflen);

    while (ctr < N) {
        stream256_squeezeblocks(buf, 1, state);
        ctr += rej_eta(a->coeffs + ctr, N - ctr, buf, STREAM256_BLOCKBYTES);
    }
}

static unsigned int unpack_gamma1(int32_t *r, unsigned int len, const uint8_t *buf, unsigned int buflen){
    unsigned int i;

    for (i = 0; (i < len/4) && (i < buflen/9); ++i){
        r[4 * i + 0]  = buf[9 * i + 0];
        r[4 * i + 0] |= (uint32_t)buf[9 * i + 1] << 8;
        r[4 * i + 0] |= (uint32_t)buf[9 * i + 2] << 16;
        r[4 * i + 0] &= 0x3FFFF;

        r[4 * i + 1]  = buf[9 * i + 2] >> 2;
        r[4 * i + 1] |= (uint32_t)buf[9 * i + 3] << 6;
        r[4 * i + 1] |= (uint32_t)buf[9 * i + 4] << 14;
        r[4 * i + 1] &= 0x3FFFF;

        r[4 * i + 2]  = buf[9 * i + 4] >> 4;
        r[4 * i + 2] |= (uint32_t)buf[9 * i + 5] << 4;
        r[4 * i + 2] |= (uint32_t)buf[9 * i + 6] << 12;
        r[4 * i + 2] &= 0x3FFFF;

        r[4 * i + 3]  = buf[9 * i + 6] >> 6;
        r[4 * i + 3] |= (uint32_t)buf[9 * i + 7] << 2;
        r[4 * i + 3] |= (uint32_t)buf[9 * i + 8] << 10;
        r[4 * i + 3] &= 0x3FFFF;

        r[4 * i + 0] = GAMMA1 - r[4 * i + 0];
        r[4 * i + 1] = GAMMA1 - r[4 * i + 1];
        r[4 * i + 2] = GAMMA1 - r[4 * i + 2];
        r[4 * i + 3] = GAMMA1 - r[4 * i + 3];
    }

    return 4*i;
}

/*static unsigned int unpack_gamma1_and_add(int32_t *r, unsigned int len, const uint8_t *buf, unsigned int buflen){
    unsigned int i;
    uint32_t t;

    for (i = 0; (i < len/4) && (i < buflen/9); ++i){
        t = buf[9 * i + 0];
        t |= (uint32_t)buf[9 * i + 1] << 8;
        t |= (uint32_t)buf[9 * i + 2] << 16;
        t &= 0x3FFFF;
        r[4*i + 0] += GAMMA1 - t;


        t  = buf[9 * i + 2] >> 2;
        t |= (uint32_t)buf[9 * i + 3] << 6;
        t |= (uint32_t)buf[9 * i + 4] << 14;
        t &= 0x3FFFF;
        r[4*i + 1] += GAMMA1 - t;

        t  = buf[9 * i + 4] >> 4;
        t |= (uint32_t)buf[9 * i + 5] << 4;
        t |= (uint32_t)buf[9 * i + 6] << 12;
        t &= 0x3FFFF;
        r[4*i + 2] += GAMMA1 - t;

        t  = buf[9 * i + 6] >> 6;
        t |= (uint32_t)buf[9 * i + 7] << 2;
        t |= (uint32_t)buf[9 * i + 8] << 10;
        t &= 0x3FFFF;
        r[4*i + 3] += GAMMA1 - t;
    }

    return 4*i;
}*/

/*************************************************
* Name:        poly_uniform_gamma1m1
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-(GAMMA1 - 1), GAMMA1] by unpacking output stream
*              of SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length CRHBYTES
*              - uint16_t nonce: 16-bit nonce
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_uniform_gamma1(poly *a,
        const uint8_t seed[CRHBYTES],
        uint16_t nonce,
        uint8_t buf[sizeof(stream256_state) + 16]) {
    unsigned int ctr, off, k;
    unsigned int buflen = STREAM256_BLOCKBYTES;
    stream256_state *state = (stream256_state *) (buf + 16);
    buf = buf + 7;

    stream256_init(state, seed, nonce);
    stream256_squeezeblocks(buf + 9, 1, state);

    ctr = unpack_gamma1(a->coeffs, N, buf + 9, buflen);

    while (ctr < N) {
        off = buflen % 9;
        for (k = 9 - off; k < 9; ++k) {
            buf[k] = buf[STREAM256_BLOCKBYTES + k];
        }

        stream256_squeezeblocks(buf + 9, 1, state);
        buflen = STREAM256_BLOCKBYTES + off;
        ctr += unpack_gamma1(a->coeffs + ctr, N - ctr, buf + 9 - off, buflen);
    }
}


/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_stream_add_uniform_gamma1(poly *a,
        const uint8_t seed[CRHBYTES],
        uint16_t nonce,
        uint8_t buf[sizeof(stream256_state) + 16]) {
    
    unsigned int ctr, off, k;
    unsigned int buflen = STREAM256_BLOCKBYTES;
    stream256_state *state = (stream256_state *) (buf + 16);
    buf = buf + 7;

    stream256_init(state, seed, nonce);
    stream256_squeezeblocks(buf + 9, 1, state);

    ctr = unpack_gamma1_and_add(a->coeffs, N, buf + 9, buflen);

    while (ctr < N) {
        off = buflen % 9;
        for (k = 9 - off; k < 9; ++k) {
            buf[k] = buf[STREAM256_BLOCKBYTES + k];
        }

        stream256_squeezeblocks(buf + 9, 1, state);
        buflen = STREAM256_BLOCKBYTES + off;
        ctr += unpack_gamma1_and_add(a->coeffs + ctr, N - ctr, buf + 9 - off, buflen);
    }

}*/
/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_challenge
*
* Description: Implementation of H. Samples polynomial with TAU nonzero
*              coefficients in {-1,1} using the output stream of
*              SHAKE256(seed).
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const uint8_t mu[]: byte array containing seed of length SEEDBYTES
**************************************************/
/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_challenge(poly *c, const uint8_t seed[SEEDBYTES], uint8_t buf[sizeof(shake256incctx) + SHAKE256_RATE]) {
    unsigned int i, b, pos;
    uint64_t signs;
    shake256incctx *state = (shake256incctx *) buf;
    buf = buf + sizeof(shake256incctx);

    shake256_inc_init(state);
    shake256_inc_absorb(state, seed, SEEDBYTES);
    shake256_inc_finalize(state);
    shake256_inc_squeeze(buf, SHAKE256_RATE, state);

    signs = 0;
    for (i = 0; i < 8; ++i) {
        signs |= (uint64_t)buf[i] << 8 * i;
    }
    pos = 8;

    for (i = 0; i < N; ++i) {
        c->coeffs[i] = 0;
    }
    for (i = N - TAU; i < N; ++i) {
        do {
            if (pos >= SHAKE256_RATE) {
                shake256_inc_squeeze(buf, SHAKE256_RATE, state);
                pos = 0;
            }

            b = buf[pos++];
        } while (b > i);

        c->coeffs[i] = c->coeffs[b];
        c->coeffs[b] = 1 - 2 * (signs & 1);
        signs >>= 1;
    }
}*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compact_challenge(chall_poly *cp, const uint8_t seed[SEEDBYTES], 
        uint8_t buf[POLY_BYTES], uint8_t hash_buf[sizeof(shake256incctx)]) {
    unsigned int i, b, pos;
    poly *c = (poly*) buf;
    /* All bytes are always used up so we can use the same address for the ouput buffer as for the state */
    shake256incctx *state = (shake256incctx *) hash_buf;

    shake256_inc_init(state);
    shake256_inc_absorb(state, seed, SEEDBYTES);
    shake256_inc_finalize(state);
    shake256_inc_squeeze(hash_buf, SHAKE256_RATE, state);

    cp->signs = 0;
    for (i = 0; i < 8; ++i) {
        cp->signs |= (uint64_t)hash_buf[i] << 8 * i;
    }
    pos = 8;

    for (i = 0; i < N; ++i) {
        c->coeffs[i] = 0;
    }
    for (i = N - TAU; i < N; ++i) {
        do {
            if (pos >= SHAKE256_RATE) {
                shake256_inc_squeeze(hash_buf, SHAKE256_RATE, state);
                pos = 0;
            }

            b = hash_buf[pos++];
        } while (b > i);

        c->coeffs[i] = c->coeffs[b];
        c->coeffs[b] = 1 - 2 * (cp->signs & 1);
        cp->signs >>= 1;
    }
    
    pos = 0;
    cp->signs = 0;
    for (i = 0; i < N; ++i){
        if (c->coeffs[i] != 0) {
            cp->signs |= ((uint64_t)((c->coeffs[i] >> 1) & 0x1) << pos);
            cp->index[pos] = i;
            pos++;
        }
    }
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_CLEAN_polyeta_pack
*
* Description: Bit-pack polynomial with coefficients in [-ETA,ETA].
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYETA_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyeta_pack(uint8_t *r, const poly *a) {
    unsigned int i;
    uint8_t t[8];
    DBENCH_START();

    for (i = 0; i < N / 8; ++i) {
        t[0] = (uint8_t) (ETA - a->coeffs[8 * i + 0]);
        t[1] = (uint8_t) (ETA - a->coeffs[8 * i + 1]);
        t[2] = (uint8_t) (ETA - a->coeffs[8 * i + 2]);
        t[3] = (uint8_t) (ETA - a->coeffs[8 * i + 3]);
        t[4] = (uint8_t) (ETA - a->coeffs[8 * i + 4]);
        t[5] = (uint8_t) (ETA - a->coeffs[8 * i + 5]);
        t[6] = (uint8_t) (ETA - a->coeffs[8 * i + 6]);
        t[7] = (uint8_t) (ETA - a->coeffs[8 * i + 7]);

        r[3 * i + 0]  = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
        r[3 * i + 1]  = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
        r[3 * i + 2]  = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
    }

    DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyeta_unpack
*
* Description: Unpack polynomial with coefficients in [-ETA,ETA].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyeta_unpack(poly *r, const uint8_t *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N / 8; ++i) {
        r->coeffs[8 * i + 0] =  (a[3 * i + 0] >> 0) & 7;
        r->coeffs[8 * i + 1] =  (a[3 * i + 0] >> 3) & 7;
        r->coeffs[8 * i + 2] = ((a[3 * i + 0] >> 6) | (a[3 * i + 1] << 2)) & 7;
        r->coeffs[8 * i + 3] =  (a[3 * i + 1] >> 1) & 7;
        r->coeffs[8 * i + 4] =  (a[3 * i + 1] >> 4) & 7;
        r->coeffs[8 * i + 5] = ((a[3 * i + 1] >> 7) | (a[3 * i + 2] << 1)) & 7;
        r->coeffs[8 * i + 6] =  (a[3 * i + 2] >> 2) & 7;
        r->coeffs[8 * i + 7] =  (a[3 * i + 2] >> 5) & 7;

        r->coeffs[8 * i + 0] = ETA - r->coeffs[8 * i + 0];
        r->coeffs[8 * i + 1] = ETA - r->coeffs[8 * i + 1];
        r->coeffs[8 * i + 2] = ETA - r->coeffs[8 * i + 2];
        r->coeffs[8 * i + 3] = ETA - r->coeffs[8 * i + 3];
        r->coeffs[8 * i + 4] = ETA - r->coeffs[8 * i + 4];
        r->coeffs[8 * i + 5] = ETA - r->coeffs[8 * i + 5];
        r->coeffs[8 * i + 6] = ETA - r->coeffs[8 * i + 6];
        r->coeffs[8 * i + 7] = ETA - r->coeffs[8 * i + 7];
    }

    DBENCH_STOP(*tpack);
}*/

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt1_pack
*
* Description: Bit-pack polynomial t1 with coefficients fitting in 10 bits.
*              Input coefficients are assumed to be standard representatives.
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYT1_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt1_pack(uint8_t *r, const poly *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N / 4; ++i) {
        r[5 * i + 0] = (uint8_t) (a->coeffs[4 * i + 0] >> 0);
        r[5 * i + 1] = (uint8_t) ((a->coeffs[4 * i + 0] >> 8) | (a->coeffs[4 * i + 1] << 2));
        r[5 * i + 2] = (uint8_t) ((a->coeffs[4 * i + 1] >> 6) | (a->coeffs[4 * i + 2] << 4));
        r[5 * i + 3] = (uint8_t) ((a->coeffs[4 * i + 2] >> 4) | (a->coeffs[4 * i + 3] << 6));
        r[5 * i + 4] = (uint8_t) (a->coeffs[4 * i + 3] >> 2);
    }

    DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt1_unpack
*
* Description: Unpack polynomial t1 with 10-bit coefficients.
*              Output coefficients are standard representatives.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt1_unpack(poly *r, const uint8_t *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N / 4; ++i) {
        r->coeffs[4 * i + 0] = ((a[5 * i + 0] >> 0) | ((uint32_t)a[5 * i + 1] << 8)) & 0x3FF;
        r->coeffs[4 * i + 1] = ((a[5 * i + 1] >> 2) | ((uint32_t)a[5 * i + 2] << 6)) & 0x3FF;
        r->coeffs[4 * i + 2] = ((a[5 * i + 2] >> 4) | ((uint32_t)a[5 * i + 3] << 4)) & 0x3FF;
        r->coeffs[4 * i + 3] = ((a[5 * i + 3] >> 6) | ((uint32_t)a[5 * i + 4] << 2)) & 0x3FF;
    }

    DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt0_pack
*
* Description: Bit-pack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYT0_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt0_pack(uint8_t *r, const poly *a) {
    unsigned int i;
    uint32_t t[8];
    DBENCH_START();

    for (i = 0; i < N / 8; ++i) {
        t[0] = (1 << (D - 1)) - a->coeffs[8 * i + 0];
        t[1] = (1 << (D - 1)) - a->coeffs[8 * i + 1];
        t[2] = (1 << (D - 1)) - a->coeffs[8 * i + 2];
        t[3] = (1 << (D - 1)) - a->coeffs[8 * i + 3];
        t[4] = (1 << (D - 1)) - a->coeffs[8 * i + 4];
        t[5] = (1 << (D - 1)) - a->coeffs[8 * i + 5];
        t[6] = (1 << (D - 1)) - a->coeffs[8 * i + 6];
        t[7] = (1 << (D - 1)) - a->coeffs[8 * i + 7];

        r[13 * i + 0]  =  (uint8_t) t[0];
        r[13 * i + 1]  =  (uint8_t) (t[0] >>  8);
        r[13 * i + 1] |=  (uint8_t) (t[1] <<  5);
        r[13 * i + 2]  =  (uint8_t) (t[1] >>  3);
        r[13 * i + 3]  =  (uint8_t) (t[1] >> 11);
        r[13 * i + 3] |=  (uint8_t) (t[2] <<  2);
        r[13 * i + 4]  =  (uint8_t) (t[2] >>  6);
        r[13 * i + 4] |=  (uint8_t) (t[3] <<  7);
        r[13 * i + 5]  =  (uint8_t) (t[3] >>  1);
        r[13 * i + 6]  =  (uint8_t) (t[3] >>  9);
        r[13 * i + 6] |=  (uint8_t) (t[4] <<  4);
        r[13 * i + 7]  =  (uint8_t) (t[4] >>  4);
        r[13 * i + 8]  =  (uint8_t) (t[4] >> 12);
        r[13 * i + 8] |=  (uint8_t) (t[5] <<  1);
        r[13 * i + 9]  =  (uint8_t) (t[5] >>  7);
        r[13 * i + 9] |=  (uint8_t) (t[6] <<  6);
        r[13 * i + 10]  =  (uint8_t) (t[6] >>  2);
        r[13 * i + 11]  =  (uint8_t) (t[6] >> 10);
        r[13 * i + 11] |=  (uint8_t) (t[7] <<  3);
        r[13 * i + 12]  =  (uint8_t) (t[7] >>  5);
    }

    DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt0_unpack
*
* Description: Unpack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt0_unpack(poly *r, const uint8_t *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N / 8; ++i) {
        r->coeffs[8 * i + 0]  = a[13 * i + 0];
        r->coeffs[8 * i + 0] |= (uint32_t)a[13 * i + 1] << 8;
        r->coeffs[8 * i + 0] &= 0x1FFF;

        r->coeffs[8 * i + 1]  = a[13 * i + 1] >> 5;
        r->coeffs[8 * i + 1] |= (uint32_t)a[13 * i + 2] << 3;
        r->coeffs[8 * i + 1] |= (uint32_t)a[13 * i + 3] << 11;
        r->coeffs[8 * i + 1] &= 0x1FFF;

        r->coeffs[8 * i + 2]  = a[13 * i + 3] >> 2;
        r->coeffs[8 * i + 2] |= (uint32_t)a[13 * i + 4] << 6;
        r->coeffs[8 * i + 2] &= 0x1FFF;

        r->coeffs[8 * i + 3]  = a[13 * i + 4] >> 7;
        r->coeffs[8 * i + 3] |= (uint32_t)a[13 * i + 5] << 1;
        r->coeffs[8 * i + 3] |= (uint32_t)a[13 * i + 6] << 9;
        r->coeffs[8 * i + 3] &= 0x1FFF;

        r->coeffs[8 * i + 4]  = a[13 * i + 6] >> 4;
        r->coeffs[8 * i + 4] |= (uint32_t)a[13 * i + 7] << 4;
        r->coeffs[8 * i + 4] |= (uint32_t)a[13 * i + 8] << 12;
        r->coeffs[8 * i + 4] &= 0x1FFF;

        r->coeffs[8 * i + 5]  = a[13 * i + 8] >> 1;
        r->coeffs[8 * i + 5] |= (uint32_t)a[13 * i + 9] << 7;
        r->coeffs[8 * i + 5] &= 0x1FFF;

        r->coeffs[8 * i + 6]  = a[13 * i + 9] >> 6;
        r->coeffs[8 * i + 6] |= (uint32_t)a[13 * i + 10] << 2;
        r->coeffs[8 * i + 6] |= (uint32_t)a[13 * i + 11] << 10;
        r->coeffs[8 * i + 6] &= 0x1FFF;

        r->coeffs[8 * i + 7]  = a[13 * i + 11] >> 3;
        r->coeffs[8 * i + 7] |= (uint32_t)a[13 * i + 12] << 5;
        r->coeffs[8 * i + 7] &= 0x1FFF;

        r->coeffs[8 * i + 0] = (1 << (D - 1)) - r->coeffs[8 * i + 0];
        r->coeffs[8 * i + 1] = (1 << (D - 1)) - r->coeffs[8 * i + 1];
        r->coeffs[8 * i + 2] = (1 << (D - 1)) - r->coeffs[8 * i + 2];
        r->coeffs[8 * i + 3] = (1 << (D - 1)) - r->coeffs[8 * i + 3];
        r->coeffs[8 * i + 4] = (1 << (D - 1)) - r->coeffs[8 * i + 4];
        r->coeffs[8 * i + 5] = (1 << (D - 1)) - r->coeffs[8 * i + 5];
        r->coeffs[8 * i + 6] = (1 << (D - 1)) - r->coeffs[8 * i + 6];
        r->coeffs[8 * i + 7] = (1 << (D - 1)) - r->coeffs[8 * i + 7];
    }

    DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyz_pack
*
* Description: Bit-pack polynomial with coefficients
*              in [-(GAMMA1 - 1), GAMMA1].
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYZ_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyz_pack(uint8_t *r, const poly *a) {
    unsigned int i;
    uint32_t t[4];
    DBENCH_START();

    for (i = 0; i < N / 4; ++i) {
        t[0] = GAMMA1 - a->coeffs[4 * i + 0];
        t[1] = GAMMA1 - a->coeffs[4 * i + 1];
        t[2] = GAMMA1 - a->coeffs[4 * i + 2];
        t[3] = GAMMA1 - a->coeffs[4 * i + 3];

        r[9 * i + 0]  = (uint8_t) t[0];
        r[9 * i + 1]  = (uint8_t) (t[0] >> 8);
        r[9 * i + 2]  = (uint8_t) (t[0] >> 16);
        r[9 * i + 2] |= (uint8_t) (t[1] << 2);
        r[9 * i + 3]  = (uint8_t) (t[1] >> 6);
        r[9 * i + 4]  = (uint8_t) (t[1] >> 14);
        r[9 * i + 4] |= (uint8_t) (t[2] << 4);
        r[9 * i + 5]  = (uint8_t) (t[2] >> 4);
        r[9 * i + 6]  = (uint8_t) (t[2] >> 12);
        r[9 * i + 6] |= (uint8_t) (t[3] << 6);
        r[9 * i + 7]  = (uint8_t) (t[3] >> 2);
        r[9 * i + 8]  = (uint8_t) (t[3] >> 10);
    }

    DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyz_unpack
*
* Description: Unpack polynomial z with coefficients
*              in [-(GAMMA1 - 1), GAMMA1].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyz_unpack(poly *r, const uint8_t *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N / 4; ++i) {
        r->coeffs[4 * i + 0]  = a[9 * i + 0];
        r->coeffs[4 * i + 0] |= (uint32_t)a[9 * i + 1] << 8;
        r->coeffs[4 * i + 0] |= (uint32_t)a[9 * i + 2] << 16;
        r->coeffs[4 * i + 0] &= 0x3FFFF;

        r->coeffs[4 * i + 1]  = a[9 * i + 2] >> 2;
        r->coeffs[4 * i + 1] |= (uint32_t)a[9 * i + 3] << 6;
        r->coeffs[4 * i + 1] |= (uint32_t)a[9 * i + 4] << 14;
        r->coeffs[4 * i + 1] &= 0x3FFFF;

        r->coeffs[4 * i + 2]  = a[9 * i + 4] >> 4;
        r->coeffs[4 * i + 2] |= (uint32_t)a[9 * i + 5] << 4;
        r->coeffs[4 * i + 2] |= (uint32_t)a[9 * i + 6] << 12;
        r->coeffs[4 * i + 2] &= 0x3FFFF;

        r->coeffs[4 * i + 3]  = a[9 * i + 6] >> 6;
        r->coeffs[4 * i + 3] |= (uint32_t)a[9 * i + 7] << 2;
        r->coeffs[4 * i + 3] |= (uint32_t)a[9 * i + 8] << 10;
        r->coeffs[4 * i + 3] &= 0x3FFFF;

        r->coeffs[4 * i + 0] = GAMMA1 - r->coeffs[4 * i + 0];
        r->coeffs[4 * i + 1] = GAMMA1 - r->coeffs[4 * i + 1];
        r->coeffs[4 * i + 2] = GAMMA1 - r->coeffs[4 * i + 2];
        r->coeffs[4 * i + 3] = GAMMA1 - r->coeffs[4 * i + 3];
    }

    DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyw1_pack
*
* Description: Bit-pack polynomial w1 with coefficients in [0,15] or [0,43].
*              Input coefficients are assumed to be standard representatives.
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYW1_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyw1_pack(uint8_t *r, const poly *a) {
    unsigned int i;
    DBENCH_START();

    for (i = 0; i < N / 4; ++i) {
        r[3 * i + 0]  = (uint8_t) a->coeffs[4 * i + 0];
        r[3 * i + 0] |= (uint8_t) (a->coeffs[4 * i + 1] << 6);
        r[3 * i + 1]  = (uint8_t) (a->coeffs[4 * i + 1] >> 2);
        r[3 * i + 1] |= (uint8_t) (a->coeffs[4 * i + 2] << 4);
        r[3 * i + 2]  = (uint8_t) (a->coeffs[4 * i + 2] >> 4);
        r[3 * i + 2] |= (uint8_t) (a->coeffs[4 * i + 3] << 2);
    }

    DBENCH_STOP(*tpack);
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress_w_pack_w1(comp_w0_poly *b, uint8_t *r, const poly *a) {
    unsigned int i;
    int32_t w1, w0;

    for (i = 0; i < N / 4; ++i) {
        w1 = PQCLEAN_DILITHIUM2_MASK_OPTMEM_decompose(&w0, a->coeffs[4*i + 0]);
        compress_w(&b->coeffs[10*i + 0], w1, w0, 0);
        r[3 * i + 0]  = (uint8_t) w1;

        w1 = PQCLEAN_DILITHIUM2_MASK_OPTMEM_decompose(&w0, a->coeffs[4*i + 1]);
        compress_w(&b->coeffs[10*i + 0], w1, w0, 1);
        r[3 * i + 0] |= (uint8_t) (w1 << 6);
        r[3 * i + 1]  = (uint8_t) (w1 >> 2);

        w1 = PQCLEAN_DILITHIUM2_MASK_OPTMEM_decompose(&w0, a->coeffs[4*i + 2]);
        compress_w(&b->coeffs[10*i + 5], w1, w0, 0);
        r[3 * i + 1] |= (uint8_t) (w1 << 4);
        r[3 * i + 2]  = (uint8_t) (w1 >> 4);

        w1 = PQCLEAN_DILITHIUM2_MASK_OPTMEM_decompose(&w0, a->coeffs[4*i + 3]);
        compress_w(&b->coeffs[10*i + 5], w1, w0, 1);
        r[3 * i + 2] |= (uint8_t) (w1 << 2);
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress(comp_poly *b, const poly *a){
    unsigned int i;

    for (i = 0; i < N; ++i) {
        compress_q(&b->coeffs[3*i], a->coeffs[i]);
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompress(poly *b, const comp_poly *a){
    unsigned int i;

    for (i = 0; i < N; ++i) {
        b->coeffs[i] = decompress_q(&a->coeffs[3*i]);
    }
}


void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress_w0(comp_w0_poly *b, const poly *a){
    unsigned int i;

    for (i = 0; i < N; ++i) {
        compress_w0(&b->coeffs[5*(i/2)], a->coeffs[i], (i & 0x1));
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompress_chall(poly *p, chall_poly *c) {
    unsigned int i;

    for (i = 0; i < N; ++i) {
        p->coeffs[i] = 0;
    }

    for (i = 0; i < TAU; ++i) {
        p->coeffs[c->index[i]] = 1 - 2 * ((c->signs >> i) & 0x1);
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_add_cpoly(poly *c, const poly *a, const comp_poly *b){
    unsigned int i;
    int32_t t;

    for (i = 0; i < N; ++i) {
        t = decompress_q(&b->coeffs[3*i]);
        c->coeffs[i] = a->coeffs[i] + t; 
    }
}

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_sub_from_cpoly(poly *c, const comp_poly *a, const poly *b){
    unsigned int i;
    int32_t t;

    for (i = 0; i < N; ++i) {
        t = decompress_q(&a->coeffs[3*i]);
        c->coeffs[i] = t - b->coeffs[i]; 
    }
}*/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_bitslice2dense(poly *v, uint32_t *u) {
    unsigned i, j;
    uint32_t buff[N];
    
    for (i = 0; i < N/32; ++i) {
        for (j = 0; j < COEF_257_NBITS; ++j) {
            v->coeffs[i * 32 + j] = u[i * COEF_257_NBITS * NSHARES + j];
        }
        for (j = COEF_257_NBITS; j < 32; ++j) {
            v->coeffs[i * 32 + j] = 0;
        }

        transpose32((uint32_t *) &v->coeffs[i * 32]);
    }

    for (i = 0; i < N; ++i) {
        buff[i] = v->coeffs[(i/64)*64 + 32*(i%2) + (i % 64)/2];
    }

    for (i = 0; i < N; ++i) {
        v->coeffs[i] = buff[i];
    }
}
