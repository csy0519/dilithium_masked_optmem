#include "masked_smallpoly.h"


static const int16_t zetas_qp[128] = {
    0, 16, 4, 64, 2, 32, 8, 128, -60, 68, 17, 15, -120, -121, 34, 30,
    -35, -46, 117, 73, -70, -92, -23, -111, 44, -67, -81, -11, 88, 123, 95, -22,
    -42, 99, 89, -118, -84, -59, -79, 21, -50, -29, 57, -116, -100, -58, 114, 25,
    -72, -124, -31, 18, 113, 9, -62, 36, -49, -13, 61, -52, -98, -26, 122, -104,
    27, -82, 108, -71, 54, 93, -41, 115, -78, 37, -55, -109, 101, 74, -110, 39,
    83, 43, 75, -85, -91, 86, -107, 87, -97, -10, 126, -40, 63, -20, -5, -80,
    -106, 103, 90, -102, 45, -51, -77, 53, -65, -12, -3, -48, 127, -24, -6, -96,
    112, -7, -66, -28, -33, -14, 125, -56, -38, -94, 105, -119, -76, 69, -47, 19
};


void PQCLEAN_DILITHIUM2_MASK_OPTMEM_fnt_257(int16_t a[N]) {
    unsigned int len, start, j, k;
    int16_t zeta, t;

    k = 0;
    for (len = 128; len > 1; len >>= 1) {
        for (start = 0; start < N; start = j + len) {
            zeta = zetas_qp[++k];
            for (j = start; j < start + len; ++j) {
                t = PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce_257((int32_t)zeta * a[j + len]);
                a[j + len] = a[j] - t;
                a[j] = a[j] + t;
            }
        }
    }
}


void PQCLEAN_DILITHIUM2_MASK_OPTMEM_invfnt_257(int16_t a[N]) {
    unsigned int start, len, j, k;
    int16_t t, zeta;
    const int32_t f = -2; // mont^2/128

    k = 128;
    for (len = 2; len < N/2; len <<= 1) {
        for (start = 0; start < N; start = j + len) {
            zeta = -zetas_qp[--k];
            for (j = start; j < start + len; ++j) {
                t = a[j];
                a[j] = t + a[j + len];
                a[j + len] = t - a[j + len];
                a[j + len] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce_257((int32_t)zeta * a[j + len]);
            }
        }
    }

    /* Last stage merged with multiplication by f to prevent potential overflow*/
    len = 128;
    for (start = 0; start < N; start = j + len) {
        zeta = -zetas_qp[--k];
        for (j = start; j < start + len; ++j) {
            t = a[j];
            a[j] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce_257(f * ((int32_t) t + (int32_t)a[j + len]));
            a[j + len] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce_257(f * (int32_t)zeta * ((int32_t) t - (int32_t) a[j + len]));
        }
    }
}

static int16_t fqmul(int16_t a, int16_t b) {
    return PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce_257((int32_t) (a * b));
}

static void basemul(int16_t *r, int16_t *a, int16_t *b, int16_t zeta) {
    int16_t t[2];

    t[0]  = fqmul(a[1], b[1]);
    t[0]  = fqmul(t[0], zeta);
    t[0] += fqmul(a[0], b[0]);
    t[1]  = fqmul(a[0], b[1]);
    t[1] += fqmul(a[1], b[0]);
    r[0] = t[0];
    r[1] = t[1];
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_basemul(int16_t *r, int16_t *a, int16_t *b) {
    for (int i = 0; i < N/4; ++i){
        basemul(&r[4 * i], &a[4 * i], &b[4 * i], zetas_qp[64 + i]);
        basemul(&r[4 * i + 2], &a[4 * i + 2], &b[4 * i + 2], -zetas_qp[64 + i]);
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_caddqp(int16_t *a){
    for (int i = 0; i < N; ++i){
        a[i] = PQCLEAN_DILITHIUM2_MASK_OPTMEM_caddqp(a[i]);
    }
}


static void load_poly(uint32_t *r, uint32_t *p, size_t mask_stride) {
    for (int i = 0; i < N/BSSIZE; ++i){
        for (int j = 0; j < COEF_257_NBITS; ++j){
            r[16*i+j] = p[i*mask_stride + j];
        }
        for (int j = COEF_257_NBITS; j < 16; ++j){
            r[16*i+j] = 0;
        }
    }
    
    for (int i = 0; i < N/(2*BSSIZE); ++i){
        transpose32(&r[BSSIZE*i]);
    }
}

static void store_poly(uint32_t *r, uint16_t *p, size_t mask_stride, uint16_t buff[2*BSSIZE]) {
    for (int i = 0; i < 4; ++i){
        for (int j = 0; j < 64; ++j) {
            buff[j/32 + 2*(j%32)] = p[64*i + j];
        }
        transpose32((uint32_t *)buff);
        for (int j = 0; j < 2; ++j){
            for (int k = 0; k < 9; ++k){
                r[2*i*mask_stride + j * mask_stride + k] = ((uint32_t *)buff)[16*j+k]; 
            }
        }
    }
}


void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_csi(int16_t* c, uint32_t* p, uint16_t buff[N + 2*BSSIZE]){
    int16_t *curr_poly;
    uint16_t *swap_buff;
    curr_poly = (int16_t *) buff;
    swap_buff = (uint16_t *) (buff + N);

    for (int n = 0; n < NSHARES; ++n){
        
        load_poly((uint32_t *) curr_poly, &p[n*COEF_257_NBITS], COEF_257_NBITS*NSHARES);
        if (n == 0){
            PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_sub_eta((halfpoly *) curr_poly, (halfpoly *) curr_poly);
        }
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_fnt_257(curr_poly);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_basemul(curr_poly, c, curr_poly);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_invfnt_257(curr_poly);
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_caddqp(curr_poly);

        store_poly(&p[n*COEF_257_NBITS], (uint16_t *) curr_poly, NSHARES * COEF_257_NBITS, swap_buff);
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_decompress_chall(halfpoly *c, chall_poly *cp) {
    unsigned int i;
    for (i = 0; i < N; ++i){
        c->coeffs[i] = 0;
    }
    for (i = 0; i < TAU; ++i) {
        c->coeffs[cp->index[i]] = (1 - 2 * ((cp->signs >> i) & 0x1));
    }
}


void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_smallpoly_sub_eta(halfpoly *c, const halfpoly *b){
    unsigned int i;
    for (i = 0; i < N; ++i) {
        c->coeffs[i] = b->coeffs[i] - ETA;
    }
}