#ifndef PQCLEAN_DILITHIUM2_MASK_OPTMEM_POLY_H
#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_POLY_H
#include "params.h"
#include "symmetric.h"
#include <stdint.h>

#define POLY_UNIFORM_NBLOCKS ((768 + STREAM128_BLOCKBYTES - 1)/STREAM128_BLOCKBYTES)
#define POLY_UNIFORM_ETA_NBLOCKS ((136 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#define POLY_UNIFORM_GAMMA1_NBLOCKS ((POLYZ_PACKEDBYTES + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)

typedef struct {
    int32_t coeffs[N];
} poly;

typedef struct {
   uint8_t coeffs[3*N];
} comp_poly;

typedef struct {
   uint64_t signs;
   uint8_t index[TAU_MAX]; 
} chall_poly;

typedef struct {
   uint8_t coeffs[5*(N/2)];
} comp_w0_poly;

void compress_q(uint8_t *b, const uint32_t a);

int32_t decompress_q(const uint8_t *a);

void compress_w(uint8_t *b, const int32_t w1, const int32_t w0, const size_t odd);
int32_t decompress_w0(const uint8_t *a, const size_t odd);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_reduce(poly *a);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_caddq(poly *a);
/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_freeze(poly *a);*/

//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_add(poly *c, const poly *a, const poly *b);
//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_sub(poly *c, const poly *a, const poly *b);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_shiftl(poly *a);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_sub_eta(poly *c, const poly *b);

//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_masked_poly_negate(poly *a);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_negate_chall_c(chall_poly *c);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_ntt(poly *a);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_invntt_tomont(poly *a);
//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compressed_pointwise_montgomery(poly *c, const poly *a, const comp_poly *b);

//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compute_ct1(poly *r, const chall_poly *c, const uint8_t *a);
//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compute_ct0(poly *r, const chall_poly *c, const uint8_t *a);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_stream_acc_matrix_pointwise(comp_poly *b, 
        const poly *a, const uint8_t rho[SEEDBYTES], const size_t i, const size_t j, 
        uint8_t buf[sizeof(stream128_state) + 2]);

//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_power2round(poly *a1, poly *a0, const poly *a);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_power2round_a1(poly *a1, const poly *a);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_power2round_a0(poly *a0, const poly *a1, const comp_poly *a);
//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompose(poly *a1, poly *a0, const poly *a);
//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompose_a1(poly *a1, const poly *a);
//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompose_a0(poly *a0, const poly *a1, const comp_poly *a);
//unsigned int PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_make_hint(poly *h, const poly *a0, const poly *a1);
unsigned int PQCLEAN_DILITHIUM2_MASK_OPTMEM_cpoly_make_hint(uint8_t *h, const poly *a, const comp_poly *b);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_cpoly_make_hint_and_pack(uint8_t *h, const poly *a, const comp_w0_poly *b, uint32_t *n, size_t index);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_use_hint(poly *b, const poly *a, const uint8_t *h, const int hlen);

int PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_chknorm(const poly *a, int32_t B);
//int PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_lowbits_chknorm(const poly *a, int32_t B);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_cpoly_init(comp_poly *a, const uint32_t val);
/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_uniform(poly *a,
        const uint8_t seed[SEEDBYTES],
        uint16_t nonce,
        uint8_t buf[sizeof(stream128_state) + POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES + 2]);*/
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_uniform_eta(poly *a,
        const uint8_t seed[CRHBYTES],
        uint16_t nonce,
        uint8_t buf[sizeof(stream256_state) + POLY_UNIFORM_ETA_NBLOCKS * STREAM256_BLOCKBYTES]);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_uniform_gamma1(poly *a,
        const uint8_t seed[CRHBYTES],
        uint16_t nonce,
        uint8_t buf[sizeof(stream256_state) + 16]);

/*void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_stream_add_uniform_gamma1(poly *a,
        const uint8_t seed[CRHBYTES],
        uint16_t nonce,
        uint8_t buf[sizeof(stream256_state) + 16]);*/

//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_challenge(poly *c, const uint8_t seed[SEEDBYTES], uint8_t buf[sizeof(shake256incctx) + SHAKE256_RATE]);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compact_challenge(chall_poly *cp, const uint8_t seed[SEEDBYTES], 
        uint8_t buf[POLY_BYTES], uint8_t hash_buf[sizeof(shake256incctx)]);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyeta_pack(uint8_t *r, const poly *a);
//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyeta_unpack(poly *r, const uint8_t *a);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt1_pack(uint8_t *r, const poly *a);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt1_unpack(poly *r, const uint8_t *a);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt0_pack(uint8_t *r, const poly *a);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyt0_unpack(poly *r, const uint8_t *a);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyz_pack(uint8_t *r, const poly *a);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyz_unpack(poly *r, const uint8_t *a);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyw1_pack(uint8_t *r, const poly *a);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress_w_pack_w1(comp_w0_poly *b, uint8_t *r, const poly *a);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress(comp_poly *b, const poly *a);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompress(poly *b, const comp_poly *a);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_compress_w0(comp_w0_poly *b, const poly *a);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_decompress_chall(poly *p, chall_poly *c);
void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_add_cpoly(poly *c, const poly *a, const comp_poly *b);
//void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_sub_from_cpoly(poly *c, const comp_poly *a, const poly *b);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_bitslice2dense(poly *v, uint32_t *u);

#endif
