#ifndef PQCLEAN_DILITHIUM2_MASK_OPTMEM_POLYVEC_H
#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_POLYVEC_H
#include "params.h"
#include "poly.h"
#include <stdint.h>

/* Vectors of polynomials of length L */
typedef struct {
    poly vec[L];
} polyvecl;


/* Vectors of polynomials of length K */
typedef struct {
    poly vec[K];
} polyveck;

typedef struct {
    comp_poly vec[K];
} comp_polyveck;

typedef struct {
    comp_w0_poly vec[K];
} comp_w0_polyveck;

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_cpolyveck_init(comp_polyveck *v, const uint32_t val);

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyvec_stream_acc_matrix_pointwise(comp_polyveck *v, 
        const poly *u, 
        const uint8_t rho[SEEDBYTES], 
        const uint32_t i, 
        uint8_t buf[sizeof(stream128_state) + 2]);

#endif
