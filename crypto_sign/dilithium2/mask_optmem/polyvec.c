#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include <stdint.h>

/**************************************************************/
/************ Vectors of polynomials of length K **************/
/**************************************************************/

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_cpolyveck_init(comp_polyveck *v, const uint32_t val) {
    unsigned int i;

    for (i = 0; i < K; ++i) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_cpoly_init(&v->vec[i], val);
    }
}

void PQCLEAN_DILITHIUM2_MASK_OPTMEM_polyvec_stream_acc_matrix_pointwise(comp_polyveck *v, 
        const poly *u, 
        const uint8_t rho[SEEDBYTES], 
        const uint32_t i, 
        uint8_t buf[sizeof(stream128_state) + 2]){
    unsigned int j;

    for (j = 0; j < K; ++j) {
        PQCLEAN_DILITHIUM2_MASK_OPTMEM_poly_stream_acc_matrix_pointwise(&v->vec[j], u, rho, j, i, buf);
    }
}
