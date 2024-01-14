#include "params.h"
#include "reduce.h"
#include <stdint.h>

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce
*
* Description: For finite field element a with -2^{31}Q <= a <= Q*2^31,
*              compute r \equiv a*2^{-32} (mod Q) such that -Q < r < Q.
*
* Arguments:   - int64_t: finite field element a
*
* Returns r.
**************************************************/
int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce(int64_t a) {
    int32_t t;

    t = (int32_t)((uint64_t)a * (uint64_t)QINV);
    t = (a - (int64_t)t * Q) >> 32;
    return t;
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_reduce32
*
* Description: For finite field element a with a <= 2^{31} - 2^{22} - 1,
*              compute r \equiv a (mod Q) such that -6283009 <= r <= 6283007.
*
* Arguments:   - int32_t: finite field element a
*
* Returns r.
**************************************************/
int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_reduce32(int32_t a) {
    int32_t t;

    t = (a + (1 << 22)) >> 23;
    t = a - t * Q;
    return t;
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_caddq
*
* Description: Add Q if input coefficient is negative.
*
* Arguments:   - int32_t: finite field element a
*
* Returns r.
**************************************************/
int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_caddq(int32_t a) {
    a += (a >> 31) & Q;
    return a;
}

/*************************************************
* Name:        PQCLEAN_DILITHIUM2_MASK_OPTMEM_freeze
*
* Description: For finite field element a, compute standard
*              representative r = a mod^+ Q.
*
* Arguments:   - int32_t: finite field element a
*
* Returns r.
**************************************************/
int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_freeze(int32_t a) {
    a = PQCLEAN_DILITHIUM2_MASK_OPTMEM_reduce32(a);
    a = PQCLEAN_DILITHIUM2_MASK_OPTMEM_caddq(a);
    return a;
}


int16_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce_257(int32_t a) {
    int16_t t;

    t = (int16_t)a * QPINV;
    t = (a - (int16_t)t * QP) >> 16;
    return t;
}

int16_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_caddqp(int16_t a) {
    a += (a >> 15) & QP;
    return a;
}
