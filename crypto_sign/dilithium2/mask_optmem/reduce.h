#ifndef PQCLEAN_DILITHIUM2_MASK_OPTMEM_REDUCE_H
#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_REDUCE_H
#include "params.h"
#include <stdint.h>

#define MONT (-4186625) // 2^32 % Q
#define QINV 58728449 // q^(-1) mod 2^32
#define QPINV (-255) // qp^(-1) mod 2^16 

int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce(int64_t a);

int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_reduce32(int32_t a);

int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_caddq(int32_t a);

int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_freeze(int32_t a);

int16_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_montgomery_reduce_257(int32_t a);

int16_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_caddqp(int16_t a);

#endif
