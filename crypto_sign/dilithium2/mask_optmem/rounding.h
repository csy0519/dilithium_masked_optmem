#ifndef PQCLEAN_DILITHIUM2_MASK_OPTMEM_ROUNDING_H
#define PQCLEAN_DILITHIUM2_MASK_OPTMEM_ROUNDING_H
#include "params.h"
#include <stdint.h>

int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_power2round(int32_t *a0, int32_t a);

int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_power2round_a1(int32_t a);

int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_power2round_a0(int32_t a1, int32_t a);

int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_decompose(int32_t *a0, int32_t a);

unsigned int PQCLEAN_DILITHIUM2_MASK_OPTMEM_make_hint(int32_t a0, int32_t a1);

int32_t PQCLEAN_DILITHIUM2_MASK_OPTMEM_use_hint(int32_t a, unsigned int hint);

#endif
