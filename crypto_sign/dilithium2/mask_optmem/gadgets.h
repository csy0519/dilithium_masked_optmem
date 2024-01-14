/* Copyright 2022 UCLouvain, Belgium and PQM4 contributors
 *
 * This file is part of pqm4_masked.
 *
 * pqm4_masked is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, version 3.
 *
 * pqm4_masked is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * pqm4_masked. If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef GADGETS_H
#define GADGETS_H

#include <stddef.h>
#include <stdint.h>
#include "masked.h"

// Atomic gadgets
void masked_xor_c(size_t nshares, uint32_t *out, size_t out_stride,
                  const uint32_t *ina, size_t ina_stride, const uint32_t *inb,
                  size_t inb_stride);
void masked_and_c(size_t nshares, uint32_t *z, size_t z_stride,
                  const uint32_t *a, size_t a_stride, const uint32_t *b,
                  size_t b_stride);
void copy_sharing_c(size_t nshares, uint32_t *out, size_t out_stride,
                  const uint32_t *in, size_t in_stride);


#define masked_and(nshares, z, z_stride, a, a_stride, b, b_stride)             \
  masked_and_c(nshares, z, z_stride, a, a_stride, b, b_stride)
#define masked_xor(nshares, z, z_stride, a, a_stride, b, b_stride)             \
  masked_xor_c(nshares, z, z_stride, a, a_stride, b, b_stride)
#define copy_sharing(nshares, out, out_stride, in, in_stride)                   \
  copy_sharing_c(nshares,out,out_stride,in, in_stride)

void RefreshIOS_rec(size_t nshares, size_t d, uint32_t *x, size_t x_msk_stride);
// Adders
void secadd(size_t nshares, size_t kbits, size_t in_kbits, size_t kbits_out, uint32_t *out,
            size_t out_msk_stride, size_t out_data_stride, const uint32_t *in1,
            size_t in1_msk_stride, size_t in1_data_stride, const uint32_t *in2,
            size_t in2_msk_stride, size_t in2_data_stride, uint32_t buf[3 * nshares]);
void secadd_modp(size_t nshares, size_t kbits, size_t in_kbits, uint32_t p, uint32_t *out,
                 size_t out_msk_stride, size_t out_data_stride,
                 const uint32_t *in1, size_t in1_msk_stride,
                 size_t in1_data_stride, const uint32_t *in2,
                 size_t in2_msk_stride, size_t in2_data_stride, uint32_t buf[4 * nshares]);
                 
void secadd_constant_bmsk(size_t nshares, size_t kbits, size_t kbits_out,
                          uint32_t *out, size_t out_msk_stride,
                          size_t out_data_stride, const uint32_t *in1,
                          size_t in1_msk_stride, size_t in1_data_stride,
                          uint32_t constant, const uint32_t *bmsk,
                          size_t bmsk_msk_stride, uint32_t buf[4 * nshares]);
void secadd_constant(size_t nshares, size_t kbits, size_t kbits_out,
                     uint32_t *out, size_t out_msk_stride,
                     size_t out_data_stride, const uint32_t *in1,
                     size_t in1_msk_stride, size_t in1_data_stride,
                     uint32_t constant, uint32_t buf[3 * nshares]);

void secadd_splitted(size_t nshares, size_t kbits, size_t kbits_out, uint32_t *out,
            size_t out_msk_stride, size_t out_data_stride, const uint32_t *in1,
            size_t in1_msk_stride, size_t in1_data_stride, uint32_t buf[3 * nshares]);

// Conversions
void seca2b_modp(size_t nshares, size_t kbits, uint32_t p, uint32_t *in,
                 size_t in_msk_stride, size_t in_data_stride, uint32_t buf[4 * nshares]);
void secb2a_modp_centered_in(size_t nshares,
                size_t in_kbits, uint32_t p, 
                uint32_t *in, size_t in_msk_stride, size_t in_data_stride, 
                uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
                uint32_t constant, const uint32_t *bmsk, size_t bmsk_msk_stride, 
                uint32_t buf[COEF_NBITS * (nshares - 1) + (COEF_NBITS + 1) * nshares + 4 * nshares]);

void secb2a_mod257(size_t nshares, uint32_t p, uint32_t *in, size_t in_msk_stride,
                 size_t in_data_stride,
                 uint32_t buf[COEF_257_NBITS * (nshares - 1) + (COEF_257_NBITS + 1) * nshares + 4 * nshares +  4]);
                 
void seca2a_centered_modp(size_t nshares,
                 //   size_t kbits, // MUST BE EQUAL TO COEF_NBITS
                 uint32_t in_p, uint32_t *in, size_t in_msk_stride,
                 size_t in_data_stride, uint32_t out_p, uint32_t *out, size_t out_msk_stride,
                 size_t out_data_stride, 
                 uint32_t buf[COEF_NBITS * (nshares - 1) + (COEF_NBITS + 1) * nshares + 4 * nshares]);

#endif // GADGETS_H
