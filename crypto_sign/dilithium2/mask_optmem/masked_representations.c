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
#include "masked_representations.h"
#include "masked.h"
#include "params.h"
#include <stdint.h>


/*************************************************
 * Name: transpose32
 *
 * Description: Transpose 32 x 32 bit matrix
 *
 * Arguments:
 * - uint32_t[32] a: the matrix, updated such that the new value of (a[i] >> j)
 *   & 0x1 is the initial value of (a[j] >> i) & 0x1.
 *
 * Adapted from "Hacker's Delight" by Henry S. Warren, Jr.
 * at http://www.icodeguru.com/Embedded/Hacker's-Delight/
 * **************************************************/
void transpose32(uint32_t a[32]) {
  int j;
  unsigned int k;
  uint32_t m, t;
  m = 0x0000FFFF;
  for (j = 16; j != 0; j = j >> 1, m = m ^ (m << j)) {
    for (k = 0; k < 32; k = (k + j + 1) & ~j) {
      t = (a[k + j] ^ (a[k] >> j)) & m;
      a[k + j] = a[k + j] ^ t;
      a[k] = a[k] ^ (t << j);
    }
  }
}

void polyeta_canonical2bitslice_opt(uint32_t out[32], const poly *p){
  unsigned int i;

  for (i = 0; i < N / 8; ++i) {
      out[i] = (p->coeffs[2*i + 0] + ETA);
      out[i] |= (p->coeffs[2*i + 1] + ETA) << 3;
      out[i] |= (p->coeffs[2*i + 64] + ETA) << 6;
      out[i] |= (p->coeffs[2*i + 65] + ETA) << 9;
      out[i] |= (p->coeffs[2*i + 128] + ETA) << 12;
      out[i] |= (p->coeffs[2*i + 129] + ETA) << 15;
      out[i] |= (p->coeffs[2*i + 192] + ETA) << 18;
      out[i] |= (p->coeffs[2*i + 193] + ETA) << 21;
  }

  transpose32(out);
}

/*void polyeta_bitslice2canonical_opt(poly *r, const BitslicePolyEta *in, uint32_t buf[32]){
  unsigned int i, j;
  uint32_t *a = buf;

  for (i = 0; i < ETA_NBITS; ++i){
    for (j = 0; j < (N/32); ++j){
      a[j*ETA_NBITS + i] = in->coeffs[i][j];
    }
  } 

  for (i = ETA_NBITS*(N/32); i < 32; ++i) {
    a[i] = 0;
  }

  transpose32(a); 

  for (i = 0; i < N / 8; ++i) {
      r->coeffs[i + 0] = ((a[i] >> 0) & 0x7) - ETA;
      r->coeffs[i + 32] = ((a[i] >> 3) & 0x7) - ETA;
      r->coeffs[i + 64] = ((a[i] >> 6) & 0x7) - ETA;
      r->coeffs[i + 96] = ((a[i] >> 9) & 0x7) - ETA;
      r->coeffs[i + 128] = ((a[i] >> 12) & 0x7) - ETA;
      r->coeffs[i + 160] = ((a[i] >> 15) & 0x7) - ETA;
      r->coeffs[i + 192] = ((a[i] >> 18) & 0x7) - ETA;
      r->coeffs[i + 224] = ((a[i] >> 21) & 0x7) - ETA;
  }
}




void masked_dense2bitslice_opt(
    size_t nshares, size_t coeffs_size, uint32_t *bitslice,
    size_t bitslice_msk_stride, size_t bitslice_data_stride,
    const uint32_t *dense, size_t dense_msk_stride, size_t dense_data_stride,
    uint32_t buf[32]) {
  uint32_t *a = buf;
  for (size_t d = 0; d < nshares; ++d) {
    for (size_t i = 0; i < 32; ++i) {
      a[i] = dense[i * dense_data_stride + d * dense_msk_stride];
    }
    transpose32(a);
    for (size_t i = 0; i < coeffs_size; ++i) {
      bitslice[d * bitslice_msk_stride + i * bitslice_data_stride] = a[i];
    }
  }
}

void masked_bitslice2dense_opt(size_t nshares, size_t coeffs_size,
                               uint32_t *dense, size_t dense_msk_stride,
                               size_t dense_data_stride,
                               const uint32_t *bitslice,
                               size_t bitslice_msk_stride,
                               size_t bitslice_data_stride,
                               uint32_t buf[32]) {
  uint32_t *a = buf;
  for (size_t d = 0; d < nshares; ++d) {
    for (size_t i = 0; i < coeffs_size; ++i) {
      a[i] = bitslice[d * bitslice_msk_stride + i * bitslice_data_stride];
    }
    // Avoid uninitialized vars -> UB :(
    for (size_t i = coeffs_size; i < 32; ++i) {
      a[i] = 0;
    }
    transpose32(a);
    for (size_t i = 0; i < 32; ++i) {
      dense[d * dense_msk_stride + i * dense_data_stride] = a[i];
    }
  }
}*/


void masked_inplace_bitslice2dense_opt(size_t nshares, size_t coeffs_size,
                                  uint32_t *bitslice){

  for (int d = (int) nshares - 1; d >= 0; --d) {
    for (int i = (int) coeffs_size - 1; i >= 0; --i) {
      bitslice[d*32 + i] = bitslice[d * coeffs_size + i];
    }
    // Avoid uninitialized vars -> UB :(
    for (int i = (int) coeffs_size; i < 32; ++i) {
      bitslice[d*32 + i] = 0;
    }
    transpose32(&bitslice[d*32]);
  }
}