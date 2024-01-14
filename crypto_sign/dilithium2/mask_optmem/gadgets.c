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
#include "gadgets.h"
#include "masked.h"
#include "masked_representations.h"
#include "masked_utils.h"
#include <stdint.h>

static inline uint32_t pini_and_core(uint32_t a, uint32_t b, uint32_t r) {
  uint32_t temp;
  uint32_t s;
  /*__asm__("eor %[temp], %[b], %[r]\n\t"
      "and %[temp], %[a], %[temp]\n\t"
      "bic %[s], %[r], %[a]\n\t"
      "eor %[s], %[s], %[temp]"
      : [ s ] "=r"(s),
        [ temp ] "=&r"(
            temp) outputs, use temp as an arbitrary-location clobber
      : [ a ] "r"(a), [ b ] "r"(b), [ r ] "r"(r) inputs
  );*/
  temp = b ^ r;
  temp = a & temp;
  s = r & ~a;
  s = s ^ temp;
  return s;
}

/*************************************************
 * Name:        masked_and
 *
 * Description: Performs masked AND (z = a & b ) gate with nshares.
 *
 * Arguments:   - size_t nshares: number of shares
 *            - uint32_t *z: output buffer
 *            - size_t z_stride: output buffer stride
 *            - uint32_t *a: first input buffer
 *            - size_t a_stride: a buffer stride
 *            - uint32_t *b: second input buffer
 *            - size_t b_stride: b buffer stride
 **************************************************/
void masked_and_c(size_t nshares, uint32_t *z, size_t z_stride,
                  const uint32_t *a, size_t a_stride, const uint32_t *b,
                  size_t b_stride) {
  uint32_t ztmp[nshares];
  uint32_t r;
  uint32_t i, j;

  for (i = 0; i < nshares; ++i) {
    ztmp[i] = a[i * a_stride] & b[i * b_stride];
  }

  for (i = 0; i < (nshares - 1); ++i) {
    for (j = i + 1; j < nshares; ++j) {
      r = rand32();
      // PINI
      ztmp[i] ^= pini_and_core(a[i * a_stride], b[j * b_stride], r);
      ztmp[j] ^= pini_and_core(a[j * a_stride], b[i * b_stride], r);
    }
  }
  for (i = 0; i < nshares; ++i) {
    z[i * z_stride] = ztmp[i];
  }
}

/*************************************************
 * Name:        masked_xor
 *
 * Description: Performs masked XOR (z = a ^ b ) gate with nshares.
 *
 * Arguments:   - size_t nshares: number of shares
 *            - uint32_t *z: output buffer
 *            - size_t z_stride: output buffer stride
 *            - uint32_t *a: first input buffer
 *            - size_t a_stride: a buffer stride
 *            - uint32_t *b: second input buffer
 *            - size_t b_stride: b buffer stride
 **************************************************/
void masked_xor_c(size_t nshares, uint32_t *out, size_t out_stride,
                  const uint32_t *ina, size_t ina_stride, const uint32_t *inb,
                  size_t inb_stride) {
  for (size_t i = 0; i < nshares; ++i) {
    out[i * out_stride] = ina[i * ina_stride] ^ inb[i * inb_stride];
  }
}

/*************************************************
 * Name:        copy_sharing
 *
 * Description: Copy input sharing to output sharing
 *
 * Arguments: - size_t nshares: number of shares
 *            - uint32_t *out: output buffer
 *            - size_t out_stride: out buffer stride
 *            - uint32_t *in: input buffer
 *            - size_t in_stride: in buffer stride
 **************************************************/
void copy_sharing_c(size_t nshares, uint32_t *out, size_t out_stride,
                  const uint32_t *in, size_t in_stride) {

  // TODO generate the ASM too
  for (size_t i = 0; i < nshares; ++i) {
    out[i * out_stride] = in[i * in_stride];
  }
}

/*************************************************
 * Name:        secadd
 *
 * Description: Performs masked addition of two bitslice words.
 *              out = (in1 + in2)%(2**kbits_out)
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words
 *            - size_t kbits_out: number of bits in the output word.
 *                kbits = kbits_out or kbits = kbits_out - 1
 *            - uint32_t *out: output buffer
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - uint32_t *in1: first input buffer
 *            - size_t in1_msk_stride: stride between shares
 *            - size_t in1_data_stride: stride between masked bits
 *            - uint32_t *in2: second input buffer
 *            - size_t in2_msk_stride: stride between shares
 *            - size_t in2_data_stride: stride between masked bits
 **************************************************/
void secadd(size_t nshares, size_t kbits, size_t in_kbits, size_t kbits_out, uint32_t *out,
            size_t out_msk_stride, size_t out_data_stride, const uint32_t *in1,
            size_t in1_msk_stride, size_t in1_data_stride, const uint32_t *in2,
            size_t in2_msk_stride, size_t in2_data_stride, uint32_t buf[4 * nshares]) {

  size_t i;
  uint32_t *carry, *xpy, *xpc, *dummy;

  carry = buf;
  xpy = buf + nshares;
  xpc = buf + 2*nshares;
  dummy = buf + 3*nshares;

  for (i = 0; i < nshares; ++i){
    dummy[i] = 0;
  }

  masked_and(nshares, carry, 1, &in1[0 * in1_data_stride], in1_msk_stride,
             &in2[0 * in2_data_stride], in2_msk_stride);

  masked_xor(nshares, &out[0 * out_data_stride], out_msk_stride,
             &in1[0 * in1_data_stride], in1_msk_stride,
             &in2[0 * in2_data_stride], in2_msk_stride);

  for (i = 1; i < in_kbits; ++i) {
    // xpy = in2 ^ in1
    // xpc = in1 ^ carry
    // out = xpy ^ carry

    masked_xor(nshares, xpy, 1, &in1[i * in1_data_stride], in1_msk_stride,
               &in2[i * in2_data_stride], in2_msk_stride);
    masked_xor(nshares, xpc, 1, &in1[i * in1_data_stride], in1_msk_stride,
               carry, 1);
    masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, xpy, 1,
               carry, 1);

    if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
      break;
    } else if (i == (kbits - 1)) {
      masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
      masked_xor(nshares, &out[(kbits)*out_data_stride], out_msk_stride, carry,
                 1, &in1[i * in1_data_stride], in1_msk_stride);
      break;
    }

    masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
    masked_xor(nshares, carry, 1, carry, 1, &in1[i * in1_data_stride],
               in1_msk_stride);
  }

  for (i = in_kbits; i < kbits; ++i) {
    // xpy = in2 ^ in1
    // xpc = in1 ^ carry
    // out = xpy ^ carry

    masked_xor(nshares, xpy, 1, dummy, 1,
               &in2[i * in2_data_stride], in2_msk_stride);
    masked_xor(nshares, xpc, 1, dummy, 1,
               carry, 1);
    masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, xpy, 1,
               carry, 1);

    if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
      break;
    } else if (i == (kbits - 1)) {
      masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
      masked_xor(nshares, &out[(kbits)*out_data_stride], out_msk_stride, carry,
                 1, dummy, 1);
      break;
    }

    masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
    masked_xor(nshares, carry, 1, carry, 1, dummy, 1);
  }
}

/*************************************************
 * Name:        secadd_constant_bmsk
 *
 * Description: Performs masked addition of an input with a masked
 *              constant sucht that:
 *              out = (in + bmsk*constant)%(2**kbits_out)
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words
 *            - size_t kbits_out: number of bits in the output word.
 *                kbits = kbits_out or kbits = kbits_out - 1
 *            - uint32_t *out: output buffer
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - uint32_t *in: first input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 *            - uint32_t constant: public constant
 *            - uint32_t *bmsk: masked bit buffer
 *            - size_t bmsk_msk_stride: shares stride of masked bit
 **************************************************/
void secadd_constant_bmsk(size_t nshares, size_t kbits, size_t kbits_out,
                          uint32_t *out, size_t out_msk_stride,
                          size_t out_data_stride, const uint32_t *in1,
                          size_t in1_msk_stride, size_t in1_data_stride,
                          uint32_t constant, const uint32_t *bmsk,
                          size_t bmsk_msk_stride, uint32_t buf[4 * nshares]) {
  size_t i, d;
  uint32_t *carry, *xpy, *xpc, *xpz;

  carry = buf;
  xpy = buf + nshares;
  xpc = buf + 2*nshares;
  xpz = buf + 3*nshares;

  if (constant & 0x1) {
    masked_and(nshares, carry, 1, &in1[0 * in1_data_stride], in1_msk_stride,
               bmsk, bmsk_msk_stride);
    masked_xor(nshares, &out[0 * out_data_stride], out_msk_stride,
               &in1[0 * in1_data_stride], in1_msk_stride, bmsk,
               bmsk_msk_stride);
  } else {
    for (d = 0; d < nshares; ++d) {
      carry[d] = 0;
    }
    copy_sharing(nshares, out, out_msk_stride, in1, in1_msk_stride);
  }
  for (i = 1; i < kbits; ++i) {
    // xpy = in2 ^ in1
    // xpc = in1 ^ carry
    // out = xpy ^ carry
    copy_sharing(nshares, xpz, 1, &in1[i * in1_data_stride], in1_msk_stride);
    if ((constant >> i) & 0x1) {
      masked_xor(nshares, xpy, 1, xpz, 1, bmsk, bmsk_msk_stride);
      masked_xor(nshares, xpc, 1, xpz, 1, carry, 1);
      masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, xpy, 1,
                 carry, 1);

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
        masked_xor(nshares, &out[(kbits)*out_data_stride], out_msk_stride,
                   carry, 1, xpz, 1);
        return;
      }

      masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
      masked_xor(nshares, carry, 1, carry, 1, xpz, 1);
    } else {
      // compute the carry
      masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, carry, 1, xpz, 1);

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, &out[(kbits)*out_data_stride], out_msk_stride, carry, 1, xpz, 1);
        return;
      }
      masked_and(nshares, carry, 1, carry, 1, xpz, 1);
    }
  }
}

/*************************************************
 * Name:        secadd_constant
 *
 * Description: Performs masked addition of an input with a
 *              public constant such that
 *              out = (in + constant)%(2**kbits_out)
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words
 *            - size_t kbits_out: number of bits in the output word.
 *                kbits = kbits_out or kbits = kbits_out - 1
 *            - uint32_t *out: output buffer
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - uint32_t *in: first input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 *            - uint32_t constant: public constant
 **************************************************/
void secadd_constant(size_t nshares, size_t kbits, size_t kbits_out,
                     uint32_t *out, size_t out_msk_stride,
                     size_t out_data_stride, const uint32_t *in1,
                     size_t in1_msk_stride, size_t in1_data_stride,
                     uint32_t constant, uint32_t buf[3 * nshares]) {

  size_t i, d;
  uint32_t *carry, *xpy, *xpc;

  carry = buf;
  xpy = buf + nshares;
  xpc = buf + 2*nshares;

  if (constant & 0x1) {
    for (d = 0; d < nshares; ++d) {
      carry[d] = in1[d * in1_msk_stride];
    }
    copy_sharing(nshares, out, out_msk_stride, in1, in1_msk_stride);
    out[0] ^= 0xFFFFFFFF;
  } else {
    for (d = 0; d < nshares; ++d) {
      carry[d] = 0;
    }
    copy_sharing(nshares, out, out_msk_stride, in1, in1_msk_stride);
  }

  for (i = 1; i < kbits; ++i) {
    copy_sharing(nshares, xpy, 1, &in1[i * in1_data_stride], in1_msk_stride);

    if ((constant >> i) & 0x1) {
      masked_xor(nshares, xpc, 1, xpy, 1, carry, 1);
      masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, xpy, 1,
                 carry, 1);

      xpy[0] ^= 0xFFFFFFFF;
      out[i * out_data_stride] ^= 0xFFFFFFFF;

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
        xpy[0] ^= 0xFFFFFFFF;
        masked_xor(nshares, &out[(kbits)*out_data_stride], out_msk_stride,
                   carry, 1, xpy, 1);
        xpy[0] ^= 0xFFFFFFFF;

        // add the kbits_out of the constant
        out[(kbits)*out_data_stride] ^=
            0xFFFFFFFF * ((constant >> kbits) & 0x1);

        return;
      }

      masked_and(nshares, carry, 1, xpy, 1, xpc, 1);

      xpy[0] ^= 0xFFFFFFFF;
      
      masked_xor(nshares, carry, 1, carry, 1, xpy, 1);
    } else {
      // compute the carry
      masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, carry, 1, xpy, 1);

      if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
        return;
      } else if (i == (kbits - 1)) {
        masked_and(nshares, &out[(kbits)*out_data_stride], out_msk_stride,
                   carry, 1, xpy, 1);

        // add the kbits_out of the constant
        out[(kbits)*out_data_stride] ^=
            0xFFFFFFFF * ((constant >> kbits) & 0x1);
        return;
      }
      masked_and(nshares, carry, 1, carry, 1, xpy, 1);
    }
  }
}

/*************************************************
 * Name:        secadd_modp
 *
 * Description: Performs masked addition of two bitslice words with
 *              arbitrary modulo such that
 *              out = (in1 + in2)%(p)
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words:
 *                requires kbits = ceil(log(p))
 *            - uint32_t: modulo for the reduction
 *            - uint32_t *out: output buffer
 *            - size_t out_msk_stride: stride between shares
 *            - size_t out_data_stride: stride between masked bits
 *            - uint32_t *in1: first input buffer
 *            - size_t in1_msk_stride: stride between shares
 *            - size_t in1_data_stride: stride between masked bits
 *            - uint32_t *in2: second input buffer
 *            - size_t in2_msk_stride: stride between shares
 *            - size_t in2_data_stride: stride between masked bits
 **************************************************/
void secadd_modp(size_t nshares, size_t kbits, size_t in_kbits, uint32_t q, uint32_t *out,
                 size_t out_msk_stride, size_t out_data_stride,
                 const uint32_t *in1, size_t in1_msk_stride,
                 size_t in1_data_stride, const uint32_t *in2,
                 size_t in2_msk_stride, size_t in2_data_stride, uint32_t buf[4*nshares]) {

  secadd(nshares, kbits + 1, in_kbits, kbits + 1, out, out_msk_stride, out_data_stride, in1, in1_msk_stride,
         in1_data_stride, in2, in2_msk_stride, in2_data_stride, buf);

  secadd_constant(nshares, kbits + 1, kbits + 1, out, out_msk_stride, out_data_stride, out, out_msk_stride, out_data_stride,
                  (1 << (kbits + 1)) - q, buf);

  secadd_constant_bmsk(nshares, kbits, kbits, out, out_msk_stride,
                       out_data_stride, out, out_msk_stride, out_data_stride, q, &out[kbits * out_data_stride],
                       out_msk_stride, buf);
}

static void masked_splited_and(size_t nshares, uint32_t *z, size_t z_stride,
                  const uint32_t *a, size_t a_stride) {
  uint32_t ztmp[nshares];
  uint32_t r;
  uint32_t i, j;

  for (i = 0; i < nshares; ++i) {
    ztmp[i] = 0;
  }

  for (i = 0; i < nshares/2; ++i) {
    for (j = nshares/2; j < nshares; ++j) {
      r = rand32();
      // PINI
      ztmp[i] ^= pini_and_core(a[i * a_stride], a[j * a_stride], r);
      ztmp[j] ^= r;
    }
  }
  for (i = 0; i < nshares; ++i) {
    z[i * z_stride] = ztmp[i];
  }
}

static void masked_half_xor(size_t nshares, uint32_t *out, size_t out_stride,
                  const uint32_t *ina, size_t ina_stride, const uint32_t *inb,
                  size_t inb_stride) {
  unsigned int i;

  for (i = 0; i < nshares/2; ++i) {
    out[i * out_stride] = ina[i * ina_stride] ^ inb[i * inb_stride];
  }

  for (i = nshares/2; i < nshares; ++i) {
    out[i * out_stride] = inb[i * inb_stride];
  }
}

void secadd_splitted(size_t nshares, size_t kbits, size_t kbits_out, uint32_t *out,
            size_t out_msk_stride, size_t out_data_stride, const uint32_t *in1,
            size_t in1_msk_stride, size_t in1_data_stride, uint32_t buf[3 * nshares]) {

  size_t i;
  uint32_t *carry, *xpy, *xpc;

  carry = buf;
  xpy = buf + nshares;
  xpc = buf + 2*nshares;

  masked_splited_and(nshares, carry, 1, &in1[0 * in1_data_stride], in1_msk_stride);

  copy_sharing(nshares, out, out_msk_stride, in1, in1_msk_stride);

  for (i = 1; i < kbits; ++i) {
    // xpy = in2 ^ in1
    // xpc = in1 ^ carry
    // out = xpy ^ carry

    copy_sharing(nshares, xpy, 1, &in1[i * in1_data_stride], in1_msk_stride);

    masked_half_xor(nshares, xpc, 1, xpy, 1, carry, 1);
    masked_xor(nshares, &out[i * out_data_stride], out_msk_stride, xpy, 1,
               carry, 1);

    if ((i == (kbits - 1)) && (i == (kbits_out - 1))) {
      break;
    } else if (i == (kbits - 1)) {
      masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
      masked_half_xor(nshares, &out[(kbits)*out_data_stride], out_msk_stride, xpy, 1, carry, 1);
      break;
    }

    masked_and(nshares, carry, 1, xpy, 1, xpc, 1);
    masked_half_xor(nshares, carry, 1, xpy, 1, carry, 1);
  }

}

/*************************************************
 * Name:        seca2b_modp
 *
 * Description: Inplace arithmetic to boolean masking conversion:
 *            sum(in_i)%(p) = XOR(in_i)
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words. kbits =
 *ceil(log(p))
 *            - uint32_t p: modulus of the arithmetic masking
 *            - uint32_t *in: input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 **************************************************/
void seca2b_modp(size_t nshares, size_t kbits, uint32_t p, uint32_t *in,
                 size_t in_msk_stride, size_t in_data_stride, uint32_t buf[4 * nshares]) {

  size_t d;

  if (nshares == 1) {
    return;
  }

  size_t nshares_low = nshares / 2;
  size_t nshares_high = nshares - nshares_low;

  seca2b_modp(nshares_low, kbits, p, in, in_msk_stride, in_data_stride, buf);
  seca2b_modp(nshares_high, kbits, p, &in[nshares_low * in_msk_stride],
              in_msk_stride, in_data_stride, buf);

  secadd_constant(nshares_low, kbits, kbits + 1, in, in_msk_stride, in_data_stride, in,
                  in_msk_stride, in_data_stride, (1 << (kbits + 1)) - p, buf);

  for (d = nshares_low; d < nshares; ++d) {
    in[kbits * in_data_stride + d * in_msk_stride] = 0;
  }

  secadd_splitted(nshares, kbits + 1, kbits + 1, in, in_msk_stride, in_data_stride, in, in_msk_stride, in_data_stride, buf);

  secadd_constant_bmsk(nshares, kbits, kbits, in, in_msk_stride, in_data_stride,
                       in, in_msk_stride, in_data_stride, p, &in[kbits * in_data_stride], in_msk_stride, buf);
}

/*************************************************
 * Name:        secb2a_modp
 *
 * Description: Inplace boolean to arithmetic masking conversion:
 *            sum(in_i)%(p) = XOR(in_i)
 *
 * /!\ This function is operating on two slices such that 64 conversions are
 * performed in parallel.
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t kbits: number of bits in input words. kbits =
 *ceil(log(p))
 *            - uint32_t p: modulus of the arithmetic masking
 *            - uint32_t *in: input buffer
 *            - size_t in_msk_stride: stride between shares
 *            - size_t in_data_stride: stride between masked bits
 **************************************************/
void secb2a_modp_centered_in(size_t nshares,
                size_t in_kbits, uint32_t p, 
                uint32_t *in, size_t in_msk_stride, size_t in_data_stride, 
                uint32_t *out, size_t out_msk_stride, size_t out_data_stride,
                uint32_t constant, const uint32_t *bmsk, size_t bmsk_msk_stride, 
                uint32_t buf[COEF_NBITS * (nshares - 1) + (COEF_NBITS + 1) * nshares + 4 * nshares]) {

  uint32_t *z_dense, *b_str, *add_buf;
  size_t d, i;

  z_dense = buf;
  b_str = z_dense + COEF_NBITS * (nshares - 1);
  add_buf = b_str + (COEF_NBITS + 1) * nshares;

  // generate uniform sharing for z
  for (d = 0; d < nshares - 1; ++d) {
    for (i = 0; i < BSSIZE; ++i) {
      z_dense[d * BSSIZE + i] = rand_q();
    }
  }

  // map z to bitslice representation
  for (d = 0; d < nshares - 1; ++d) {
    transpose32(&z_dense[d * BSSIZE]);
  }
  // compress z
  for (d = 0; d < nshares - 1; ++d){
    for (i = 0; i < COEF_NBITS; ++i) {
      z_dense[d * COEF_NBITS + i] = z_dense[d * BSSIZE + i];
    }
  }

  // copy z to b
  for (i = 0; i < COEF_NBITS; ++i) {
    for (d = 0; d < nshares - 1; ++d){
      b_str[d * (COEF_NBITS + 1) + i] = z_dense[d * COEF_NBITS + i];
    }
  }

  // last shares of b set to zero
  for (i = 0; i < COEF_NBITS; ++i) {
    b_str[(nshares - 1) * (COEF_NBITS + 1) + i] = 0;
  }

  seca2b_modp(nshares, COEF_NBITS, p, b_str, COEF_NBITS + 1, 1, add_buf);

  secadd_constant_bmsk(nshares, COEF_NBITS, COEF_NBITS + 1, b_str, COEF_NBITS + 1, 1, 
              b_str, COEF_NBITS + 1, 1, constant, bmsk, bmsk_msk_stride, add_buf);
  
  secadd_modp(nshares, COEF_NBITS, in_kbits, p, b_str, COEF_NBITS + 1, 1, in, in_msk_stride,
              in_data_stride, b_str, COEF_NBITS + 1, 1, add_buf);

  /*secadd(nshares, COEF_NBITS + 1, in_kbits, COEF_NBITS + 1, b_str, COEF_NBITS + 1, 1, in, in_msk_stride,
         in_data_stride, b_str, COEF_NBITS + 1, 1, add_buf);

  secadd_constant(nshares, COEF_NBITS + 1, COEF_NBITS + 1, b_str, COEF_NBITS + 1, 1, b_str, COEF_NBITS + 1, 1,
                  (1 << (COEF_NBITS + 1)) - p, add_buf);

  secadd_constant_bmsk(nshares, COEF_NBITS, COEF_NBITS, b_str, COEF_NBITS + 1,
                       1, b_str, COEF_NBITS + 1, 1, p, &b_str[COEF_NBITS], COEF_NBITS + 1, add_buf);*/

  // do p - z for each share and coefficient
  // one's complement
  for (i = 0; i < COEF_NBITS * (nshares - 1); ++i){
    z_dense[i] = ~z_dense[i];
  }

  // add p + 1
  for (d = 0; d < nshares - 1; ++d) {
    secadd_constant(1, COEF_NBITS, COEF_NBITS, &z_dense[d * COEF_NBITS], COEF_NBITS, 1, &z_dense[d * COEF_NBITS], COEF_NBITS, 1, p + 1, add_buf);
  }

  // copy z_dense to out 
  for (d = 0; d < nshares - 1; ++d){
    for (i = 0; i < COEF_NBITS; ++i) {
      out[i*out_data_stride + d * out_msk_stride] = z_dense[d * COEF_NBITS + i];
    }
  }

  // unmask b_str
  for (i = 0; i < COEF_NBITS; ++i) {
    RefreshIOS_rec(nshares, nshares, &b_str[i], COEF_NBITS + 1);

    out[(nshares - 1) * out_msk_stride + i*out_data_stride] = b_str[i];
    for (d = 1; d < nshares; ++d) {
      out[(nshares - 1) * out_msk_stride + i*out_data_stride] ^= b_str[d * (COEF_NBITS + 1) + i];
    }
  } 
}


void secb2a_mod257(size_t nshares,
                 //   size_t kbits, // MUST BE EQUAL TO COEF_257_NBITS
                 uint32_t p, uint32_t *in, size_t in_msk_stride,
                 size_t in_data_stride, 
                 uint32_t buf[COEF_257_NBITS * (nshares - 1) + (COEF_257_NBITS + 1) * nshares + 4 * nshares]) {

  uint32_t *z_dense, *b_str, *add_buf;
  size_t d, i;

  if (nshares == 1) {
    return;
  }

  z_dense = buf;
  b_str = z_dense + COEF_257_NBITS * (nshares - 1);
  add_buf = b_str + (COEF_257_NBITS + 1) * nshares;

  // generate uniform sharing for z
  for (d = 0; d < nshares - 1; ++d) {
    for (i = 0; i < BSSIZE; ++i) {
      z_dense[(nshares - 2)*COEF_257_NBITS + i] = rand_257();
    }
    transpose32(&z_dense[(nshares - 2)*COEF_257_NBITS]);
    for (i = 0; i < COEF_257_NBITS; ++i) {
      z_dense[d*COEF_257_NBITS + i] = z_dense[(nshares - 2)*COEF_257_NBITS + i];
    }
  }

  // copy z to b
  for (d = 0; d < nshares - 1; ++d){
    for (i = 0; i < COEF_257_NBITS; ++i) {
      b_str[d * (COEF_257_NBITS + 1) + i] = z_dense[d * COEF_257_NBITS + i];
    }
  }

  // last shares of b set to zero
  for (i = 0; i < COEF_257_NBITS; ++i) {
    b_str[(nshares - 1) * (COEF_257_NBITS + 1) + i] = 0;
  }

  seca2b_modp(nshares, COEF_257_NBITS, p, b_str, COEF_257_NBITS + 1, 1, add_buf);

  for (d = 0; d < nshares; ++d) {
    b_str[d * (COEF_257_NBITS + 1) + COEF_257_NBITS] = 0;
  }

  secadd_modp(nshares, COEF_257_NBITS, COEF_257_NBITS, p, b_str, COEF_257_NBITS + 1, 1, in, in_msk_stride,
              in_data_stride, b_str, COEF_257_NBITS + 1, 1, add_buf);

  // do p - z for each share and coefficient
  // one's complement
  for (i = 0; i < COEF_257_NBITS * (nshares - 1); ++i){
    z_dense[i] = ~z_dense[i];
  }

  // add p + 1
  for (d = 0; d < nshares - 1; ++d) {
    secadd_constant(1, COEF_257_NBITS, COEF_257_NBITS, &z_dense[d * COEF_257_NBITS], COEF_257_NBITS, 1, &z_dense[d * COEF_257_NBITS], COEF_257_NBITS, 1, p + 1, add_buf);
  }

  // copy z_dense to in 
  for (d = 0; d < nshares - 1; ++d){
    for (i = 0; i < COEF_257_NBITS; ++i) {
      in[i*in_data_stride + d * in_msk_stride] = z_dense[d * COEF_257_NBITS + i];
    }
  }

  // unmask b_str
  for (i = 0; i < COEF_257_NBITS; ++i) {
    RefreshIOS_rec(nshares, nshares, &b_str[i], COEF_257_NBITS + 1);

    in[(nshares - 1) * in_msk_stride + i] = b_str[i];
    for (d = 1; d < nshares; ++d) {
      in[(nshares - 1) * in_msk_stride + i] ^= b_str[d * (COEF_257_NBITS + 1) + i];
    }
  } 
}

void seca2a_centered_modp(size_t nshares,
                 //   size_t kbits, // MUST BE EQUAL TO COEF_NBITS
                 uint32_t in_p, uint32_t *in, size_t in_msk_stride,
                 size_t in_data_stride, uint32_t out_p, uint32_t *out, size_t out_msk_stride,
                 size_t out_data_stride, 
                 uint32_t buf[COEF_NBITS * (nshares - 1) + (COEF_NBITS + 1) * nshares + 4 * nshares]) {
    uint32_t i, j;
    uint32_t bmsk[NSHARES];

    for (i = 0; i < nshares; ++i){
        for (j = 0; j < COEF_257_NBITS; ++j){
            buf[i*(COEF_257_NBITS + 1) + j] = in[i*in_msk_stride + j*in_data_stride];
        }
        buf[i*(COEF_257_NBITS + 1) + COEF_257_NBITS] = 0;
    }

    seca2b_modp(nshares, COEF_257_NBITS, in_p, buf, COEF_257_NBITS + 1, 1, &buf[nshares*(COEF_257_NBITS + 1)]);

    for (i = 0; i < nshares; ++i){
        for (j = 0; j < COEF_257_NBITS; ++j){
            in[i*in_msk_stride + j*in_data_stride] = buf[i*(COEF_257_NBITS + 1) + j] ;
        }
    }
    
    masked_xor(nshares, bmsk, 1, &in[(COEF_257_NBITS-1)*in_data_stride], in_msk_stride, &in[(COEF_257_NBITS-2)*in_data_stride], in_msk_stride);

    secb2a_modp_centered_in(nshares, COEF_257_NBITS, out_p, in, in_msk_stride, in_data_stride,
                          out, out_msk_stride, out_data_stride, Q - QP, bmsk, 1, buf);
}


/*************************************************
 * Name:        RefreshIOS_rec
 *
 * Description: IOS refresh on boolean sharing
 *
 * Arguments: - size_t nshares: number of shares
 *            - size_t d: current recursion shars:
 *            - uint32_t *x: input buffer
 *            - size_t x_msk_stride: stride between shares
 **************************************************/
void RefreshIOS_rec(size_t nshares, size_t d, uint32_t *x,
                    size_t x_msk_stride) {
  uint32_t r;
  if (d == 1) {
  } else if (d == 2) {
    r = rand32();
    x[0 * x_msk_stride] ^= r;
    x[1 * x_msk_stride] ^= r;
  } else {
    RefreshIOS_rec(nshares, d / 2, x, x_msk_stride);
    RefreshIOS_rec(nshares, d - d / 2, &x[(d / 2) * x_msk_stride],
                   x_msk_stride);
    for (unsigned int i = 0; i < d / 2; ++i) {
      r = rand32();
      x[i * x_msk_stride] ^= r;
      x[(i + d / 2) * x_msk_stride] ^= r;
    }
  }
}
