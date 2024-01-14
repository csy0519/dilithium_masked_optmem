#ifndef SHAKE_H
#define SHAKE_H

#include <stddef.h>
#include <stdint.h>

#define SHAKE128_RATE 168
#define SHAKE256_RATE 136


#define PQC_SHAKEINCCTX_BYTES (sizeof(uint64_t)*26)
#define PQC_SHAKECTX_BYTES (sizeof(uint64_t)*25)

// Context for incremental API
typedef struct {
    uint64_t ctx[PQC_SHAKEINCCTX_BYTES/sizeof(uint64_t)];
} shake128incctx;

// Context for non-incremental API
typedef struct {
    uint64_t ctx[PQC_SHAKECTX_BYTES/sizeof(uint64_t)];
} shake128ctx;

// Context for incremental API
typedef struct {
    uint64_t ctx[PQC_SHAKEINCCTX_BYTES/sizeof(uint64_t)];
} shake256incctx;

// Context for non-incremental API
typedef struct {
    uint64_t ctx[PQC_SHAKECTX_BYTES/sizeof(uint64_t)];
} shake256ctx;


/* Initialize the state and absorb the provided input.
 *
 * This function does not support being called multiple times
 * with the same state.
 */
void shake128_absorb(shake128ctx *state, const uint8_t *input, size_t inlen);
/* Squeeze output out of the sponge.
 *
 * Supports being called multiple times
 */
void shake128_squeezeblocks(uint8_t *output, size_t nblocks, shake128ctx *state);
/* Copy the state. */
void shake128_ctx_clone(shake128ctx *dest, const shake128ctx *src);

/* Initialize incremental hashing API */
void shake128_inc_init(shake128incctx *state);
/* Absorb more information into the XOF.
 *
 * Can be called multiple times.
 */
void shake128_inc_absorb(shake128incctx *state, const uint8_t *input, size_t inlen);
/* Finalize the XOF for squeezing */
void shake128_inc_finalize(shake128incctx *state);
/* Squeeze output out of the sponge.
 *
 * Supports being called multiple times
 */
void shake128_inc_squeeze(uint8_t *output, size_t outlen, shake128incctx *state);
/* Copy the context of the SHAKE128 XOF */
void shake128_inc_ctx_clone(shake128incctx* dest, const shake128incctx *src);

/* Initialize the state and absorb the provided input.
 *
 * This function does not support being called multiple times
 * with the same state.
 */
void shake256_absorb(shake256ctx *state, const uint8_t *input, size_t inlen);
/* Squeeze output out of the sponge.
 *
 * Supports being called multiple times
 */
void shake256_squeezeblocks(uint8_t *output, size_t nblocks, shake256ctx *state);
/* Copy the context held by this XOF */
void shake256_ctx_clone(shake256ctx *dest, const shake256ctx *src);

/* Initialize incremental hashing API */
void shake256_inc_init(shake256incctx *state);
void shake256_inc_absorb(shake256incctx *state, const uint8_t *input, size_t inlen);
/* Prepares for squeeze phase */
void shake256_inc_finalize(shake256incctx *state);
/* Squeeze output out of the sponge.
 *
 * Supports being called multiple times
 */
void shake256_inc_squeeze(uint8_t *output, size_t outlen, shake256incctx *state);
/* Copy the state */
void shake256_inc_ctx_clone(shake256incctx* dest, const shake256incctx *src);

/* One-stop SHAKE128 call */
void shake128(uint8_t *output, size_t outlen,
              const uint8_t *input, size_t inlen,
              uint8_t buf[sizeof(shake128ctx)]);

/* One-stop SHAKE256 call */
void shake256(uint8_t *output, size_t outlen,
              const uint8_t *input, size_t inlen, 
              uint8_t buf[sizeof(shake256ctx)]);

#endif
