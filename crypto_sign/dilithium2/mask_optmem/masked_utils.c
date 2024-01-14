#include "masked_utils.h"

/*************************************************
 * Name:        randombytes_masking
 *
 * Description: Generates random bytes. Defined to 
 *              be used as another source of randomness and thus work
 *              with tests that use deterministic randomness.
 *
 * Arguments: - size_t n: bytes written to output
 **************************************************/
void randombytes_masking(uint8_t *output, size_t n){
  size_t i;

  for(i = 0; i < n; ++i){
    output[i] = (uint8_t) rand(); 
  }
}

uint32_t rand32() {
  uint32_t r = 0;
  randombytes_masking((uint8_t *)&r, 4);
  return r; 
}

uint32_t rand_q() {
  uint32_t r;
  
  do {
    r = 0;
    randombytes_masking((uint8_t *)&r, 3);
    r &= 0x7FFFFF;
  } while (r >= Q);

  return r;
}

uint32_t rand_257() {
  uint32_t r = 0;
  
  do {
    r = 0;
    randombytes_masking((uint8_t *)&r, 2);
    r &= 0x1FF;
  } while (r >= 257);

  return r;
}