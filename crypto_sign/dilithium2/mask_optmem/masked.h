#ifndef MASKED_H
#define MASKED_H
#include "params.h"
#include <stdint.h>

#define BSSIZE 32
#define COEF_NBITS 23
#define COEF_257_NBITS 9

#define ETA_NBITS 3

typedef struct{
    uint32_t coeffs[ETA_NBITS][N/32];
} BitslicePolyEta;

typedef struct{
    BitslicePolyEta BPoly[K];
} BitslicePolyveckEta;

typedef struct{
    BitslicePolyEta BPoly[L];
} BitslicePolyveclEta;

#endif