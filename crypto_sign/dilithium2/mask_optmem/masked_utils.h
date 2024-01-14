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
#ifndef MASKED_UTILS_H
#define MASKED_UTILS_H
#include <stdint.h>
#include <stdlib.h>
#include "poly.h"
#include "params.h"
#include "randombytes.h"

void randombytes_masking(uint8_t *output, size_t n);

uint32_t rand32(void);

uint32_t rand_q(void);

uint32_t rand_257(void);

#endif
