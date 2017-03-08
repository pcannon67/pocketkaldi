// Created at 2016-11-08

#include <stdint.h>
#include <stdlib.h>
#include "pocketkaldi.h"

#ifndef POCKETKALDI_ALLOC_H_
#define POCKETKALDI_ALLOC_H_

// Tempolary it is just a wrapper of malloc & free
typedef int pk_alloc_t;

// Initialize the allocator
POCKETKALDI_EXPORT
void pk_alloc_init(pk_alloc_t *alloc);

// Allocate a memory block contains at least `size` bytes
POCKETKALDI_EXPORT
void *pk_alloc(pk_alloc_t *, size_t size);

// Remove a memory block pointed by `pointer`
POCKETKALDI_EXPORT
void pk_free(pk_alloc_t *, void *pointer);

#endif  // POCKETKALDI_ALLOC_H_
