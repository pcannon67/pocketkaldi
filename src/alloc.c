// Created at 2016-11-08

#include <stdlib.h>
#include "alloc.h"

#define UNUSED(x) (void)(x)

void pk_alloc_init(pk_alloc_t *alloc) {
  UNUSED(alloc);
}

void *pk_alloc(pk_alloc_t *alloc, size_t size) {
  UNUSED(alloc);
  return malloc(size);
}

void *pk_realloc(pk_alloc_t *alloc, void *ptr, size_t size) {
  UNUSED(alloc);
  return realloc(ptr, size);
}

void pk_free(pk_alloc_t *alloc, void *pointer) {
  UNUSED(alloc);
  free(pointer);
}
