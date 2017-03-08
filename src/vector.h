// Created at 2016-11-22

#ifndef POCKETKALDI_VECTOR_H_
#define POCKETKALDI_VECTOR_H_

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "alloc.h"

#define PKVECTOR_INITIALSIZE 16;

#define PKVECTOR_DEFINE(T)                                                     \
    typedef struct pk_vector_##T##_t {                                         \
      int size;                                                                \
      int capability;                                                          \
      T *data;                                                                 \
      pk_alloc_t *alloc;                                                       \
    } pk_vector_##T##_t;                                                       \
                                                                               \
    static inline void pk_vector_##T##_init(                                   \
        pk_vector_##T##_t *self,                                               \
        pk_alloc_t *alloc) {                                                   \
      self->alloc = alloc;                                                     \
      self->capability = PKVECTOR_INITIALSIZE;                                 \
      self->size = 0;                                                          \
      self->data = (T *)pk_alloc(self->alloc, self->capability * sizeof(T));   \
    }                                                                          \
                                                                               \
    static inline void pk_vector_##T##_push_back(                              \
        pk_vector_##T##_t *self,                                               \
        T item) {                                                              \
      if (self->size >= self->capability) {                                    \
        int new_capability = self->capability * 2;                             \
        void *new_ptr = pk_realloc(                                            \
            self->alloc,                                                       \
            self->data,                                                        \
            new_capability * sizeof(T));                                       \
        self->capability = new_capability;                                     \
        self->data = (T *)new_ptr;                                             \
      }                                                                        \
                                                                               \
      self->data[self->size] = item;                                           \
      ++(self->size);                                                          \
    }                                                                          \
                                                                               \
    static inline void pk_vector_##T##_destroy(pk_vector_##T##_t *self) {      \
      pk_free(self->alloc, self->data);                                        \
      self->size = 0;                                                          \
      self->capability = 0;                                                    \
      self->data = NULL;                                                       \
      self->alloc = NULL;                                                      \
    }                                                                          \

#endif  // POCKETKALDI_VECTOR_H_