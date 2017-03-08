// Created at 2016-11-22

#ifndef POCKETKALDI_LIST_H_
#define POCKETKALDI_LIST_H_

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include "util.h"

#define PKLIST_INITIALSIZE 16;

#define PKLIST_DEFINE(T, name)                                                 \
    typedef struct name##_t {                                                  \
      int size;                                                                \
      T *data;                                                                 \
      int _capability;                                                         \
    } name##_t;                                                                \
                                                                               \
    static inline void name##_init(name##_t *self) {                           \
      self->_capability = PKLIST_INITIALSIZE;                                  \
      self->size = 0;                                                          \
      self->data = (T *)pk_alloc(self->_capability * sizeof(T));               \
    }                                                                          \
                                                                               \
    static inline void name##_push_back(                                       \
        name##_t *self,                                                        \
        T item) {                                                              \
      if (self->size >= self->_capability) {                                   \
        int new_capability = self->_capability * 2;                            \
        void *new_ptr = pk_realloc(                                            \
            self->data,                                                        \
            new_capability * sizeof(T));                                       \
        self->_capability = new_capability;                                    \
        self->data = (T *)new_ptr;                                             \
      }                                                                        \
                                                                               \
      self->data[self->size] = item;                                           \
      ++(self->size);                                                          \
    }                                                                          \
                                                                               \
    static inline void name##_destroy(name##_t *self) {                        \
      pk_free(self->data);                                                     \
      self->size = 0;                                                          \
      self->_capability = 0;                                                   \
      self->data = NULL;                                                       \
    }                                                                          \
                                                                               \
    static inline T name##_back(name##_t *self) {                              \
      assert(self->size > 0 && "list is empty when calling _back()");          \
      return self->data[self->size - 1];                                       \
    }                                                                          \
                                                                               \
    static inline void name##_pop_back(name##_t *self) {                       \
      assert(self->size > 0 && "list is empty when calling pop_back()");       \
      --self->size;                                                            \
    }                                                                          \
                                                                               \
    static inline bool name##_empty(name##_t *self) {                          \
      return self->size == 0;                                                  \
    }                                                                          

#endif  // POCKETKALDI_LIST_H_