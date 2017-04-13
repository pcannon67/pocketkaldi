// Created at 2016-11-22

#ifndef POCKETKALDI_LIST_H_
#define POCKETKALDI_LIST_H_

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include "util.h"

#define PKLIST_INITIALSIZE 16;

#define PKLIST_DEFINE(T, name) \
    typedef struct name##_t { \
      int size; \
      T *data; \
      int _capability; \
    } name##_t; \
    /* Initialize the list */ \
    static inline void name##_init(name##_t *self) { \
      self->_capability = PKLIST_INITIALSIZE; \
      self->size = 0; \
      self->data = (T *)pk_alloc(self->_capability * sizeof(T)); \
    } \
    /* Requests that the capacity be at least enough to contain n elements.*/ \
    static inline void name##_reserve(name##_t *self, int cap) { \
        void *new_ptr = realloc(self->data, cap * sizeof(T)); \
        self->_capability = cap; \
        self->data = (T *)new_ptr; \
    } \
    /* Add an empty element into list, and returns the pointer to it */ \
    static inline T *name##_emplace_back( \
        name##_t *self) { \
      if (self->size >= self->_capability) { \
        name##_reserve(self, self->_capability * 2); \
      } \
      return &self->data[self->size++]; \
    } \
    /* Push back an element into list */ \
    static inline void name##_push_back( \
        name##_t *self, \
        T item) { \
      T *back = name##_emplace_back(self); \
      *back = item; \
    } \
    /* Clear all elements in the list */ \
    static inline void name##_clear(name##_t *self) { \
      self->size = 0; \
    } \
    /* Destroy the list */ \
    static inline void name##_destroy(name##_t *self) { \
      pk_free(self->data); \
      self->size = 0; \
      self->_capability = 0; \
      self->data = NULL; \
    } \
    /* Returns the last element */ \
    static inline T name##_back(name##_t *self) { \
      assert(self->size > 0 && "list is empty when calling _back()"); \
      return self->data[self->size - 1]; \
    } \
    /* Pop the last element from list */ \
    static inline void name##_pop_back(name##_t *self) { \
      assert(self->size > 0 && "list is empty when calling pop_back()"); \
      --self->size; \
    } \
    /* Returns true if the list is empty */ \
    static inline bool name##_empty(name##_t *self) { \
      return self->size == 0; \
    }

#endif  // POCKETKALDI_LIST_H_
