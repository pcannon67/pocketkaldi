// 2017-01-27

#include "matrix.h"

#include <assert.h>
#include <stdlib.h>
#include "util.h"

void pk_matrix_init(pk_matrix_t *self, int nrow, int ncol) {
  self->nrow = nrow;
  self->ncol = ncol;
  int total_elems = nrow * ncol;
  if (total_elems > 0) {
    self->data = (float *)pk_alloc(sizeof(float) * total_elems);
  } else {
    self->data = NULL;
  }
}

void pk_matrix_resize(pk_matrix_t *self, int nrow, int ncol) {
  self->nrow = nrow;
  self->ncol = ncol;
  int total_elems = nrow * ncol;
  if (total_elems > 0) {
    self->data = (float *)pk_realloc(self->data, sizeof(float) * total_elems);
  } else {
    free(self->data);
    self->data = NULL;
  }
}

void pk_matrix_destroy(pk_matrix_t *self) {
  self->nrow = 0;
  self->ncol = 0;
  free(self->data);
  self->data = NULL;
}

void pk_matrix_fill(pk_matrix_t *self, float val) {
  for (int i = 0; i < self->nrow * self->ncol; ++i) {
    self->data[i] = val;
  }
}

const pk_vector_t pk_matrix_getcol(const pk_matrix_t *self, int col) {
  assert(col >= 0 && col < self->ncol && "pk_matrix: index out of boundary");
  int start_idx = col * self->nrow;

  pk_vector_t column;
  column.dim = self->nrow;
  column.data = self->data + start_idx;
  return column;
}

void pk_vector_init(pk_vector_t *self, int dim) {
  self->dim = dim;
  if (dim != 0) {
    self->data = (float *)pk_alloc(sizeof(float) * dim);
  } else {
    self->data = NULL;
  }
}

void pk_vector_destroy(pk_vector_t *self) {
  if (self->data != NULL) pk_free(self->data);
  self->data = NULL;
  self->dim = 0;
}

void pk_vector_copyfrom(pk_vector_t *self, const float *source, int n) {
  assert(n <= self->dim && "pk_vector: n should be less than self->dim");
  for (int i = 0; i < n; ++i) {
    self->data[i] = source[i];
  }
}

float pk_vector_dot(const pk_vector_t *self, const pk_vector_t *vec) {
  assert(self->dim == vec->dim && "pk_vector_dot: vector size mismatch");

  float sum = 0.0;
  for (int i = 0; i < self->dim; ++i) {
    sum += self->data[i] * vec->data[i];
  }
  return sum;
}

void pk_vector_resize(pk_vector_t *self, int dim) {
  self->dim = dim;
  if (dim == 0) {
    pk_free(self->data);
    self->data = NULL;
  } else {
    self->data = pk_realloc(self->data, sizeof(float) * dim);
  }
}

void pk_vector_fill(pk_vector_t *self, float value) {
  assert(self->data != NULL && "pk_vector_fill: self-data is NULL");
  for (int i = 0; i < self->dim; ++i) {
    self->data[i] = value;
  }
}

const pk_vector_t 
pk_vector_subvector(const pk_vector_t *self, int start, int dim) {
  assert(dim > 0 && start + dim <= self->dim &&
         "pk_vector_subvector: start and dim out of boundary");
  pk_vector_t subvector;
  subvector.data = self->data + start;
  subvector.dim = dim;
  return subvector;
}
