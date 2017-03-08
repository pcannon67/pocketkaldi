// 2017-01-27

#include "matrix.h"

#include <assert.h>
#include "stdlib.h"
#include "util.h"

void pk_matrix_init(pk_matrix_t *self, int nrow, int ncol) {
  self->nrow = nrow;
  self->ncol = ncol;
  int total_elems = nrow * ncol;
  self->data = (float *)pk_alloc(sizeof(float) * total_elems);
}

void pk_matrix_fill(pk_matrix_t *self, float val) {
  for (int i = 0; i < self->nrow * self->ncol; ++i) {
    self->data[i] = val;
  }
}

void pk_matrx_setcol(pk_matrix_t *self, int col, const float *data) {
  assert(col >= 0 && col < self->ncol && "pk_matrix: index out of boundary");

  int start_idx = col * self->nrow;
  for (int i = 0; i < self->nrow; ++i) {
    self->data[start_idx + i] = data[i];
  }
}

const float *pk_matrix_getcol(const pk_matrix_t *self, int col) {
  assert(col >= 0 && col < self->ncol && "pk_matrix: index out of boundary");

  int start_idx = col * self->nrow;
  return self->data + start_idx;
}
