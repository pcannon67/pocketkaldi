// 2017-01-27

#include "matrix.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

void pk_matrix_read(pk_matrix_t *self, pk_readable_t *fd, pk_status_t *status) {
  pk_bytebuffer_t content;
  pk_vector_t col_vec;
  pk_bytebuffer_init(&content);
  pk_vector_init(&col_vec, 0, NAN);

  // Clear self first
  pk_matrix_destroy(self);

  int32_t section_size = pk_readable_readsectionhead(
    fd,
    PK_MATRIX_SECTION,
    status);
  if (!status->ok) goto pk_matrix_read_failed;
  pk_bytebuffer_reset(&content, section_size);

  // Read columns and rows
  pk_readable_readbuffer(fd, &content, status);
  if (!status->ok) goto pk_matrix_read_failed;
  int num_col = pk_bytebuffer_readint32(&content);
  int num_row = pk_bytebuffer_readint32(&content);

  // Initialize the matrix
  pk_matrix_destroy(self);
  pk_matrix_init(self, num_row, num_col);

  // Read each column vectors
  for (int col_idx = 0; col_idx < num_col; ++col_idx) {
    pk_vector_read(&col_vec, fd, status);
    if (!status->ok) goto pk_matrix_read_failed;
    if (col_vec.dim != num_row) {
      PK_STATUS_CORRUPTED(
          status,
          "pk_matrix_read: col_vec.dim == %d expected, but %d found (%s)",
          num_row,
          col_vec.dim,
          fd->filename);
      goto pk_matrix_read_failed;
    }

    for (int d = 0; d < col_vec.dim; ++d) {
      self->data[col_idx * num_row + d] = col_vec.data[d];
    }
  }

  if (false) {
pk_matrix_read_failed:
    pk_matrix_destroy(self);
  }

  pk_bytebuffer_destroy(&content);
  pk_vector_destroy(&col_vec);
}

void pk_matrix_resize(pk_matrix_t *self, int nrow, int ncol) {
  int self_totalelem = self->nrow * self->ncol;
  int total_elem = nrow * ncol;
  self->nrow = nrow;
  self->ncol = ncol;
  if (self_totalelem == total_elem) {
    return;
  } else if (total_elem > 0) {
    self->data = (float *)pk_realloc(self->data, sizeof(float) * total_elem);
  } else {
    free(self->data);
    self->data = NULL;
  }
}

void pk_matrix_copy(pk_matrix_t *dest, const pk_matrix_t *src) {
  if (dest->nrow != src->nrow || dest->ncol != src->ncol) {
    pk_matrix_resize(dest, src->nrow, src->ncol);
  }

  for (int i = 0; i < src->nrow * src->ncol; ++i) {
    dest->data[i] = src->data[i];
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

pk_vector_t pk_matrix_getcol(const pk_matrix_t *self, int col) {
  assert(col >= 0 && col < self->ncol && "pk_matrix: index out of boundary");
  int start_idx = col * self->nrow;

  pk_vector_t column;
  column.dim = self->nrow;
  column.data = self->data + start_idx;
  return column;
}
