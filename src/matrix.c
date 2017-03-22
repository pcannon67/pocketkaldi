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

const pk_vector_t pk_matrix_getcol(const pk_matrix_t *self, int col) {
  assert(col >= 0 && col < self->ncol && "pk_matrix: index out of boundary");
  int start_idx = col * self->nrow;

  pk_vector_t column;
  column.dim = self->nrow;
  column.data = self->data + start_idx;
  return column;
}

void pk_vector_init(pk_vector_t *self, int dim, float fill_with) {
  self->dim = dim;
  if (dim != 0) {
    self->data = (float *)pk_alloc(sizeof(float) * dim);
    if (fill_with != NAN) {
      pk_vector_fill(self, fill_with);
    }
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

void pk_vector_copy(pk_vector_t *dest, const pk_vector_t *src) {
  if (dest->dim != src->dim) pk_vector_resize(dest, src->dim);
  
  for (int i = 0; i < src->dim; ++i) {
    dest->data[i] = src->data[i];
  }
}

void pk_vector_floor(pk_vector_t *self, float floor) {
  for (int i = 0; i < self->dim; ++i) {
    if (self->data[i] < floor) self->data[i] = floor;
  }
}

void pk_vector_log(pk_vector_t *self) {
  for (int i = 0; i < self->dim; ++i) {
    self->data[i] = logf(self->data[i]);
  }
}

void pk_vector_scale(pk_vector_t *self, float scale) {
  for (int i = 0; i < self->dim; ++i) {
    self->data[i] = scale * self->data[i];
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

float pk_vector_dotmat(
    const pk_vector_t *self, 
    const pk_matrix_t *W,
    pk_vector_t *out) {
  assert(self->dim == W->nrow && "pk_vector_dotmat: x and W size mismatch");
  
  int nrow = W->nrow;
  int ncol = W->ncol;
  pk_vector_resize(out, ncol);
  for (int c = 0; c < ncol; ++c) {
    float sum = 0.0f;
    for (int d = 0; d < nrow; ++d) {
      sum += self->data[d] * W->data[c * nrow + d];
    }
    out->data[c] = sum;
  }
}

void pk_vector_add(const pk_vector_t *self, const pk_vector_t *vec) {
  assert(self->dim == vec->dim && "pk_vector_add: vector size mismatch");

  for (int i = 0; i < self->dim; ++i) {
    self->data[i] += vec->data[i];
  }
}

void pk_vector_addscale(
    const pk_vector_t *self,
    float scale,
    const pk_vector_t *vec) {
  assert(self->dim == vec->dim && "pk_vector_add: vector size mismatch");

  for (int i = 0; i < self->dim; ++i) {
    self->data[i] += scale * vec->data[i];
  }
}

void pk_vector_scalaradd(
    const pk_vector_t *self,
    float scalar,
    const pk_vector_t *vec) {
  assert(self->dim == vec->dim && "pk_vector_scalaradd: vector size mismatch");

  for (int i = 0; i < self->dim; ++i) {
    self->data[i] += scalar * vec->data[i];
  }
}

void pk_vector_subtract(const pk_vector_t *self, const pk_vector_t *vec) {
  assert(self->dim == vec->dim && "pk_vector_subtract: vector size mismatch");

  for (int i = 0; i < self->dim; ++i) {
    self->data[i] -= vec->data[i];
  }
}

void pk_vector_add2(const pk_vector_t *self, const pk_vector_t *vec) {
  assert(self->dim == vec->dim && "pk_vector_add2: vector size mismatch");

  for (int i = 0; i < self->dim; ++i) {
    float val = vec->data[i];
    self->data[i] += val * val;
  }
}

void pk_vector_subtract2(const pk_vector_t *self, const pk_vector_t *vec) {
  assert(self->dim == vec->dim && "pk_vector_subtract2: vector size mismatch");

  for (int i = 0; i < self->dim; ++i) {
    float val = vec->data[i];
    self->data[i] -= val * val;
  }
}

void pk_vector_resize(pk_vector_t *self, int dim) {
  if (dim == self->dim) {
    // Do nothing when two dimensions are the same
    return;
  } else if (dim == 0) {
    pk_free(self->data);
    self->data = NULL;
    self->dim = 0;
  } else {
    self->data = pk_realloc(self->data, sizeof(float) * dim);
    self->dim = dim;
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

void pk_vector_read(pk_vector_t *self, pk_readable_t *fd, pk_status_t *status) {
  int32_t section_size = pk_readable_readsectionhead(
      fd,
      PK_VECTOR_SECTION,
      status);

  // Read content
  pk_bytebuffer_t content;
  pk_bytebuffer_init(&content);
  pk_bytebuffer_reset(&content, section_size);
  if (status->ok) {
    pk_readable_readbuffer(fd, &content, status);
  }

  // Read and check dimension
  int dim = 0;
  if (status->ok) {
    dim = pk_bytebuffer_readint32(&content);

    // Check dim and content size
    int64_t expected_contentsize = dim * sizeof(float) + sizeof(int32_t);
    if (expected_contentsize != content.size) {
      PK_STATUS_CORRUPTED(
          status,
          "pk_vector_read: content_size == %lld expected, but %lld found (%s)",
          expected_contentsize,
          content.size,
          fd->filename);
    }
  }

  if (status->ok) {
    pk_vector_resize(self, dim);
    for (int i = 0; i < dim; ++i) {
      self->data[i] = pk_bytebuffer_readfloat(&content);
    }
  }

  pk_bytebuffer_destroy(&content);
}
