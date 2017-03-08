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

float pk_vector_dot(const pk_vector_t *self, const pk_vector_t *vec) {
  assert(self->dim == vec->dim && "pk_vector_dot: vector size mismatch");

  float sum = 0.0;
  for (int i = 0; i < self->dim; ++i) {
    sum += self->data[i] * vec->data[i];
  }
  return sum;
}

void pk_vector_add(const pk_vector_t *self, const pk_vector_t *vec) {
  assert(self->dim == vec->dim && "pk_vector_add: vector size mismatch");

  for (int i = 0; i < self->dim; ++i) {
    self->data[i] += vec->data[i];
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
  pk_bytebuffer_init(&content, section_size);
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
