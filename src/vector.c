
#include "vector.h"

#include <math.h>
#include <cblas.h>

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
  pk_vector_resize(out, W->ncol);
  cblas_sgemm(
      CblasColMajor,
      CblasNoTrans,
      CblasNoTrans,
      1,
      W->ncol,
      self->dim,
      1.0,
      self->data,
      1,
      W->data,
      W->nrow,
      0.0,
      out->data,
      1);
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
