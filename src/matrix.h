// 2017-01-27

#ifndef POCKETKALDI_MATRIX_H_
#define POCKETKALDI_MATRIX_H_

#include <stdio.h>
#include "pocketkaldi.h"
#include "util.h"

#define PK_VECTOR_SECTION "VEC0"
#define PK_MATRIX_SECTION "MAT0"

// pk_matrix_t is a column-major matrix type in pocketkaldi
typedef struct pk_matrix_t {
  int ncol;
  int nrow;
  float *data;
} pk_matrix_t;

// pk_vector_t is a simple column vector
typedef struct pk_vector_t {
  int dim;
  float *data;
} pk_vector_t;

// Initialize the matrix with specified rows and columns
POCKETKALDI_EXPORT
void pk_matrix_init(pk_matrix_t *self, int nrow, int ncol);

// Read matrix from fd. When any error occured, it will set status->ok to false
POCKETKALDI_EXPORT
void pk_matrix_read(pk_matrix_t *self, pk_readable_t *fd, pk_status_t *status);

// Fill the matrix with value
POCKETKALDI_EXPORT
void pk_matrix_fill(pk_matrix_t *self, float val);

// Resize the matrix
POCKETKALDI_EXPORT
void pk_matrix_resize(pk_matrix_t *self, int nrow, int ncol);

// Copy the matrix from src to dest
POCKETKALDI_EXPORT
void pk_matrix_copy(pk_matrix_t *dest, const pk_matrix_t *src);

// Destroy the matrix
POCKETKALDI_EXPORT
void pk_matrix_destroy(pk_matrix_t *self);

// Returns a borrowed vector of the specified column. The values (data) could be
// changed. And since it is just borrowed from the matrix, destory is not
// needed. 
POCKETKALDI_EXPORT
const pk_vector_t pk_matrix_getcol(const pk_matrix_t *self, int col);

// Initialize the vector and fill with a number. If fill_with == NAN, doesn't
// fill anything
POCKETKALDI_EXPORT
void pk_vector_init(pk_vector_t *self, int dim, float fill_with);

// Read vector from file. If failed, set status to fail
POCKETKALDI_EXPORT
void pk_vector_read(pk_vector_t *self, pk_readable_t *fd, pk_status_t *status);

// Resize the vector
POCKETKALDI_EXPORT
void pk_vector_resize(pk_vector_t *self, int dim);

// Fill the vector with value
POCKETKALDI_EXPORT
void pk_vector_fill(pk_vector_t *self, float value);

// Copy n values from source to the vector. n should be less than self->dim
POCKETKALDI_EXPORT
void pk_vector_copyfrom(pk_vector_t *self, const float *source, int n);

// Copy vector from src to dest
POCKETKALDI_EXPORT
void pk_vector_copy(pk_vector_t *dest, const pk_vector_t *src);

// Apply floor to the vector
POCKETKALDI_EXPORT
void pk_vector_floor(pk_vector_t *self, float floor);

// Apply log to the vector
POCKETKALDI_EXPORT
void pk_vector_log(pk_vector_t *self);

// Apply scale to the vector
POCKETKALDI_EXPORT
void pk_vector_scale(pk_vector_t *self, float scale);

// Calculates the dot product of self and vec
POCKETKALDI_EXPORT
float pk_vector_dot(const pk_vector_t *self, const pk_vector_t *vec);

// Calculate (self^T dot W)^T, snd stores the output into out
POCKETKALDI_EXPORT
float pk_vector_dotmat(
    const pk_vector_t *self, 
    const pk_matrix_t *W,
    pk_vector_t *out);

// Add vec to self
POCKETKALDI_EXPORT
void pk_vector_add(const pk_vector_t *self, const pk_vector_t *vec);

// Add scale * vec to self
POCKETKALDI_EXPORT
void pk_vector_addscale(
    const pk_vector_t *self,
    float scale,
    const pk_vector_t *vec);

// Add scalar * vec to self
POCKETKALDI_EXPORT
void pk_vector_scalaradd(
    const pk_vector_t *self,
    float scalar,
    const pk_vector_t *vec);

// Subtract vec from self
POCKETKALDI_EXPORT
void pk_vector_subtract(const pk_vector_t *self, const pk_vector_t *vec);

// Add element-wise squaring of vec to self
POCKETKALDI_EXPORT
void pk_vector_add2(const pk_vector_t *self, const pk_vector_t *vec);

// Subtract element-wise squaring of vec from self
POCKETKALDI_EXPORT
void pk_vector_subtract2(const pk_vector_t *self, const pk_vector_t *vec);

// Creates a subvector according to self. It just borrows the memory of self.
// So to avoid issues it couldn't be destroyed and modified.
POCKETKALDI_EXPORT const pk_vector_t 
pk_vector_subvector(const pk_vector_t *self, int start, int dim);

// Destroy a vector
POCKETKALDI_EXPORT
void pk_vector_destroy(pk_vector_t *self);

#endif  // POCKETKALDI_MATRIX_H_
