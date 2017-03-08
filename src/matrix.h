// 2017-01-27

#ifndef POCKETKALDI_MATRIX_H_
#define POCKETKALDI_MATRIX_H_

#include "pocketkaldi.h"

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

// Fill the matrix with value
POCKETKALDI_EXPORT
void pk_matrix_fill(pk_matrix_t *self, float val);

// Resize the matrix
POCKETKALDI_EXPORT
void pk_matrix_resize(pk_matrix_t *self, int nrow, int ncol);

// Destroy the matrix
POCKETKALDI_EXPORT
void pk_matrix_destroy(pk_matrix_t *self);

// Returns a borrowed vector of the specified column. The values (data) could be
// changed. And since it is just borrowed from the matrix, destory is not
// needed. 
POCKETKALDI_EXPORT
const pk_vector_t pk_matrix_getcol(const pk_matrix_t *self, int col);

// Initialize the vector
POCKETKALDI_EXPORT
void pk_vector_init(pk_vector_t *self, int dim);

// Resize the vector
POCKETKALDI_EXPORT
void pk_vector_resize(pk_vector_t *self, int dim);

// Fill the vector with value
POCKETKALDI_EXPORT
void pk_vector_fill(pk_vector_t *self, float value);

// Copy n values from source to the vector. n should be less than self->dim
POCKETKALDI_EXPORT
void pk_vector_copyfrom(pk_vector_t *self, const float *source, int n);

// Calculates the dot product of self and vec
POCKETKALDI_EXPORT
float pk_vector_dot(const pk_vector_t *self, const pk_vector_t *vec);

// Creates a subvector according to self. It just borrows the memory of self.
// So to avoid issues it couldn't be destroyed and modified.
POCKETKALDI_EXPORT const pk_vector_t 
pk_vector_subvector(const pk_vector_t *self, int start, int dim);

// Destroy a vector
POCKETKALDI_EXPORT
void pk_vector_destroy(pk_vector_t *self);

#endif  // POCKETKALDI_MATRIX_H_
