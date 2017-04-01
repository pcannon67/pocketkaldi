// 2017-01-27



#ifndef POCKETKALDI_MATRIX_H_
#define POCKETKALDI_MATRIX_H_

#define PK_VECTOR_SECTION "VEC0"
#define PK_MATRIX_SECTION "MAT0"

#include <stdio.h>
#include "pocketkaldi.h"
#include "util.h"
#include "vector.h"

typedef struct pk_vector_t pk_vector_t;

// pk_matrix_t is a column-major matrix type in pocketkaldi
typedef struct pk_matrix_t {
  int ncol;
  int nrow;
  float *data;
} pk_matrix_t;

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
pk_vector_t pk_matrix_getcol(const pk_matrix_t *self, int col);

#endif  // POCKETKALDI_MATRIX_H_
