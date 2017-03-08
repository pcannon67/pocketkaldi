// 2017-01-27

#ifndef POCKETKALDI_MATRIX_H_
#define POCKETKALDI_MATRIX_H_

// pk_matrix_t is a column-major matrix type in pocketkaldi
typedef struct pk_matrix_t {
  int ncol;
  int nrow;
  float *data;
} pk_matrix_t;

// Initialize the matrix with specified rows and columns
void pk_matrix_init(pk_matrix_t *self, int nrow, int ncol);

// Fill the matrix with value
void pk_matrix_fill(pk_matrix_t *self, float val);

// Set the n-th column vector of matrix. The n is specified by col and values
// are in data. It reads nrow float values in data and stores into the column
// vector.
void pk_matrx_setcol(pk_matrix_t *self, int col, const float *data);

// Returns a pointer to the n-th column vector, n is specified by col.
const float *pk_matrix_getcol(const pk_matrix_t *self, int col);

#endif  // POCKETKALDI_MATRIX_H_
