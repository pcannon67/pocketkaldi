// 2017-01-27

#include "matrix.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
extern "C" {
#include <cblas.h>
}
#include "util.h"
#include "gemm.h"

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
  int num_col;
  int num_row;
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
  num_col = pk_bytebuffer_readint32(&content);
  num_row = pk_bytebuffer_readint32(&content);

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

// Scale the elements of matrix
POCKETKALDI_EXPORT
void pk_matrix_scale(pk_matrix_t *self, float scale) {
  for (int d = 0; d < self->ncol * self->nrow; ++d) {
    self->data[d] *= scale;
  }
}

void pk_matrix_matmat(
    const pk_matrix_t *A,
    const pk_matrix_t *B,
    pk_matrix_t *C) {
  assert(A->nrow == B->nrow && A->ncol == C->nrow && B->ncol == C->ncol &&
         "pk_matrix_matmat: matrix shape mismatch");
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

namespace pocketkaldi {

template<typename Real>
void MatrixBase<Real>::CopyFromMat(const MatrixBase<Real> &M, int trans) {
  if (trans == kNoTrans) {
    assert(num_rows_ == M.NumRows() && num_cols_ == M.NumCols());
    for (int i = 0; i < num_rows_; i++) {
      (*this).Row(i).CopyFromVec(M.Row(i));
    }
  } else {
    assert(num_cols_ == M.NumRows() && num_rows_ == M.NumCols());
    int this_stride = stride_, other_stride = M.Stride();
    Real *this_data = data_;
    const Real *other_data = M.Data();
    for (int i = 0; i < num_rows_; i++) {
      for (int j = 0; j < num_cols_; j++) {
        this_data[i * this_stride + j] = other_data[j * other_stride + i];
      }
    }
  }
}

template<typename Real>
void Matrix<Real>::Swap(Matrix<Real> *other) {
  std::swap(this->data_, other->data_);
  std::swap(this->num_cols_, other->num_cols_);
  std::swap(this->num_rows_, other->num_rows_);
  std::swap(this->stride_, other->stride_);
}

template<typename Real>
void MatrixBase<Real>::Transpose() {
  assert(num_rows_ == num_cols_);
  int M = num_rows_;
  for (int i = 0;i < M;i++) {
    for (int j = 0;j < i;j++) {
      Real &a = (*this)(i, j), &b = (*this)(j, i);
      std::swap(a, b);
    }
  }
}

template<typename Real>
void Matrix<Real>::Init(
    const int rows,
    const int cols,
    const int stride_type) {
  if (rows * cols == 0) {
    assert(rows == 0 && cols == 0);
    this->num_rows_ = 0;
    this->num_cols_ = 0;
    this->stride_ = 0;
    this->data_ = NULL;
    return;
  }
  assert(rows > 0 && cols > 0);
  int skip, stride;
  size_t size;
  void *data;  // aligned memory block

  // compute the size of skip and real cols
  skip = ((16 / sizeof(Real)) - cols % (16 / sizeof(Real)))
      % (16 / sizeof(Real));
  stride = cols + skip;
  size = static_cast<size_t>(rows) * static_cast<size_t>(stride)
      * sizeof(Real);

  // allocate the memory and set the right dimensions and parameters
  if (posix_memalign(&data, 32, size) == 0) {
    MatrixBase<Real>::data_ = static_cast<Real *> (data);
    MatrixBase<Real>::num_rows_ = rows;
    MatrixBase<Real>::num_cols_ = cols;
    MatrixBase<Real>::stride_  = (stride_type == kDefaultStride ? 
      stride : cols);
  } else {
    throw std::bad_alloc();
  }
}

template<typename Real>
void Matrix<Real>::Resize(int rows,
                          int cols,
                          int resize_type,
                          int stride_type) {
  // the next block uses recursion to handle what we have to do if
  // resize_type == kCopyData.
  if (resize_type == kCopyData) {
    if (this->data_ == NULL || rows == 0) {
      // nothing to copy.
      resize_type = kSetZero;  
    } else if (rows == this->num_rows_ && cols == this->num_cols_) { 
      // nothing to do.
      return; 
    } else {
      // set tmp to a matrix of the desired size; if new matrix
      // is bigger in some dimension, zero it.
      int new_resize_type =
          (rows > this->num_rows_ || cols > this->num_cols_) ? 
              kSetZero : kUndefined;
      Matrix<Real> tmp(rows, cols, new_resize_type);
      int rows_min = std::min(rows, this->num_rows_),
          cols_min = std::min(cols, this->num_cols_);
      tmp.Range(0, rows_min, 0, cols_min).
          CopyFromMat(this->Range(0, rows_min, 0, cols_min));
      tmp.Swap(this);
      // and now let tmp go out of scope, deleting what was in *this.
      return;
    }
  }
  // At this point, resize_type == kSetZero or kUndefined.

  if (MatrixBase<Real>::data_ != NULL) {
    if (rows == MatrixBase<Real>::num_rows_
        && cols == MatrixBase<Real>::num_cols_) {
      if (resize_type == kSetZero)
        this->SetZero();
      return;
    }
    else
      Destroy();
  }
  Init(rows, cols, stride_type);
  if (resize_type == kSetZero) MatrixBase<Real>::SetZero();
}

template<typename Real>
void Matrix<Real>::Destroy() {
  // we need to free the data block if it was defined
  if (NULL != MatrixBase<Real>::data_) {
    free( MatrixBase<Real>::data_);
  }
  MatrixBase<Real>::data_ = NULL;
  MatrixBase<Real>::num_rows_ = MatrixBase<Real>::num_cols_
      = MatrixBase<Real>::stride_ = 0;
}

template<typename Real>
void VectorBase<Real>::SetZero() {
  memset(data_, 0, dim_ * sizeof(Real));
}

template<typename Real>
Status Matrix<Real>::Read(util::ReadableFile *fd) {
  static const char *kSectionName = "MAT0";

  // Read section name
  int32_t section_size;
  PK_CHECK_STATUS(fd->ReadAndVerifyString(kSectionName));
  PK_CHECK_STATUS(fd->ReadValue<int32_t>(&section_size));

  int32_t num_rows, num_cols;
  PK_CHECK_STATUS(fd->ReadValue<int32_t>(&num_rows));
  PK_CHECK_STATUS(fd->ReadValue<int32_t>(&num_cols));

  // Resize the matrix according to the new size
  Resize(num_rows, num_cols, kUndefined);

  Vector<Real> row_read;
  for (int row_idx = 0; row_idx < this->NumRows(); ++row_idx) {
    SubVector<Real> row = this->Row(row_idx);
    row_read.Read(fd);
    if (row_read.Dim() != this->NumCols()) {
      return Status::Corruption(util::Format(
          "Matrix::Read: row_read.Dim() == {} expected, but {} found: {}",
          this->NumCols(),
          row_read.Dim(),
          fd->filename()));
    }

    row.CopyFromVec(row_read);
  }

  return Status::OK();
}

// Constructor... note that this is not const-safe as it would
// be quite complicated to implement a "const SubMatrix" class that
// would not allow its contents to be changed.
template<typename Real>
SubMatrix<Real>::SubMatrix(const MatrixBase<Real> &M,
                           int ro,
                           int r,
                           int co,
                           int c) {
  if (r == 0 || c == 0) {
    // we support the empty sub-matrix as a special case.
    assert(c == 0 && r == 0);
    this->data_ = NULL;
    this->num_cols_ = 0;
    this->num_rows_ = 0;
    this->stride_ = 0;
    return;
  }
  assert(ro < (M.num_rows_) && co < M.num_cols_ &&
         r <= M.num_rows_ - ro && c <= M.num_cols_ - co);
  // point to the begining of window
  MatrixBase<Real>::num_rows_ = r;
  MatrixBase<Real>::num_cols_ = c;
  MatrixBase<Real>::stride_ = M.Stride();
  MatrixBase<Real>::data_ = M.Data_workaround() +
      static_cast<size_t>(co) +
      static_cast<size_t>(ro) * static_cast<size_t>(M.Stride());
}


template<typename Real>
SubMatrix<Real>::SubMatrix(Real *data,
                           int num_rows,
                           int num_cols,
                           int stride):
    MatrixBase<Real>(data, num_cols, num_rows, stride) { // caution: reversed order!
  if (data == NULL) {
    assert(num_rows * num_cols == 0);
    this->num_rows_ = 0;
    this->num_cols_ = 0;
    this->stride_ = 0;
  } else {
    assert(this->stride_ >= this->num_cols_);
  }
}

template<typename Real>
void MatrixBase<Real>::SetZero() {
  if (num_cols_ == stride_)
    memset(data_, 0, sizeof(Real)*num_rows_*num_cols_);
  else
    for (int row = 0; row < num_rows_; row++)
      memset(data_ + row*stride_, 0, sizeof(Real)*num_cols_);
}

template<typename Real>
void MatrixBase<Real>::SetRand() {
  for (int row = 0; row < num_rows_; row++) {
    for (int col = 0; col < num_cols_; col++) {
      (*this)(row, col) = static_cast<Real>(rand()) / RAND_MAX;
    }
  }
}

template class Matrix<float>;
template class Matrix<double>;
template class MatrixBase<float>;
template class MatrixBase<double>;
template class SubMatrix<float>;
template class SubMatrix<double>;

// C <- A * B
template<typename Real>
void SimpleMatMat(
    const MatrixBase<Real> &A,
    const MatrixBase<Real> &B,
    MatrixBase<Real> *C) {
  assert(B.NumCols() == C->NumCols());
  assert(A.NumRows() == C->NumRows());
  assert(A.NumCols() == B.NumRows());

  C->SetZero();
  for (int row = 0; row < A.NumRows(); ++row) {
    for (int col = 0; col < B.NumCols(); ++col) {
      for (int k = 0; k < A.NumCols(); ++k) {
        (*C)(row, col) += A(row, k) * B(k, col);
      }
    }
  }
}

template 
void SimpleMatMat<float>(
    const MatrixBase<float> &A,
    const MatrixBase<float> &B,
    MatrixBase<float> *C);

void MatMat(
    const MatrixBase<float> &A,
    const MatrixBase<float> &B,
    MatrixBase<float> *C,
    GEMM<float> *sgemm) {
  assert(B.NumCols() == C->NumCols());
  assert(A.NumRows() == C->NumRows());
  assert(A.NumCols() == B.NumRows());

  sgemm->Gemm(
      A.NumRows(),
      B.NumCols(),
      A.NumCols(),
      1.0f,
      A.Data(), A.Stride(), 1,
      B.Data(), B.Stride(), 1,
      0.0f,
      C->Data(), C->Stride(), 1);
}

}  // namespace pocketkaldi
