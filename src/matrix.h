// 2017-01-27



#ifndef POCKETKALDI_MATRIX_H_
#define POCKETKALDI_MATRIX_H_

#define PK_MATRIX_SECTION "MAT0"

#include <stdio.h>
#include <math.h>
#include "pocketkaldi.h"
#include "util.h"
#include "vector.h"
#include "gemm.h"

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

// Multiplies two matrices C <- A^T dot B
POCKETKALDI_EXPORT
void pk_matrix_matmat(
    const pk_matrix_t *A,
    const pk_matrix_t *B,
    pk_matrix_t *C);

// Resize the matrix
POCKETKALDI_EXPORT
void pk_matrix_resize(pk_matrix_t *self, int nrow, int ncol);

// Copy the matrix from src to dest
POCKETKALDI_EXPORT
void pk_matrix_copy(pk_matrix_t *dest, const pk_matrix_t *src);

// Scale the elements of matrix
POCKETKALDI_EXPORT
void pk_matrix_scale(pk_matrix_t *self, float scale);

// Destroy the matrix
POCKETKALDI_EXPORT
void pk_matrix_destroy(pk_matrix_t *self);

// Returns a borrowed vector of the specified column. The values (data) could be
// changed. And since it is just borrowed from the matrix, destory is not
// needed. 
POCKETKALDI_EXPORT
pk_vector_t pk_matrix_getcol(const pk_matrix_t *self, int col);

namespace pocketkaldi {

template<typename Real>
class SubMatrix;

/// Base class which provides matrix operations not involving resizing
/// or allocation.   Classes Matrix and SubMatrix inherit from it and take care
/// of allocation and resizing.
template<typename Real>
class MatrixBase {
 public:
  enum {
    kTrans,
    kNoTrans
  };

  /// Returns number of rows (or zero for emtpy matrix).
  inline int NumRows() const { return num_rows_; }

  /// Returns number of columns (or zero for emtpy matrix).
  inline int NumCols() const { return num_cols_; }

  /// Gives pointer to raw data (const).
  inline const Real* Data() const {
    return data_;
  }

  /// Sets matrix to zero.
  void SetZero();

  /// Gives pointer to raw data (non-const).
  inline Real* Data() { return data_; }

  /// Stride (distance in memory between each row).  Will be >= NumCols.
  inline int Stride() const {  return stride_; }

  /// Sets to random values between 0 and 1
  void SetRand();

  /// Indexing operator, non-const
  /// (only checks sizes if compiled with -DKALDI_PARANOID)
  inline Real&  operator() (int r, int c) {
    assert(r < num_rows_ && c < num_cols_);
    return *(data_ + r * stride_ + c);
  }

  /// Indexing operator, const
  /// (only checks sizes if compiled with -DKALDI_PARANOID)
  inline const Real operator() (int r, int c) const {
    assert(r < num_rows_ && c < num_cols_);
    return *(data_ + r * stride_ + c);
  }

  /// Return specific row of matrix [const].
  inline const SubVector<Real> Row(int i) const {
    assert(i < num_rows_);
    return SubVector<Real>(data_ + (i * stride_), NumCols());
  }

  /// Return specific row of matrix.
  inline SubVector<Real> Row(int i) {
    assert(i < num_rows_);
    return SubVector<Real>(data_ + (i * stride_), NumCols());
  }

  /// Copy given matrix. (no resize is done).
  void CopyFromMat(const MatrixBase<Real> &M, int trans = kNoTrans);

  /// Transpose the matrix.  Works for non-square
  /// matrices as well as square ones.
  void Transpose();

  /// Return a sub-part of matrix.
  inline SubMatrix<Real> Range(int row_offset,
                               int num_rows,
                               int col_offset,
                               int num_cols) const {
    return SubMatrix<Real>(*this, row_offset, num_rows,
                           col_offset, num_cols);
  }


  friend class SubMatrix<Real>;
 protected:
  ///  Initializer, callable only from child.
  explicit MatrixBase(Real *data, int cols, int rows, int stride) :
    data_(data), num_cols_(cols), num_rows_(rows), stride_(stride) {
  }

  ///  Initializer, callable only from child.
  /// Empty initializer, for un-initialized matrix.
  explicit MatrixBase(): data_(NULL) {
  }

  // Make sure pointers to MatrixBase cannot be deleted.
  ~MatrixBase() { }

  /// A workaround that allows SubMatrix to get a pointer to non-const data
  /// for const Matrix. Unfortunately C++ does not allow us to declare a
  /// "public const" inheritance or anything like that, so it would require
  /// a lot of work to make the SubMatrix class totally const-correct--
  /// we would have to override many of the Matrix functions.
  inline Real*  Data_workaround() const {
    return data_;
  }

  /// data memory area
  Real* data_;

  /// these atributes store the real matrix size as it is stored in memory
  /// including memalignment
  int num_cols_;   /// < Number of columns
  int num_rows_;   /// < Number of rows
  int stride_;
};

/// A class for storing matrices.
template<typename Real>
class Matrix : public MatrixBase<Real> {
 public:
  enum {
    kSetZero,
    kUndefined,
    kCopyData
  };
  enum {
    kDefaultStride,
    kStrideEqualNumCols,
  };

  /// Empty constructor.
  Matrix() {
    this->data_ = nullptr;
    this->num_cols_ = 0;
    this->num_rows_ = 0;
    this->stride_ = 0;
  }

  /// Basic constructor.
  Matrix(int r,
         int c,
         int resize_type = kSetZero,
         int stride_type = kStrideEqualNumCols) { 
    Resize(r, c, resize_type, stride_type); 
  }

  /// Swaps the contents of *this and *other.  Shallow swap.
  void Swap(Matrix<Real> *other);

  /// Distructor to free matrices.
  ~Matrix() { Destroy(); }

  /// Sets matrix to a specified size (zero is OK as long as both r and c are
  /// zero).  The value of the new data depends on resize_type:
  ///   -if kSetZero, the new data will be zero
  ///   -if kUndefined, the new data will be undefined
  ///   -if kCopyData, the new data will be the same as the old data in any
  ///      shared positions, and zero elsewhere.
  ///
  /// You can set stride_type to kStrideEqualNumCols to force the stride
  /// to equal the number of columns; by default it is set so that the stride
  /// in bytes is a multiple of 16.
  ///
  /// This function takes time proportional to the number of data elements.
  void Resize(int r,
              int c,
              int resize_type = kSetZero,
              int stride_type = kStrideEqualNumCols);

  // Read matrix from ReadableFile
  Status Read(util::ReadableFile *fd);

 private:
  /// Deallocates memory and sets to empty matrix (dimension 0, 0).
  void Destroy();

  /// Init assumes the current class contents are invalid (i.e. junk or have
  /// already been freed), and it sets the matrix to newly allocated memory with
  /// the specified number of rows and columns.  r == c == 0 is acceptable.  The data
  /// memory contents will be undefined.
  void Init(int r,
            int c,
            int stride_type);

};

template<typename Real>
class SubMatrix : public MatrixBase<Real> {
 public:
  // Initialize a SubMatrix from part of a matrix; this is
  // a bit like A(b:c, d:e) in Matlab.
  // This initializer is against the proper semantics of "const", since
  // SubMatrix can change its contents.  It would be hard to implement
  // a "const-safe" version of this class.
  SubMatrix(const MatrixBase<Real> &T,
            int ro,  // row offset, 0 < ro < NumRows()
            int r,   // number of rows, r > 0
            int co,  // column offset, 0 < co < NumCols()
            int c);   // number of columns, c > 0

  // This initializer is mostly intended for use in CuMatrix and related
  // classes.  Be careful!
  SubMatrix(Real *data,
            int num_rows,
            int num_cols,
            int stride);

  ~SubMatrix<Real>() {}
};

// C <- A * B
template<typename Real>
void SimpleMatMat(
    const MatrixBase<Real> &A,
    const MatrixBase<Real> &B,
    MatrixBase<Real> *C);

// C <- A * B
void MatMat(
    const MatrixBase<float> &A,
    const MatrixBase<float> &B,
    MatrixBase<float> *C,
    GEMM<float> *sgemm);

}  // namespace pocketkaldi

#endif  // POCKETKALDI_MATRIX_H_
