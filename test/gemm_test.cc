// gemm.cc -- Created at 2017-05-29

#include <math.h>
#include <tuple>
#include "matrix.h"
#include "gemm.h"

using pocketkaldi::Matrix;
using pocketkaldi::MatrixBase;
using pocketkaldi::GEMM;
using pocketkaldi::SimpleMatMat;
using pocketkaldi::MatMat;

// Returns the maximun difference of dimensions
float CompareMatrix(
    const MatrixBase<float> &A,
    const MatrixBase<float> &B) {
  assert(A.NumCols() == B.NumCols());
  assert(A.NumRows() == B.NumRows());

  float max_diff = 0.0f;
  for (int row_idx = 0; row_idx < A.NumRows(); ++row_idx) {
    for (int col_idx = 0; col_idx < A.NumCols(); ++col_idx) {
      float diff = fabs(A(row_idx, col_idx) - B(row_idx, col_idx));
      if (diff > max_diff) max_diff = diff;
    }
  }

  return max_diff;
}

void TestSgemm() {
  GEMM<float> sgemm;
  Matrix<float> A;
  Matrix<float> B;
  Matrix<float> C;
  Matrix<float> CRef;

  std::vector<std::tuple<int, int, int>> test_sizes {
    std::make_tuple(512, 512, 512),
    std::make_tuple(100, 100, 1),
    std::make_tuple(1, 1, 1),
    std::make_tuple(121, 233, 17)
  };
  for (const std::tuple<int, int, int> &test_size : test_sizes) {
    int m = std::get<0>(test_size);
    int n = std::get<1>(test_size);
    int k = std::get<2>(test_size);

    A.Resize(m, k);
    A.SetRand();
    B.Resize(k, n);
    B.SetRand();
    C.Resize(m, n);
    CRef.Resize(m, n);

    SimpleMatMat(A, B, &CRef);
    MatMat(A, B, &C, &sgemm);

    assert(CompareMatrix(C, CRef) < 0.01);
  }
}

int main() {
  TestSgemm();
  return 0;
}
