// gemm.cc - Created at 2017-05-29, use the code of ulmBLAS
// 
// The ulmBLAS framework uses code, concepts and ideas of the great BLIS
// framework:
//
// ---
//
// The BLIS framework is licensed under the following license, typically
// known as the "new" or "modified" or "3-clause" BSD license.
//
// 
// Copyright (C) 2014, The University of Texas
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//  - Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  - Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  - Neither the name of The University of Texas nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "gemm.h"

#include <stdlib.h>
#include <unistd.h>
#include <new>

namespace pocketkaldi {


template<typename T>
GEMM<T>::GEMM() : A_(nullptr), B_(nullptr), C_(nullptr) {
  int packed_size = (MC * KC + KC * NC + MR * NR) * sizeof(T);
  if (0 != posix_memalign(
      reinterpret_cast<void **>(&packed_buffer_),
      256,
      packed_size)) {
    throw std::bad_alloc();
  }
  A_ = packed_buffer_;
  B_ = A_ + MC * KC;
  C_ = B_ + KC * NC;
}

template<typename T>
GEMM<T>::~GEMM() {
  free(packed_buffer_);
  packed_buffer_ = nullptr;
  A_ = nullptr;
  B_ = nullptr;
  C_ = nullptr;
}

template<typename T>
void GEMM<T>::Gemm(
    int m, int n, int k,
    T alpha,
    const T *A, int incRowA, int incColA,
    const T *B, int incRowB, int incColB,
    T beta,
    T *C, int incRowC, int incColC) {
  int mb = (m + MC - 1) / MC;
  int nb = (n + NC - 1) / NC;
  int kb = (k + KC - 1) / KC;

  int _mc = m % GEMM<T>::MC;
  int _nc = n % NC;
  int _kc = k % KC;

  int mc, nc, kc;
  int i, j, l;

  T _beta;

  if (alpha == 0.0 || k == 0) {
    Gescal(m, n, beta, C, incRowC, incColC);
    return;
  }

  for (j = 0; j < nb; ++j) {
    nc = (j != nb - 1 || _nc == 0) ? NC : _nc;

    for (l = 0; l < kb; ++l) {
      kc = (l != kb - 1 || _kc == 0) ? KC: _kc;
      _beta = (l == 0) ? beta : 1.0;

      PackB(
          kc, nc,
          &B[l * KC * incRowB + j * NC * incColB],
          incRowB, incColB,
          B_);

      for (i = 0; i < mb; ++i) {
        mc = (i != mb-1 || _mc == 0) ? MC : _mc;

        PackA(
            mc, kc,
            &A[i * MC * incRowA + l * KC * incColA],
            incRowA, incColA,
            A_);

        MacroKernel(
            mc, nc, kc,
            alpha, _beta,
            &C[i * MC * incRowC + j * NC * incColC],
            incRowC, incColC);
      }
    }
  }
}

template<typename T>
void GEMM<T>::Gescal(
    int m, int n,
    T alpha,
    T *X, int incRowX, int incColX) {
  int i, j;

  if (alpha != 0.0) {
    for (j = 0; j < n; ++j) {
      for (i = 0; i < m; ++i) {
        X[i * incRowX + j * incColX] *= alpha;
      }
    }
  } else {
    for (j = 0; j < n; ++j) {
      for (i = 0; i < m; ++i) {
        X[i * incRowX + j * incColX] = 0.0;
      }
    }
  }
}

template<typename T>
void GEMM<T>::Geaxpy(
    int m, int n,
    T alpha,
    const T *X, int incRowX, int incColX,
    T *Y, int incRowY, int incColY) {
  int i, j;

  if (alpha != 1.0) {
    for (j = 0; j < n; ++j) {
      for (i = 0; i < m; ++i) {
        Y[i * incRowY + j * incColY] += alpha * X[i * incRowX + j * incColX];
      }
    }
  } else {
    for (j = 0; j < n; ++j) {
      for (i = 0; i < m; ++i) {
        Y[i * incRowY + j * incColY] += X[i * incRowX + j * incColX];
      }
    }
  }
}

template<typename T>
void GEMM<T>::PackMRxk(int k, const T *A, int incRowA, int incColA, T *buffer) {
  int i, j;

  for (j = 0; j < k; ++j) {
    for (i = 0; i < MR; ++i) {
      buffer[i] = A[i * incRowA];
    }
    buffer += MR;
    A += incColA;
  }
}

template<typename T>
void GEMM<T>::PackA(
    int mc, int kc,
    const T *A, int incRowA, int incColA,
    T *buffer) {
  int mp = mc / MR;
  int _mr = mc % MR;
  int i, j;

  for (i = 0; i < mp; ++i) {
    PackMRxk(kc, A, incRowA, incColA, buffer);
    buffer += kc * MR;
    A += MR * incRowA;
  }
  if (_mr > 0) {
    for (j = 0; j < kc; ++j) {
      for (i = 0; i < _mr; ++i) {
        buffer[i] = A[i * incRowA];
      }
      for (i = _mr; i < MR; ++i) {
        buffer[i] = 0.0;
      }
      buffer += MR;
      A += incColA;
    }
  }
}

template<typename T>
void GEMM<T>::PackkxNR(int k, const T *B, int incRowB, int incColB, T *buffer) {
  int i, j;

  for (i = 0; i < k; ++i) {
    for (j = 0; j < NR; ++j) {
        buffer[j] = B[j * incColB];
    }
    buffer += NR;
    B += incRowB;
  }
}

template<typename T>
void GEMM<T>::PackB(
    int kc, int nc,
    const T *B, int incRowB, int incColB,
    T *buffer) {
  int np  = nc / NR;
  int _nr = nc % NR;
  int i, j;

  for (j = 0; j < np; ++j) {
    PackkxNR(kc, B, incRowB, incColB, buffer);
    buffer += kc * NR;
    B += NR * incColB;
  }
  if (_nr > 0) {
    for (i = 0; i < kc; ++i) {
      for (j = 0; j < _nr; ++j) {
        buffer[j] = B[j * incColB];
      }
      for (j = _nr; j < NR; ++j) {
        buffer[j] = 0.0;
      }
      buffer += NR;
      B += incRowB;
    }
  }
}

template<typename T>
void GEMM<T>::MacroKernel(
    int mc, int nc, int kc,
    T alpha, T beta,
    T *C, int incRowC, int incColC) {
  int mp = (mc + MR - 1) / MR;
  int np = (nc + NR - 1) / NR;
  int _mr = mc % MR;
  int _nr = nc % NR;
  int mr, nr;
  int i, j;

  for (j = 0; j < np; ++j) {
    nr = (j != np - 1 || _nr == 0) ? NR : _nr;

    for (i = 0; i < mp; ++i) {
      mr = (i != mp - 1 || _mr == 0) ? MR : _mr;

      if (mr == MR && nr == NR) {
        Kernel(
          kc,
          &alpha,
          &A_[i * kc * MR],
          &B_[j * kc * NR],
          &beta,
          &C[i * MR * incRowC + j * NR * incColC],
          incRowC, incColC);
      } else {
        T zero_beta = 0.0;
        Kernel(
            kc,
            &alpha,
            &A_[i * kc * MR],
            &B_[j * kc * NR],
            &zero_beta,
            C_, 1, MR);
        Gescal(
            mr, nr,
            beta,
            &C[i * MR * incRowC + j * NR * incColC],
            incRowC, incColC);
        Geaxpy(
            mr, nr,
            1.0,
            C_, 1, MR,
            &C[i * MR * incRowC + j * NR * incColC],
            incRowC, incColC);
      }
    }
  }
}

template class GEMM<float>;

}  // namespace pocketkaldi
