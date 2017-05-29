// gemm.h - Created at 2017-05-27, use the code of ulmBLAS
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

#include <stdint.h>

#ifndef POCKETKALDI_GEMM_H_
#define POCKETKALDI_GEMM_H_

namespace pocketkaldi {

// Constants for fp32
template<typename D>
struct GEMMConst {};
template<>
struct GEMMConst<float> {
  static constexpr int MC = 288;
  static constexpr int KC = 512;
  static constexpr int NC = 4096;
  static constexpr int MR = 6;
  static constexpr int NR = 16;
};

template<typename T>
class GEMM {
 public:
  GEMM();
  ~GEMM();

  static constexpr int MC = GEMMConst<T>::MC;
  static constexpr int KC = GEMMConst<T>::KC;
  static constexpr int NC = GEMMConst<T>::NC;
  static constexpr int MR = GEMMConst<T>::MR;
  static constexpr int NR = GEMMConst<T>::NR;

  // Compute C <- beta * C + alpha * A * B
  void Gemm(
      int m, int n, int k,
      T alpha,
      const T *A, int incRowA, int incColA,
      const T *B, int incRowB, int incColB,
      T beta,
      T *C, int incRowC, int incColC);

  // Compute X *= alpha
  void Gescal(
      int m, int n,
      T alpha,
      T *X, int incRowX, int incColX);

  // Compute Y += alpha * X
  void Geaxpy(
      int m, int n,
      T alpha,
      const T *X, int incRowX, int incColX,
      T *Y, int incRowY, int incColY);

 private:
  T *packed_buffer_;
  T *A_;
  T *B_;
  T *C_;

  // Packing complete panels from A (i.e. without padding)
  void PackMRxk(int k, const T *A, int incRowA, int incColA, T *buffer);
  
  // Packing panels from A with padding if required
  void PackA(
      int mc, int kc,
      const T *A, int incRowA, int incColA,
      T *buffer);

  // Packing complete panels from B (i.e. without padding)
  void PackkxNR(int k, const T *B, int incRowB, int incColB, T *buffer);

  // Packing panels from B with padding if required
  void PackB(
      int kc, int nc,
      const T *B, int incRowB, int incColB,
      T *buffer);

  // Macro Kernel for the multiplication of blocks of A and B.  We assume that
  // these blocks were previously packed to buffers A_ and B_.
  void MacroKernel(
      int mc, int nc, int kc,
      T alpha, T beta,
      T *C, int incRowC, int incColC);

  // Computes C <- beta * C + alpha * A * B
  // Where A is MR * k, B is k * NR
  void Kernel(
      int64_t k,
      T *alpha,
      T *a,
      T *b,
      T *beta,
      T *c, int64_t rs_c, int64_t cs_c);
};

}  // namespace pocketkaldi

#endif  // POCKETKALDI_GEMM_H_
