// gemm_haswell.cc -- Created at 2017-05-29
// Use the code of blis
/*

   BLIS    
   An object-based framework for developing high-performance BLAS-like
   libraries.

   Copyright (C) 2014, The University of Texas at Austin

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    - Neither the name of The University of Texas at Austin nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "gemm.h"

namespace pocketkaldi {


#define SGEMM_INPUT_GS_BETA_NZ \
  "vmovlps    (%%rcx        ),  %%xmm0,  %%xmm0  \n\t" \
  "vmovhps    (%%rcx,%%rsi,1),  %%xmm0,  %%xmm0  \n\t" \
  "vmovlps    (%%rcx,%%rsi,2),  %%xmm1,  %%xmm1  \n\t" \
  "vmovhps    (%%rcx,%%r13  ),  %%xmm1,  %%xmm1  \n\t" \
  "vshufps    $0x88,   %%xmm1,  %%xmm0,  %%xmm0  \n\t" \
  "vmovlps    (%%rcx,%%rsi,4),  %%xmm2,  %%xmm2  \n\t" \
  "vmovhps    (%%rcx,%%r15  ),  %%xmm2,  %%xmm2  \n\t" \
  "vmovlps    (%%rcx,%%r13,2),  %%xmm1,  %%xmm1  \n\t" \
  "vmovhps    (%%rcx,%%r10  ),  %%xmm1,  %%xmm1  \n\t" \
  "vshufps    $0x88,   %%xmm1,  %%xmm2,  %%xmm2  \n\t" \
  "vperm2f128 $0x20,   %%ymm2,  %%ymm0,  %%ymm0  \n\t"

#define SGEMM_OUTPUT_GS_BETA_NZ \
  "vextractf128  $1, %%ymm0,  %%xmm2           \n\t" \
  "vmovss            %%xmm0, (%%rcx        )   \n\t" \
  "vpermilps  $0x39, %%xmm0,  %%xmm1           \n\t" \
  "vmovss            %%xmm1, (%%rcx,%%rsi,1)   \n\t" \
  "vpermilps  $0x39, %%xmm1,  %%xmm0           \n\t" \
  "vmovss            %%xmm0, (%%rcx,%%rsi,2)   \n\t" \
  "vpermilps  $0x39, %%xmm0,  %%xmm1           \n\t" \
  "vmovss            %%xmm1, (%%rcx,%%r13  )   \n\t" \
  "vmovss            %%xmm2, (%%rcx,%%rsi,4)   \n\t" \
  "vpermilps  $0x39, %%xmm2,  %%xmm1           \n\t" \
  "vmovss            %%xmm1, (%%rcx,%%r15  )   \n\t" \
  "vpermilps  $0x39, %%xmm1,  %%xmm2           \n\t" \
  "vmovss            %%xmm2, (%%rcx,%%r13,2)   \n\t" \
  "vpermilps  $0x39, %%xmm2,  %%xmm1           \n\t" \
  "vmovss            %%xmm1, (%%rcx,%%r10  )   \n\t"

void bli_sgemm_asm_6x16(
    int64_t k,
    float *alpha,
    float *a,
    float *b,
    float *beta,
    float *c, int64_t rs_c, int64_t cs_c) {
  //void*   a_next = bli_auxinfo_next_a( data );
  //void*   b_next = bli_auxinfo_next_b( data );

  uint64_t   k_iter = k / 4;
  uint64_t   k_left = k % 4;

  __asm__ __volatile__
  (
  "                                            \n\t"
  "vzeroall                                    \n\t" // zero all xmm/ymm registers.
  "                                            \n\t"
  "                                            \n\t"
  "movq                %2, %%rax               \n\t" // load address of a.
  "movq                %3, %%rbx               \n\t" // load address of b.
  //"movq                %9, %%r15               \n\t" // load address of b_next.
  "                                            \n\t"
  "addq           $32 * 4, %%rbx               \n\t"
  "                                            \n\t" // initialize loop by pre-loading
  "vmovaps           -4 * 32(%%rbx), %%ymm0    \n\t"
  "vmovaps           -3 * 32(%%rbx), %%ymm1    \n\t"
  "                                            \n\t"
  "movq                %6, %%rcx               \n\t" // load address of c
  "movq                %7, %%rdi               \n\t" // load rs_c
  "leaq        (,%%rdi,4), %%rdi               \n\t" // rs_c *= sizeof(float)
  "                                            \n\t"
  "leaq   (%%rdi,%%rdi,2), %%r13               \n\t" // r13 = 3*rs_c;
  "leaq   (%%rcx,%%r13,1), %%rdx               \n\t" // rdx = c + 3*rs_c;
  "prefetcht0   7 * 8(%%rcx)                   \n\t" // prefetch c + 0*rs_c
  "prefetcht0   7 * 8(%%rcx,%%rdi)             \n\t" // prefetch c + 1*rs_c
  "prefetcht0   7 * 8(%%rcx,%%rdi,2)           \n\t" // prefetch c + 2*rs_c
  "prefetcht0   7 * 8(%%rdx)                   \n\t" // prefetch c + 3*rs_c
  "prefetcht0   7 * 8(%%rdx,%%rdi)             \n\t" // prefetch c + 4*rs_c
  "prefetcht0   7 * 8(%%rdx,%%rdi,2)           \n\t" // prefetch c + 5*rs_c
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "movq      %0, %%rsi                         \n\t" // i = k_iter;
  "testq  %%rsi, %%rsi                         \n\t" // check i via logical AND.
  "je     1f                                   \n\t" // if i == 0, jump to code that
  "                                            \n\t" // contains the k_left loop.
  "                                            \n\t"
  "                                            \n\t"
  "2:                                          \n\t" // MAIN LOOP
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t" // iteration 0
  "prefetcht0   64 * 4(%%rax)                  \n\t"
  "                                            \n\t"
  "vbroadcastss       0 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss       1 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm4    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm5    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm6    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm7    \n\t"
  "                                            \n\t"
  "vbroadcastss       2 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss       3 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm8    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm9    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm10   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm11   \n\t"
  "                                            \n\t"
  "vbroadcastss       4 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss       5 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm12   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm13   \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm14   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm15   \n\t"
  "                                            \n\t"
  "vmovaps           -2 * 32(%%rbx), %%ymm0    \n\t"
  "vmovaps           -1 * 32(%%rbx), %%ymm1    \n\t"
  "                                            \n\t"
  "                                            \n\t" // iteration 1
  "vbroadcastss       6 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss       7 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm4    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm5    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm6    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm7    \n\t"
  "                                            \n\t"
  "vbroadcastss       8 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss       9 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm8    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm9    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm10   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm11   \n\t"
  "                                            \n\t"
  "vbroadcastss      10 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss      11 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm12   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm13   \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm14   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm15   \n\t"
  "                                            \n\t"
  "vmovaps            0 * 32(%%rbx), %%ymm0    \n\t"
  "vmovaps            1 * 32(%%rbx), %%ymm1    \n\t"
  "                                            \n\t"
  "                                            \n\t" // iteration 2
  "prefetcht0   76 * 4(%%rax)                  \n\t"
  "                                            \n\t"
  "vbroadcastss      12 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss      13 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm4    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm5    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm6    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm7    \n\t"
  "                                            \n\t"
  "vbroadcastss      14 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss      15 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm8    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm9    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm10   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm11   \n\t"
  "                                            \n\t"
  "vbroadcastss      16 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss      17 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm12   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm13   \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm14   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm15   \n\t"
  "                                            \n\t"
  "vmovaps            2 * 32(%%rbx), %%ymm0    \n\t"
  "vmovaps            3 * 32(%%rbx), %%ymm1    \n\t"
  "                                            \n\t"
  "                                            \n\t" // iteration 3
  "vbroadcastss      18 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss      19 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm4    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm5    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm6    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm7    \n\t"
  "                                            \n\t"
  "vbroadcastss      20 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss      21 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm8    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm9    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm10   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm11   \n\t"
  "                                            \n\t"
  "vbroadcastss      22 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss      23 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm12   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm13   \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm14   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm15   \n\t"
  "                                            \n\t"
  "addq          $4 *  6 * 4, %%rax            \n\t" // a += 4*6  (unroll x mr)
  "addq          $4 * 16 * 4, %%rbx            \n\t" // b += 4*16 (unroll x nr)
  "                                            \n\t"
  "vmovaps           -4 * 32(%%rbx), %%ymm0    \n\t"
  "vmovaps           -3 * 32(%%rbx), %%ymm1    \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "decq   %%rsi                                \n\t" // i -= 1;
  "jne    2b                                   \n\t" // iterate again if i != 0.
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "1:                                          \n\t"
  "                                            \n\t"
  "movq      %1, %%rsi                         \n\t" // i = k_left;
  "testq  %%rsi, %%rsi                         \n\t" // check i via logical AND.
  "je     4f                                   \n\t" // if i == 0, we're done; jump to end.
  "                                            \n\t" // else, we prepare to enter k_left loop.
  "                                            \n\t"
  "                                            \n\t"
  "3:                                          \n\t" // EDGE LOOP
  "                                            \n\t"
  "prefetcht0   64 * 4(%%rax)                  \n\t"
  "                                            \n\t"
  "vbroadcastss       0 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss       1 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm4    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm5    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm6    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm7    \n\t"
  "                                            \n\t"
  "vbroadcastss       2 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss       3 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm8    \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm9    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm10   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm11   \n\t"
  "                                            \n\t"
  "vbroadcastss       4 *  4(%%rax), %%ymm2    \n\t"
  "vbroadcastss       5 *  4(%%rax), %%ymm3    \n\t"
  "vfmadd231ps       %%ymm0, %%ymm2, %%ymm12   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm2, %%ymm13   \n\t"
  "vfmadd231ps       %%ymm0, %%ymm3, %%ymm14   \n\t"
  "vfmadd231ps       %%ymm1, %%ymm3, %%ymm15   \n\t"
  "                                            \n\t"
  "addq          $1 *  6 * 4, %%rax            \n\t" // a += 1*6  (unroll x mr)
  "addq          $1 * 16 * 4, %%rbx            \n\t" // b += 1*16 (unroll x nr)
  "                                            \n\t"
  "vmovaps           -4 * 32(%%rbx), %%ymm0    \n\t"
  "vmovaps           -3 * 32(%%rbx), %%ymm1    \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "decq   %%rsi                                \n\t" // i -= 1;
  "jne    3b                                   \n\t" // iterate again if i != 0.
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "4:                                          \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "movq         %4, %%rax                      \n\t" // load address of alpha
  "movq         %5, %%rbx                      \n\t" // load address of beta 
  "vbroadcastss    (%%rax), %%ymm0             \n\t" // load alpha and duplicate
  "vbroadcastss    (%%rbx), %%ymm3             \n\t" // load beta and duplicate
  "                                            \n\t"
  "vmulps           %%ymm0,  %%ymm4,  %%ymm4   \n\t" // scale by alpha
  "vmulps           %%ymm0,  %%ymm5,  %%ymm5   \n\t"
  "vmulps           %%ymm0,  %%ymm6,  %%ymm6   \n\t"
  "vmulps           %%ymm0,  %%ymm7,  %%ymm7   \n\t"
  "vmulps           %%ymm0,  %%ymm8,  %%ymm8   \n\t"
  "vmulps           %%ymm0,  %%ymm9,  %%ymm9   \n\t"
  "vmulps           %%ymm0,  %%ymm10, %%ymm10  \n\t"
  "vmulps           %%ymm0,  %%ymm11, %%ymm11  \n\t"
  "vmulps           %%ymm0,  %%ymm12, %%ymm12  \n\t"
  "vmulps           %%ymm0,  %%ymm13, %%ymm13  \n\t"
  "vmulps           %%ymm0,  %%ymm14, %%ymm14  \n\t"
  "vmulps           %%ymm0,  %%ymm15, %%ymm15  \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "movq                %8, %%rsi               \n\t" // load cs_c
  "leaq        (,%%rsi,4), %%rsi               \n\t" // rsi = cs_c * sizeof(float)
  "                                            \n\t"
  "leaq   (%%rcx,%%rsi,8), %%rdx               \n\t" // load address of c +  8*cs_c;
  "                                            \n\t"
  "leaq   (%%rsi,%%rsi,2), %%r13               \n\t" // r13 = 3*cs_c;
  "leaq   (%%rsi,%%rsi,4), %%r15               \n\t" // r15 = 5*cs_c;
  "leaq   (%%r13,%%rsi,4), %%r10               \n\t" // r10 = 7*cs_c;
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t" // now avoid loading C if beta == 0
  "                                            \n\t"
  "vxorps    %%ymm0,  %%ymm0,  %%ymm0          \n\t" // set ymm0 to zero.
  "vucomiss  %%xmm0,  %%xmm3                   \n\t" // set ZF if beta == 0.
  "je      7f                                  \n\t" // if ZF = 1, jump to beta == 0 case
  "                                            \n\t"
  "                                            \n\t"
  "cmpq       $4, %%rsi                        \n\t" // set ZF if (4*cs_c) == 4.
  "jz      6f                                  \n\t" // jump to row storage case
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "5:                                          \n\t"
  "                                            \n\t"
  "                                            \n\t"
  SGEMM_INPUT_GS_BETA_NZ
  "vfmadd213ps      %%ymm4,  %%ymm3,  %%ymm0   \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  SGEMM_INPUT_GS_BETA_NZ
  "vfmadd213ps      %%ymm6,  %%ymm3,  %%ymm0   \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  SGEMM_INPUT_GS_BETA_NZ
  "vfmadd213ps      %%ymm8,  %%ymm3,  %%ymm0   \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  SGEMM_INPUT_GS_BETA_NZ
  "vfmadd213ps      %%ymm10, %%ymm3,  %%ymm0   \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  SGEMM_INPUT_GS_BETA_NZ
  "vfmadd213ps      %%ymm12, %%ymm3,  %%ymm0   \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  SGEMM_INPUT_GS_BETA_NZ
  "vfmadd213ps      %%ymm14, %%ymm3,  %%ymm0   \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  //"addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "movq      %%rdx, %%rcx                      \n\t" // rcx = c + 8*cs_c
  "                                            \n\t"
  "                                            \n\t"
  SGEMM_INPUT_GS_BETA_NZ
  "vfmadd213ps      %%ymm5,  %%ymm3,  %%ymm0   \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  SGEMM_INPUT_GS_BETA_NZ
  "vfmadd213ps      %%ymm7,  %%ymm3,  %%ymm0   \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  SGEMM_INPUT_GS_BETA_NZ
  "vfmadd213ps      %%ymm9,  %%ymm3,  %%ymm0   \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  SGEMM_INPUT_GS_BETA_NZ
  "vfmadd213ps      %%ymm11, %%ymm3,  %%ymm0   \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  SGEMM_INPUT_GS_BETA_NZ
  "vfmadd213ps      %%ymm13, %%ymm3,  %%ymm0   \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  SGEMM_INPUT_GS_BETA_NZ
  "vfmadd213ps      %%ymm15, %%ymm3,  %%ymm0   \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  //"addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "jmp    10f                                  \n\t" // jump to end.
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "6:                                          \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "vfmadd231ps      (%%rcx), %%ymm3,  %%ymm4   \n\t"
  "vmovups          %%ymm4,  (%%rcx)           \n\t"
  "addq      %%rdi, %%rcx                      \n\t"
  "vfmadd231ps      (%%rdx), %%ymm3,  %%ymm5   \n\t"
  "vmovups          %%ymm5,  (%%rdx)           \n\t"
  "addq      %%rdi, %%rdx                      \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "vfmadd231ps      (%%rcx), %%ymm3,  %%ymm6   \n\t"
  "vmovups          %%ymm6,  (%%rcx)           \n\t"
  "addq      %%rdi, %%rcx                      \n\t"
  "vfmadd231ps      (%%rdx), %%ymm3,  %%ymm7   \n\t"
  "vmovups          %%ymm7,  (%%rdx)           \n\t"
  "addq      %%rdi, %%rdx                      \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "vfmadd231ps      (%%rcx), %%ymm3,  %%ymm8   \n\t"
  "vmovups          %%ymm8,  (%%rcx)           \n\t"
  "addq      %%rdi, %%rcx                      \n\t"
  "vfmadd231ps      (%%rdx), %%ymm3,  %%ymm9   \n\t"
  "vmovups          %%ymm9,  (%%rdx)           \n\t"
  "addq      %%rdi, %%rdx                      \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "vfmadd231ps      (%%rcx), %%ymm3,  %%ymm10  \n\t"
  "vmovups          %%ymm10, (%%rcx)           \n\t"
  "addq      %%rdi, %%rcx                      \n\t"
  "vfmadd231ps      (%%rdx), %%ymm3,  %%ymm11  \n\t"
  "vmovups          %%ymm11, (%%rdx)           \n\t"
  "addq      %%rdi, %%rdx                      \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "vfmadd231ps      (%%rcx), %%ymm3,  %%ymm12  \n\t"
  "vmovups          %%ymm12, (%%rcx)           \n\t"
  "addq      %%rdi, %%rcx                      \n\t"
  "vfmadd231ps      (%%rdx), %%ymm3,  %%ymm13  \n\t"
  "vmovups          %%ymm13, (%%rdx)           \n\t"
  "addq      %%rdi, %%rdx                      \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "vfmadd231ps      (%%rcx), %%ymm3,  %%ymm14  \n\t"
  "vmovups          %%ymm14, (%%rcx)           \n\t"
  //"addq      %%rdi, %%rcx                      \n\t"
  "vfmadd231ps      (%%rdx), %%ymm3,  %%ymm15  \n\t"
  "vmovups          %%ymm15, (%%rdx)           \n\t"
  //"addq      %%rdi, %%rdx                      \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "jmp    10f                                  \n\t" // jump to end.
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "7:                                          \n\t"
  "                                            \n\t"
  "cmpq       $4, %%rsi                        \n\t" // set ZF if (4*cs_c) == 4.
  "jz      8f                                  \n\t" // jump to row storage case
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "8:                                          \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "vmovaps           %%ymm4,  %%ymm0           \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "vmovaps           %%ymm6,  %%ymm0           \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "vmovaps           %%ymm8,  %%ymm0           \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "vmovaps           %%ymm10, %%ymm0           \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "vmovaps           %%ymm12, %%ymm0           \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "vmovaps           %%ymm14, %%ymm0           \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  //"addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "movq      %%rdx, %%rcx                      \n\t" // rcx = c + 8*cs_c
  "                                            \n\t"
  "                                            \n\t"
  "vmovaps           %%ymm5,  %%ymm0           \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "vmovaps           %%ymm7,  %%ymm0           \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "vmovaps           %%ymm9,  %%ymm0           \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "vmovaps           %%ymm11, %%ymm0           \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "vmovaps           %%ymm13, %%ymm0           \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  "addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "vmovaps           %%ymm15, %%ymm0           \n\t"
  SGEMM_OUTPUT_GS_BETA_NZ
  //"addq      %%rdi, %%rcx                      \n\t" // c += rs_c;
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "jmp    10f                                  \n\t" // jump to end.
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "9:                                          \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "vmovups          %%ymm4,  (%%rcx)           \n\t"
  "addq      %%rdi, %%rcx                      \n\t"
  "vmovups          %%ymm5,  (%%rdx)           \n\t"
  "addq      %%rdi, %%rdx                      \n\t"
  "                                            \n\t"
  "vmovups          %%ymm6,  (%%rcx)           \n\t"
  "addq      %%rdi, %%rcx                      \n\t"
  "vmovups          %%ymm7,  (%%rdx)           \n\t"
  "addq      %%rdi, %%rdx                      \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "vmovups          %%ymm8,  (%%rcx)           \n\t"
  "addq      %%rdi, %%rcx                      \n\t"
  "vmovups          %%ymm9,  (%%rdx)           \n\t"
  "addq      %%rdi, %%rdx                      \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "vmovups          %%ymm10, (%%rcx)           \n\t"
  "addq      %%rdi, %%rcx                      \n\t"
  "vmovups          %%ymm11, (%%rdx)           \n\t"
  "addq      %%rdi, %%rdx                      \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "vmovups          %%ymm12, (%%rcx)           \n\t"
  "addq      %%rdi, %%rcx                      \n\t"
  "vmovups          %%ymm13, (%%rdx)           \n\t"
  "addq      %%rdi, %%rdx                      \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "vmovups          %%ymm14, (%%rcx)           \n\t"
  //"addq      %%rdi, %%rcx                      \n\t"
  "vmovups          %%ymm15, (%%rdx)           \n\t"
  //"addq      %%rdi, %%rdx                      \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "                                            \n\t"
  "10:                                         \n\t"
  "                                            \n\t"

  : // output operands (none)
  : // input operands
    "m" (k_iter), // 0
    "m" (k_left), // 1
    "m" (a),      // 2
    "m" (b),      // 3
    "m" (alpha),  // 4
    "m" (beta),   // 5
    "m" (c),      // 6
    "m" (rs_c),   // 7
    "m" (cs_c)/*,   // 8
    "m" (b_next), // 9
    "m" (a_next)*/  // 10
  : // register clobber list
    "rax", "rbx", "rcx", "rdx", "rsi", "rdi", 
    "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15",
    "xmm0", "xmm1", "xmm2", "xmm3",
    "xmm4", "xmm5", "xmm6", "xmm7",
    "xmm8", "xmm9", "xmm10", "xmm11",
    "xmm12", "xmm13", "xmm14", "xmm15",
    "memory"
  );
}

template<>
void GEMM<float>::Kernel(
    int64_t k,
    float *alpha,
    float *a,
    float *b,
    float *beta,
    float *c, int64_t rs_c, int64_t cs_c) {
  bli_sgemm_asm_6x16(k, alpha, a, b, beta, c, rs_c, cs_c);
}

}  // namespace pocketkaldi
