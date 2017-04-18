// matrix/srfft.cc

// Copyright 2009-2011  Microsoft Corporation;  Go Vivace Inc.

// See ../../COPYING for clarification regarding multiple authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.
//

// This file includes a modified version of code originally published in Malvar,
// H., "Signal processing with lapped transforms," Artech House, Inc., 1992. The
// current copyright holder of the original code, Henrique S. Malvar, has given
// his permission for the release of this modified version under the Apache
// License v2.0.


#include "srfft.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


#ifndef M_2PI
#define M_2PI 6.283185307179586476925286766559005
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

static void compute_tables(pk_srfft_t *self) {
  int imax, lg2, i, j;
  int m, m2, m4, m8, nel, n;
  float *cn, *spcn, *smcn, *c3n, *spc3n, *smc3n;
  float ang, c, s;

  lg2 = self->logn >> 1;
  if (self->logn & 1) lg2++;
  self->brseed = (int *)malloc(sizeof(int) * (1 << lg2));
  self->brseed[0] = 0;
  self->brseed[1] = 1;
  for (j = 2; j <= lg2; j++) {
    imax = 1 << (j - 1);
    for (i = 0; i < imax; i++) {
      self->brseed[i] <<= 1;
      self->brseed[i + imax] = self->brseed[i] + 1;
    }
  }

  if (self->logn < 4) {
    self->tab = NULL;
  } else {
    self->tab = (float **)malloc(sizeof(float *) * self->logn - 3);
    for (i = self->logn; i >= 4 ;i--) {
      /* Compute a few constants */
      m = 1 << i; m2 = m / 2; m4 = m2 / 2; m8 = m4 /2;

      /* Allocate memory for tables */
      nel = m4 - 2;

      self->tab[i - 4] = (float *)malloc(sizeof(float) * (6 * nel));

      /* Initialize pointers */
      cn = self->tab[i - 4]; spcn  = cn + nel;  smcn  = spcn + nel;
      c3n = smcn + nel;  spc3n = c3n + nel; smc3n = spc3n + nel;

      /* Compute tables */
      for (n = 1; n < m4; n++) {
        if (n == m8) continue;
        ang = n * M_2PI / m;
        c = cos(ang); s = sin(ang);
        *cn++ = c; *spcn++ = - (s + c); *smcn++ = s - c;
        ang = 3 * n * M_2PI / m;
        c = cos(ang); s = sin(ang);
        *c3n++ = c; *spc3n++ = - (s + c); *smc3n++ = s - c;
      }
    }
  }
}

static void complexfft_compute_recursive(
    const pk_srfft_t *self,
    float *xr,
    float *xi,
    int logn) {
  int m, m2, m4, m8, nel, n;
  float *xr1, *xr2, *xi1, *xi2;
  float *cn, *spcn, *smcn, *c3n, *spc3n, *smc3n;
  float tmp1, tmp2;
  float sqhalf = M_SQRT1_2;

  /* Check range of logn */
  assert(logn >= 0 && "Error: logn is out of bounds in SRFFT");

  /* Compute trivial cases */
  if (logn < 3) {
    if (logn == 2) {  /* length m = 4 */
      xr2  = xr + 2;
      xi2  = xi + 2;
      tmp1 = *xr + *xr2;
      *xr2 = *xr - *xr2;
      *xr  = tmp1;
      tmp1 = *xi + *xi2;
      *xi2 = *xi - *xi2;
      *xi  = tmp1;
      xr1  = xr + 1;
      xi1  = xi + 1;
      xr2++;
      xi2++;
      tmp1 = *xr1 + *xr2;
      *xr2 = *xr1 - *xr2;
      *xr1 = tmp1;
      tmp1 = *xi1 + *xi2;
      *xi2 = *xi1 - *xi2;
      *xi1 = tmp1;
      xr2  = xr + 1;
      xi2  = xi + 1;
      tmp1 = *xr + *xr2;
      *xr2 = *xr - *xr2;
      *xr  = tmp1;
      tmp1 = *xi + *xi2;
      *xi2 = *xi - *xi2;
      *xi  = tmp1;
      xr1  = xr + 2;
      xi1  = xi + 2;
      xr2  = xr + 3;
      xi2  = xi + 3;
      tmp1 = *xr1 + *xi2;
      tmp2 = *xi1 + *xr2;
      *xi1 = *xi1 - *xr2;
      *xr2 = *xr1 - *xi2;
      *xr1 = tmp1;
      *xi2 = tmp2;
      return;
    }
    else if (logn == 1) {   /* length m = 2 */
      xr2  = xr + 1;
      xi2  = xi + 1;
      tmp1 = *xr + *xr2;
      *xr2 = *xr - *xr2;
      *xr  = tmp1;
      tmp1 = *xi + *xi2;
      *xi2 = *xi - *xi2;
      *xi  = tmp1;
      return;
    }
    else if (logn == 0) return;   /* length m = 1 */
  }

  /* Compute a few constants */
  m = 1 << logn; m2 = m / 2; m4 = m2 / 2; m8 = m4 /2;


  /* Step 1 */
  xr1 = xr; xr2 = xr1 + m2;
  xi1 = xi; xi2 = xi1 + m2;
  for (n = 0; n < m2; n++) {
    tmp1 = *xr1 + *xr2;
    *xr2 = *xr1 - *xr2;
    xr2++;
    *xr1++ = tmp1;
    tmp2 = *xi1 + *xi2;
    *xi2 = *xi1 - *xi2;
    xi2++;
    *xi1++ = tmp2;
  }

  /* Step 2 */
  xr1 = xr + m2; xr2 = xr1 + m4;
  xi1 = xi + m2; xi2 = xi1 + m4;
  for (n = 0; n < m4; n++) {
    tmp1 = *xr1 + *xi2;
    tmp2 = *xi1 + *xr2;
    *xi1 = *xi1 - *xr2;
    xi1++;
    *xr2++ = *xr1 - *xi2;
    *xr1++ = tmp1;
    *xi2++ = tmp2;
    // xr1++; xr2++; xi1++; xi2++;
  }

  /* Steps 3 & 4 */
  xr1 = xr + m2; xr2 = xr1 + m4;
  xi1 = xi + m2; xi2 = xi1 + m4;
  if (logn >= 4) {
    nel = m4 - 2;
    cn  = self->tab[logn - 4]; spcn  = cn + nel;  smcn  = spcn + nel;
    c3n = smcn + nel;  spc3n = c3n + nel; smc3n = spc3n + nel;
  }
  xr1++; xr2++; xi1++; xi2++;
  // xr1++; xi1++;
  for (n = 1; n < m4; n++) {
    if (n == m8) {
      tmp1 =  sqhalf * (*xr1 + *xi1);
      *xi1 =  sqhalf * (*xi1 - *xr1);
      *xr1 =  tmp1;
      tmp2 =  sqhalf * (*xi2 - *xr2);
      *xi2 = -sqhalf * (*xr2 + *xi2);
      *xr2 =  tmp2;
    } else {
      tmp2 = *cn++ * (*xr1 + *xi1);
      tmp1 = *spcn++ * *xr1 + tmp2;
      *xr1 = *smcn++ * *xi1 + tmp2;
      *xi1 = tmp1;
      tmp2 = *c3n++ * (*xr2 + *xi2);
      tmp1 = *spc3n++ * *xr2 + tmp2;
      *xr2 = *smc3n++ * *xi2 + tmp2;
      *xi2 = tmp1;
    }
    xr1++; xr2++; xi1++; xi2++;
  }

  /* Call ssrec again with half DFT length */
  complexfft_compute_recursive(self, xr, xi, logn - 1);

  /* Call ssrec again twice with one quarter DFT length.
     Constants have to be recomputed, because they are static! */
  // m = 1 << logn; m2 = m / 2;
  complexfft_compute_recursive(self, xr + m2, xi + m2, logn - 2);
  // m = 1 << logn;
  m4 = 3 * (m / 4);
  complexfft_compute_recursive(self, xr + m4, xi + m4, logn - 2);
}

static void 
complexfft_bit_reverse_permute(
    const pk_srfft_t *self,
    float *x,
    int logn)  {
  int i, j, lg2, n;
  int off, fj, gno, *brp;
  float tmp, *xp, *xq;

  lg2 = logn >> 1;
  n = 1 << lg2;
  if (logn & 1) lg2++;

  /* Unshuffling loop */
  for (off = 1; off < n; off++) {
    fj = n * self->brseed[off]; i = off; j = fj;
    tmp = x[i]; x[i] = x[j]; x[j] = tmp;
    xp = &x[i];
    brp = &(self->brseed[1]);
    for (gno = 1; gno < self->brseed[off]; gno++) {
      xp += n;
      j = fj + *brp++;
      xq = x + j;
      tmp = *xp; *xp = *xq; *xq = tmp;
    }
  }
}

static void complexfft_compute2(
    const pk_srfft_t *self,
    float *xr,
    float *xi,
    bool forward) {
  if (!forward) {  // reverse real and imaginary parts for complex FFT.
    float *tmp = xr;
    xr = xi;
    xi = tmp;
  }
  complexfft_compute_recursive(self, xr, xi, self->logn);
  if (self->logn > 1) {
    complexfft_bit_reverse_permute(self, xr, self->logn);
    complexfft_bit_reverse_permute(self, xi, self->logn);
  }
}


static void
complexfft_compute(
    const pk_srfft_t *self,
    float *x,
    int xsize,
    bool forward,
    float *buffer,
    int buffer_size) {
  assert(xsize == self->N * 2 && "complexfft_compute: invalid xsize");
  assert(buffer_size >= self->N && "complexfft_compute: buffer_size too small");

  for (int i = 0; i < self->N; i++) {
    x[i] = x[i * 2];  // put the real part in the first half of x.
    buffer[i] = x[i * 2 + 1];  // put the imaginary part in temp_buffer.
  }
  // copy the imaginary part back to the second half of x.
  memcpy((void *)(x + self->N),
         (void *)buffer,
         sizeof(float) * self->N);

  complexfft_compute2(self, x, x + self->N, forward);
  // Now change the format back to interleaved.
  memcpy((void *)buffer,
         (void *)(x + self->N),
         sizeof(float) * self->N);
  for (int i = self->N - 1; i > 0; i--) {  // don't include 0,
    // in case MatrixIndexT (now it's int) is unsigned, the loop would not
    // terminate. Treat it as a special case.
    x[i * 2] = x[i];
    x[i * 2 + 1] = buffer[i];
  }
  x[1] = buffer[0];  // special case of i = 0.
}

inline void complex_im_exp(float x, float *a_re, float *a_im) {
  *a_re = cos(x);
  *a_im = sin(x);
}

//! ComplexMul implements, inline, the complex multiplication b *= a.
inline void complex_mul(float a_re, float a_im, float *b_re, float *b_im) {
  float tmp_re = (*b_re * a_re) - (*b_im * a_im);
  *b_im = *b_re * a_im + *b_im * a_re;
  *b_re = tmp_re;
}

inline void complex_add_product(
    float a_re,
    float a_im,
    float b_re,
    float b_im,
    float *c_re,
    float *c_im) {
  *c_re += b_re * a_re - b_im * a_im;
  *c_im += b_re * a_im + b_im * a_re;
}


void pk_srfft_init(pk_srfft_t *self, int N) {
  N /= 2;
  if ( (N & (N-1)) != 0 || N <= 1) {
    assert(false && "pk_srfft_init called with invalid number of points.");
  }

  self->N = N;
  self->logn = 0;
  while (N > 1) {
    N >>= 1;
    ++self->logn;
  }
  compute_tables(self);
}

void pk_srfft_destroy(pk_srfft_t *self) {
  free(self->brseed);
  if (self->tab != NULL) {
    for (int i = 0; i < self->logn - 3; i++) {
      free(self->tab[i]);
    }
    free(self->tab);
  }
}


// This code is mostly the same as the RealFft function.  It would be
// possible to replace it with more efficient code from Rico's book.
void pk_srfft_compute(
    const pk_srfft_t *self,
    float *data,
    int data_size,
    bool forward,
    float *buffer,
    int buffer_size) {
  int N = self->N * 2, N2 = self->N;
  if (forward) {
    // call to base class
    complexfft_compute(self, data, data_size, true, buffer, buffer_size);
  } 
    
  // exp(-2pi/N), forward; exp(2pi/N), backward
  float rootN_re, rootN_im;  
  int forward_sign = forward ? -1 : 1;
  complex_im_exp((float)(M_2PI / N * forward_sign), &rootN_re, &rootN_im);
  
  // exp(-2pik/N), forward; exp(-2pik/N), backward
  float kN_re = -forward_sign, kN_im = 0.0;

  // kN starts out as 1.0 for forward algorithm but -1.0 for backward.
  for (int k = 1; 2 * k <= N2; k++) {
    complex_mul(rootN_re, rootN_im, &kN_re, &kN_im);

    float Ck_re, Ck_im, Dk_re, Dk_im;
    // C_k = 1/2 (B_k + B_{N/2 - k}^*) :
    Ck_re = 0.5 * (data[2 * k] + data[N - 2 * k]);
    Ck_im = 0.5 * (data[2 * k + 1] - data[N - 2 * k + 1]);
    // re(D_k)= 1/2 (im(B_k) + im(B_{N/2-k})):
    Dk_re = 0.5 * (data[2 * k + 1] + data[N - 2 * k + 1]);
    // im(D_k) = -1/2 (re(B_k) - re(B_{N/2-k}))
    Dk_im =-0.5 * (data[2 * k] - data[N - 2 * k]);
    // A_k = C_k + 1^(k/N) D_k:
    data[2 * k] = Ck_re;  // A_k <-- C_k
    data[2 * k + 1] = Ck_im;
    // now A_k += D_k 1^(k/N)
    complex_add_product(
        Dk_re,
        Dk_im,
        kN_re,
        kN_im,
        &(data[2 * k]),
        &(data[2 * k + 1]));

    int kdash = N2 - k;
    if (kdash != k) {
      // Next we handle the index k' = N/2 - k.  This is necessary
      // to do now, to avoid invalidating data that we will later need.
      // The quantities C_{k'} and D_{k'} are just the conjugates of C_k
      // and D_k, so the equations are simple modifications of the above,
      // replacing Ck_im and Dk_im with their negatives.
      data[2 * kdash] = Ck_re;  // A_k' <-- C_k'
      data[2 * kdash + 1] = -Ck_im;
      // now A_k' += D_k' 1^(k'/N)
      // We use 1^(k'/N) = 1^((N/2 - k) / N) = 1^(1/2) 1^(-k/N) = -1 * (1^(k/N))^*
      // so it's the same as 1^(k/N) but with the real part negated.
      complex_add_product(
          Dk_re,
          -Dk_im,
          -kN_re,
          kN_im,
          &(data[2 * kdash]),
          &(data[2 * kdash + 1]));
    }
  }

  // Now handle k = 0.
  // In simple terms: after the complex fft, data[0] becomes the sum of real
  // parts input[0], input[2]... and data[1] becomes the sum of imaginary
  // pats input[1], input[3]...
  // "zeroth" [A_0] is just the sum of input[0]+input[1]+input[2]..
  // and "n2th" [A_{N/2}] is input[0]-input[1]+input[2]... .
  float zeroth = data[0] + data[1],
        n2th = data[0] - data[1];
  data[0] = zeroth;
  data[1] = n2th;
  if (!forward) {
    data[0] /= 2;
    data[1] /= 2;
  }
  
  if (!forward) {  // call to base class
    complexfft_compute(self, data, data_size, false, buffer, buffer_size);
    for (int i = 0; i < N; i++)
      data[i] *= 2.0;
    // This is so we get a factor of N increase, rather than N/2 which we would
    // otherwise get from [ComplexFft, forward] + [ComplexFft, backward] in dimension N/2.
    // It's for consistency with our normal FFT convensions.
  }
}
