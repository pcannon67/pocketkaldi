// 2016-12-16

#ifndef POCKETKALDI_SRFFT_H_
#define POCKETKALDI_SRFFT_H_

#include <stdbool.h>
#include "pocketkaldi.h"

typedef struct pk_srfft_t {
  int N;
  int logn;

  // brseed is Evans' seed table, ref:  (Ref: D. M. W.
  // Evans, "An improved digit-reversal permutation algorithm ...",
  // IEEE Trans. ASSP, Aug. 1987, pp. 1120-1125).
  int *brseed;

  // Tables of butterfly coefficients.
  float **tab;
} pk_srfft_t;

// N is the number of complex points (must be a power of two, or this
// will crash).  Note that the constructor does some work so it's best to
// initialize the object once and do the computation many times.
POCKETKALDI_EXPORT
void pk_srfft_init(pk_srfft_t *self, int N);

// Destroy and free all resources
POCKETKALDI_EXPORT
void pk_srfft_destroy(pk_srfft_t *self);

// This version of Compute is const; it operates on an array of size N*2
// containing [ r0 im0 r1 im1 ... ], but it uses the argument "temp_buffer" as
// temporary storage instead of a class-member variable.  It will allocate it if
// needed.
// `buffer` needs at least N float spaces
POCKETKALDI_EXPORT
void pk_srfft_compute(
    const pk_srfft_t *self,
    float *data,
    int data_size,
    bool forward,
    float *buffer,
    int buffer_size);

#endif  // POCKETKALDI_SRFFT_H_
