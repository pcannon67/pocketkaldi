// Created at 2017-01-29


#ifndef POCKETKALDI_FBANK_H_
#define POCKETKALDI_FBANK_H_

#define PK_SAMPLERATE 16000
#define PK_FRAMESHIFT_MS 10.0
#define PK_FRAMELENGTH_MS 25.0
#define PK_FBANK_DIM 40
#define PK_FBANK_LOWFREQ 20
#define PK_FBANK_HIGHFREQ (PK_SAMPLERATE / 2)
#define PK_PREEMPH_COEFF 0.97

#include "srfft.h"
#include "matrix.h"
#include "pcm_reader.h"
#include "pocketkaldi.h"

// To calculate mel fiterbanks over the power spectrum
typedef struct pk_melbanks_t {
  pk_vector_t bins[PK_FBANK_DIM];
  int fftbin_offset[PK_FBANK_DIM];
} pk_melbanks_t;

typedef struct pk_fbank_t {
  pk_melbanks_t melbanks;
  pk_srfft_t srfft;
  int frame_length_padded;
  pk_vector_t window_function;
} pk_fbank_t;

// Initialize fbank feature extractor
POCKETKALDI_EXPORT
void pk_fbank_init(pk_fbank_t *self);

// Computes the fbank feature from wave, and then stores the result into matrix
// fbank_feature. Each column in fbank_feature represent a frame.
POCKETKALDI_EXPORT
void pk_fbank_compute(
    pk_fbank_t *self,
    const pk_vector_t *wave,
    pk_matrix_t *fbank_feature);

// Destroy the fbank computer
POCKETKALDI_EXPORT
void pk_fbank_destroy(pk_fbank_t *self);

#endif  // POCKETKALDI_FBANK_H_
