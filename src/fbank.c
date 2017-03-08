// Created at 2017-01-29


#include "fbank.h"

#include <assert.h>
#include <float.h>
#include <stdint.h>
#include <math.h>
#include "pcm_reader.h"
#include "matrix.h"
#include "srfft.h"

#define FRAME_SHIFT (int)(PK_SAMPLERATE * 0.001 * PK_FRAMESHIFT_MS)
#define FRAME_LENGTH (int)(PK_SAMPLERATE * 0.001 * PK_FRAMELENGTH_MS)

#ifndef M_2PI
#define M_2PI 6.28318530718
#endif 

// 3 -> 4, 6 -> 8, 30 -> 32, this function just copied from
// src\base\kaldi-math.cc in Kaldi.
static int32_t round_up_to_nearest_power_of_two(int32_t n) {
  assert(n > 0);
  n--;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  return n + 1;
}

// Calculate the frame number of given audio wave, the same as NumFrames in
// Kaldi when --snip-edges=true
static int calc_num_frames(const pk_vector_t *wave) {
  int num_samples = wave->dim;
  if (num_samples < FRAME_LENGTH) {
    return 0;
  } else {
    return (1 + ((num_samples - FRAME_LENGTH) / FRAME_SHIFT));
  }
}

// Does all the windowing steps after actually extracting the windowed signal
static void
process_window(const pk_vector_t *window_function, pk_vector_t *sub_window) {
  // Subtract mean from waveform on each frame (--remove-dc-offset=true)
  float sum = 0;
  for (int i = 0; i < sub_window->dim; ++i) {
    sum += sub_window->data[i];
  }
  float mean = sum / sub_window->dim;
  for (int i = 0; i < sub_window->dim; ++i) {
    sub_window->data[i] -= mean;
  }

  // Signal preemphasis (--preemphasis-coefficient=0.97)
  for (int i = sub_window->dim - 1; i > 0; --i) {
    sub_window->data[i] -= PK_PREEMPH_COEFF * sub_window->data[i - 1];
  }
  sub_window->data[0] -= PK_PREEMPH_COEFF * sub_window->data[0];

  // Apply hamming window
  assert(window_function->dim == sub_window->dim &&
         "dimension of window_function and sub_window mismatch");
  for (int i = 0; i < window_function->dim; ++i) {
    sub_window->data[i] *= window_function->data[i];
  }
}

// Extract a frame from wave and store into vector window. The frame should be
// initialized and have more dimensions than FRAME_LENGTH (Usually the power of
// 2).
static void extract_window(
    const pk_vector_t *wave,
    int frame_idx,
    const pk_vector_t *window_function,
    pk_vector_t *window) {
  int start_sample = frame_idx * FRAME_SHIFT;
  int end_sample = start_sample + FRAME_LENGTH;

  // Some assertion of index
  assert(start_sample >= 0 && end_sample < wave->dim &&
         "start_sample or end_sample out of range");
  assert(window->dim >= FRAME_LENGTH &&
         "window->dim should be greater than FRAME_LENGTH");

  // Copy samples from wave to window
  pk_vector_copyfrom(window, wave->data + start_sample, FRAME_LENGTH);

  // Add padding zeros into the rest of vector window
  int padding_start = FRAME_LENGTH;
  for (int i = padding_start; i < window->dim; ++i) {
    window->data[i] = 0.0;
  }

  // Process window
  pk_vector_t sub_window = pk_vector_subvector(window, 0, FRAME_LENGTH);
  process_window(window_function, &sub_window);
}

// Calculates the mel scale of freq
static inline float mel_scale(float freq) {
  return 1127.0f * logf(1.0f + freq / 700.0f);
}

// Initialize the melbanks.
static void melbanks_init(pk_melbanks_t *self, int frame_length_padded) {
  assert(PK_FBANK_DIM >= 3 && "must have at least 3 mel bins");
  assert(frame_length_padded % 2 == 0);

  float sample_freq = PK_SAMPLERATE;
  int num_fft_bins = frame_length_padded / 2;

  // fft-bin width (nyquist-freq / num_fft_bins)
  float fft_bin_width = sample_freq / frame_length_padded;
  float mel_low_freq = mel_scale(PK_FBANK_LOWFREQ);
  float mel_high_freq = mel_scale(PK_FBANK_HIGHFREQ);

  // divide by num_bins+1 in next line because of end-effects where the bins
  // spread out to the sides.
  float mel_freq_delta = (mel_high_freq - mel_low_freq) / (PK_FBANK_DIM + 1);

  // Initialize the mel bins
  for (int i = 0; i < PK_FBANK_DIM; ++i) {
    pk_vector_init(&(self->bins[i]), 0);
    self->fftbin_offset[i] = 0;
  }

  // Calculate the fft-bins vector of each bin
  pk_vector_t this_fftbin;
  pk_vector_init(&this_fftbin, num_fft_bins);
  for (int bin_idx = 0; bin_idx < PK_FBANK_DIM; ++bin_idx) {
    float left_mel = mel_low_freq + bin_idx * mel_freq_delta;
    float center_mel = mel_low_freq + (bin_idx + 1) * mel_freq_delta;
    float right_mel = mel_low_freq + (bin_idx + 2) * mel_freq_delta;

    // Clear this_fftbin
    pk_vector_fill(&this_fftbin, 0.0f);

    int first_index = -1;
    int last_index = -1;
    for (int i = 0; i < num_fft_bins; ++i) {
      // Center frequency of this fft bin
      float freq = fft_bin_width * i;
      float mel = mel_scale(freq);

      if (mel > left_mel && mel < right_mel) {
        float weight;
        if (mel <= center_mel) {
          weight = (mel - left_mel) / (center_mel - left_mel);
        } else {
          weight = (right_mel - mel) / (right_mel - center_mel);
        }
        this_fftbin.data[i] = weight;
        if (first_index == -1) {
          first_index = i;
        }
        last_index = i;
      }
    }
    assert(first_index != -1 && last_index > first_index);

    // Ok, copy this_fftbin into self->bins[bin_idx]
    self->fftbin_offset[bin_idx] = first_index;
    int size = last_index + 1 - first_index;
    pk_vector_t *bin = &(self->bins[bin_idx]);
    pk_vector_resize(bin, size);
    pk_vector_copyfrom(bin, this_fftbin.data + first_index, size);
  }

  pk_vector_destroy(&this_fftbin);
}


// Compute Mel energies (note: not log enerties).
// At input, "fft_energies" contains the FFT energies (not log).
static void melbanks_compute(
    const pk_melbanks_t *self,
    const pk_vector_t *power_spectrum,
    pk_vector_t *mel_energies_out) {
  assert(mel_energies_out->dim == PK_FBANK_DIM &&
         "unexpected dim of mel_energies_out");

  for (int i = 0; i < PK_FBANK_DIM; ++i) {
    int offset = self->fftbin_offset[i];
    pk_vector_t fftbin = self->bins[i];
    pk_vector_t sub_power_spectrum = pk_vector_subvector(
        power_spectrum,
        offset,
        fftbin.dim);
    float energy = pk_vector_dot(&fftbin, &sub_power_spectrum);
    mel_energies_out->data[i] = energy;

    assert(energy != NAN);
  }
}

// Destroy the melbanks
static void melbanks_destroy(pk_melbanks_t *self) {
  for (int i = 0; i < PK_FBANK_DIM; ++i) {
    pk_vector_destroy(&(self->bins[i]));
    self->fftbin_offset[i] = 0;
  }
}

// ComputePowerSpectrum converts a complex FFT (as produced by the FFT
// functions in matrix/matrix-functions.h), and converts it into
// a power spectrum.  If the complex FFT is a vector of size n (representing
// half the complex FFT of a real signal of size n, as described there),
// this function computes in the first (n/2) + 1 elements of it, the
// energies of the fft bins from zero to the Nyquist frequency.  Contents of the
// remaining (n/2) - 1 elements are undefined at output.
void compute_power_spectrum(pk_vector_t *window) {
  int dim = window->dim;
  assert(dim > 0 && "compute_power_spectrum: dim should be greater than 0");
  // now we have in waveform, first half of complex spectrum
  // it's stored as [real0, realN/2-1, real1, im1, real2, im2, ...]
  int half_dim = dim / 2;

  // Handle this special case
  float first_energy = window->data[0] * window->data[0];
  float last_energy = window->data[1] * window->data[1];

  for (int i = 1; i < half_dim; ++i) {
    float real = window->data[i * 2]; 
    float im = window->data[i * 2 + 1];
    window->data[i] = real * real + im * im;
  }
  window->data[0] = first_energy;
  window->data[half_dim] = last_energy;
}

// Compute the fbank feature for window, and return as feature. It also needs
// a buffer_vector, which would be used in srfft computation. The size of
// buffer_vector should be no less than frame_length_padded
//
// NOTE: window will be changed during computation
//
static void pk_fbank_computeframe(
    const pk_fbank_t *self,
    pk_vector_t *window,
    pk_vector_t *feature,
    const pk_vector_t *buffer) {
  assert(buffer->dim >= self->frame_length_padded && "buffer too small");
  assert(feature->dim == PK_FBANK_DIM && "feature size mismatch");
  const pk_melbanks_t *melbanks = &(self->melbanks);

  // Compute FFT using split-radix algorithm. It just overwrites the data in
  // window 
  pk_srfft_compute(
      &(self->srfft),
      window->data,
      window->dim,
      true,
      buffer->data,
      buffer->dim);

  // Compute power spectrum
  compute_power_spectrum(window);
  pk_vector_t power_spectrum = pk_vector_subvector(
      window,
      0,
      window-> dim / 2 + 1);

  // Sum with mel fiterbanks over the power spectrum
  melbanks_compute(&(self->melbanks), window, feature);

  // Use log energy here
  for (int i = 0; i < feature->dim; ++i) {
    float val = feature->data[i];
    if (val < FLT_EPSILON) val = FLT_EPSILON;
    feature->data[i] = log(val);
  }
}

// Initialize the hamming window according to FRAME_LENGTH
static void hammingwindow_init(pk_vector_t *window) {
  pk_vector_resize(window, FRAME_LENGTH);
  float a = M_2PI / (FRAME_LENGTH - 1);
  for (int i = 0; i < FRAME_LENGTH; ++i) {
    float i_fl = (float)i;
    window->data[i] = 0.54 - 0.46 * cos(a * i_fl);
  }
}

void pk_fbank_init(pk_fbank_t *self) {
  int frame_length_padded = round_up_to_nearest_power_of_two(FRAME_LENGTH);
  melbanks_init(&(self->melbanks), frame_length_padded);
  pk_srfft_init(&(self->srfft), frame_length_padded);
  self->frame_length_padded = frame_length_padded;

  // Hamming window
  pk_vector_init(&(self->window_function), 0);
  hammingwindow_init(&(self->window_function));
}

void pk_fbank_compute(
    pk_fbank_t *self,
    const pk_vector_t *wave,
    pk_matrix_t *fbank_feature) {
  int num_frames = calc_num_frames(wave);
  if (num_frames == 0) {
    pk_matrix_resize(fbank_feature, 0, 0);
  } else {
    pk_matrix_resize(fbank_feature, PK_FBANK_DIM, num_frames);
  }

  // Extract fbank feature frame by frame
  pk_vector_t window, buffer;
  pk_vector_init(&window, self->frame_length_padded);
  pk_vector_init(&buffer, self->frame_length_padded);
  for (int i = 0; i < num_frames; ++i) {
    extract_window(wave, i, &(self->window_function), &window);
    pk_vector_t frame = pk_matrix_getcol(fbank_feature, i);
    pk_fbank_computeframe(self, &window, &frame, &buffer);
  }

  pk_vector_destroy(&window);
  pk_vector_destroy(&buffer);
}

void pk_fbank_destroy(pk_fbank_t *self) {
  melbanks_destroy(&(self->melbanks));
  pk_srfft_destroy(&(self->srfft));
  self->frame_length_padded = 0;
  pk_vector_destroy(&(self->window_function));
}
