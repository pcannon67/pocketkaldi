// Created at 2017-01-29


#include "fbank.h"

#include <assert.h>
#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "pcm_reader.h"
#include "matrix.h"
#include "srfft.h"

#define FRAME_SHIFT (int)(PK_SAMPLERATE * 0.001 * PK_FRAMESHIFT_MS)
#define FRAME_LENGTH (int)(PK_SAMPLERATE * 0.001 * PK_FRAMELENGTH_MS)

#ifndef M_2PI
#define M_2PI 6.28318530718
#endif 

namespace pocketkaldi {

int32_t Fbank::RoundUpToNearestPowerOfTwo(int32_t n) {
  assert(n > 0);
  n--;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  return n + 1;
}

int Fbank::CalcNumFrames(const VectorBase<float> &wave) {
  int num_samples = wave.Dim();
  if (num_samples < FRAME_LENGTH) {
    return 0;
  } else {
    return (1 + ((num_samples - FRAME_LENGTH) / FRAME_SHIFT));
  }
}

void Fbank::ProcessWindow(
    const VectorBase<float> &window_function,
    VectorBase<float> *sub_window) {
  // Subtract mean from waveform on each frame (--remove-dc-offset=true)
  float sum = 0;
  for (int i = 0; i < sub_window->Dim(); ++i) {
    sum += (*sub_window)(i);
  }
  float mean = sum / sub_window->Dim();
  for (int i = 0; i < sub_window->Dim(); ++i) {
    (*sub_window)(i) -= mean;
  }

  // Signal preemphasis (--preemphasis-coefficient=0.97)
  for (int i = sub_window->Dim() - 1; i > 0; --i) {
    (*sub_window)(i) -= PK_PREEMPH_COEFF * (*sub_window)(i - 1);
  }
  (*sub_window)(0) -= PK_PREEMPH_COEFF * (*sub_window)(0);

  // Apply hamming window
  assert(window_function.Dim() == sub_window->Dim() &&
         "dimension of window_function and sub_window mismatch");
  for (int i = 0; i < window_function.Dim(); ++i) {
    (*sub_window)(i) *= window_function(i);
  }
}

// Extract a frame from wave and store into vector window. The frame should be
// initialized and have more dimensions than FRAME_LENGTH (Usually the power of
// 2).
void Fbank::ExtractWindow(
    const VectorBase<float> &wave,
    int frame_idx,
    const VectorBase<float> &window_function,
    Vector<float> *window) {
  int start_sample = frame_idx * FRAME_SHIFT;
  int end_sample = start_sample + FRAME_LENGTH;

  // Some assertion of index
  assert(start_sample >= 0 && end_sample <= wave.Dim() &&
         "start_sample or end_sample out of range");
  assert(window->Dim() >= FRAME_LENGTH &&
         "window->Dim() should be greater than FRAME_LENGTH");

  // Copy samples from wave to window
  SubVector<float> sub_window = window->Range(0, FRAME_LENGTH);
  sub_window.CopyFromVec(wave.Range(start_sample, FRAME_LENGTH));

  // Add padding zeros into the rest of vector window
  int padding_start = FRAME_LENGTH;
  for (int i = padding_start; i < window->Dim(); ++i) {
    (*window)(i) = 0.0;
  }

  // Process window
  ProcessWindow(window_function, &sub_window);
}

// Initialize the melbanks.
Melbanks::Melbanks(int frame_length_padded) {
  assert(PK_FBANK_DIM >= 3 && "must have at least 3 mel bins");
  assert(frame_length_padded % 2 == 0);

  float sample_freq = PK_SAMPLERATE;
  int num_fft_bins = frame_length_padded / 2;

  // fft-bin width (nyquist-freq / num_fft_bins)
  float fft_bin_width = sample_freq / frame_length_padded;
  float mel_low_freq = MelScale(PK_FBANK_LOWFREQ);
  float mel_high_freq = MelScale(PK_FBANK_HIGHFREQ);

  // divide by num_bins+1 in next line because of end-effects where the bins
  // spread out to the sides.
  float mel_freq_delta = (mel_high_freq - mel_low_freq) / (PK_FBANK_DIM + 1);

  // Initialize the mel bins
  for (int i = 0; i < PK_FBANK_DIM; ++i) {
    fftbin_offset_[i] = 0;
  }

  // Calculate the fft-bins vector of each bin
  Vector<float> this_fftbin(num_fft_bins);
  for (int bin_idx = 0; bin_idx < PK_FBANK_DIM; ++bin_idx) {
    float left_mel = mel_low_freq + bin_idx * mel_freq_delta;
    float center_mel = mel_low_freq + (bin_idx + 1) * mel_freq_delta;
    float right_mel = mel_low_freq + (bin_idx + 2) * mel_freq_delta;

    // Clear this_fftbin
    this_fftbin.Set(0.0f);

    int first_index = -1;
    int last_index = -1;
    for (int i = 0; i < num_fft_bins; ++i) {
      // Center frequency of this fft bin
      float freq = fft_bin_width * i;
      float mel = MelScale(freq);

      if (mel > left_mel && mel < right_mel) {
        float weight;
        if (mel <= center_mel) {
          weight = (mel - left_mel) / (center_mel - left_mel);
        } else {
          weight = (right_mel - mel) / (right_mel - center_mel);
        }
        this_fftbin(i) = weight;
        if (first_index == -1) {
          first_index = i;
        }
        last_index = i;
      }
    }
    assert(first_index != -1 && last_index > first_index);

    // Ok, copy this_fftbin into self->bins[bin_idx]
    fftbin_offset_[bin_idx] = first_index;
    int size = last_index + 1 - first_index;
    bins_[bin_idx].Resize(size);
    bins_[bin_idx].CopyFromVec(this_fftbin.Range(first_index, size));
  }
}

void Melbanks::Compute(
    const Vector<float> &power_spectrum,
    Vector<float> *mel_energies_out) {
  assert(mel_energies_out->Dim() == PK_FBANK_DIM &&
         "unexpected dim of mel_energies_out");

  for (int i = 0; i < PK_FBANK_DIM; ++i) {
    int offset = fftbin_offset_[i];
    const Vector<float> &fftbin = bins_[i];
    SubVector<float> sub_power_spectrum = power_spectrum.Range(
        offset,
        fftbin.Dim());
    float energy = fftbin.VecVec(sub_power_spectrum);
    (*mel_energies_out)(i) = energy;

    if (isnan(energy)) {
      PK_WARN("some dimension of mel_energies_out is NAN");
    }
  }
}

// Destroy the melbanks
Melbanks::~Melbanks() {
  for (int i = 0; i < PK_FBANK_DIM; ++i) {
    fftbin_offset_[i] = 0;
  }
}

void Fbank::ComputePowerSpectrum(Vector<float> *window) {
  int dim = window->Dim();
  assert(dim > 0 && "compute_power_spectrum: dim should be greater than 0");
  // now we have in waveform, first half of complex spectrum
  // it's stored as [real0, realN/2-1, real1, im1, real2, im2, ...]
  int half_dim = dim / 2;

  // Handle this special case
  float first_energy = (*window)(0) * (*window)(0);
  float last_energy = (*window)(1) * (*window)(1);

  for (int i = 1; i < half_dim; ++i) {
    float real = (*window)(i * 2);
    float im = (*window)(i * 2 + 1);
    (*window)(i) = real * real + im * im;
  }
  (*window)(0) = first_energy;
  (*window)(half_dim) = last_energy;
}

// Compute the fbank feature for window, and return as feature. It also needs
// a buffer_vector, which would be used in srfft computation. The size of
// buffer_vector should be no less than frame_length_padded
//
// NOTE: window will be changed during computation
//
void Fbank::ComputeFrame(
    Vector<float> *window,
    Vector<float> *feature,
    Vector<float> *buffer) {
  assert(buffer->Dim() >= frame_length_padded_ && "buffer too small");
  assert(feature->Dim() == PK_FBANK_DIM && "feature size mismatch");

  // Compute FFT using split-radix algorithm. It just overwrites the data in
  // window 
  PK_DEBUG(util::Format("window->Dim() = {}", window->Dim()));
  pk_srfft_compute(
      &srfft_,
      window->Data(),
      window->Dim(),
      true,
      buffer->Data(),
      buffer->Dim());

  // Compute power spectrum
  ComputePowerSpectrum(window);
  SubVector<float> power_spectrum = window->Range(0, window->Dim() / 2 + 1);

  // Sum with mel fiterbanks over the power spectrum
  melbanks_.Compute(*window, feature);

  // Use log energy here
  feature->ApplyFloor(FLT_EPSILON);
  feature->ApplyLog();
}

// Initialize the hamming window according to FRAME_LENGTH
void Fbank::HammingWindowInit(Vector<float> *window) {
  window->Resize(FRAME_LENGTH);
  float a = M_2PI / (FRAME_LENGTH - 1);
  for (int i = 0; i < FRAME_LENGTH; ++i) {
    float i_fl = (float)i;
    (*window)(i) = 0.54 - 0.46 * cos(a * i_fl);
  }
}

Fbank::Fbank() :
    frame_length_padded_(RoundUpToNearestPowerOfTwo(FRAME_LENGTH)),
    melbanks_(frame_length_padded_) {
  pk_srfft_init(&srfft_, frame_length_padded_);

  // Hamming window
  HammingWindowInit(&window_function_);
}

void Fbank::Compute(
    const pk_vector_t *wave,
    pk_matrix_t *fbank_feature) {
  SubVector<float> wave_data(wave->data, wave->dim);
  int num_frames = CalcNumFrames(wave_data);
  if (num_frames == 0) {
    pk_matrix_resize(fbank_feature, 0, 0);
  } else {
    pk_matrix_resize(fbank_feature, PK_FBANK_DIM, num_frames);
  }

  // Extract fbank feature frame by frame
  Vector<float> window(frame_length_padded_, Vector<float>::kUndefined),
                buffer(frame_length_padded_, Vector<float>::kUndefined),
                frame_feat(PK_FBANK_DIM, Vector<float>::kUndefined);
  for (int i = 0; i < num_frames; ++i) {
    ExtractWindow(wave_data, i, window_function_, &window);
    ComputeFrame(&window, &frame_feat, &buffer);

    // Store the value back to frame
    pk_vector_t frame = pk_matrix_getcol(fbank_feature, i);
    for (int i = 0; i < frame_feat.Dim(); ++i) {
      frame.data[i] = frame_feat(i);
    }
  }
}

Fbank::~Fbank() {
  pk_srfft_destroy(&srfft_);
  frame_length_padded_ = 0;
}

}  //namespace pocketkaldi
