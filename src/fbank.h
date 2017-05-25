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
#include "vector.h"
#include "pcm_reader.h"
#include "pocketkaldi.h"

namespace pocketkaldi {

// To calculate mel fiterbanks over the power spectrum
class Melbanks {
 public:
  Vector<float> bins_[PK_FBANK_DIM];
  int fftbin_offset_[PK_FBANK_DIM];
  
  // Calculates the mel scale of freq
  inline float MelScale(float freq) {
    return 1127.0f * logf(1.0f + freq / 700.0f);
  }

  // Initialize the melbanks
  Melbanks(int frame_length_padded);
  ~Melbanks();

  // Compute Mel energies (note: not log enerties).
  // At input, "fft_energies" contains the FFT energies (not log)
  void Compute(
      const Vector<float> &power_spectrum,
      Vector<float> *mel_energies_out);
};

// Fbank feature extractor
class Fbank {
 public:
  Fbank();
  ~Fbank();

  // Computes the fbank feature from wave, and then stores the result into matrix
  // fbank_feature. Each row in fbank_feature represent a frame.
  void Compute(const pk_vector_t *wave, pk_matrix_t *fbank_feature);

 private:
  int frame_length_padded_;
  Melbanks melbanks_;
  pk_srfft_t srfft_;
  Vector<float> window_function_;

  // Initialize the melbanks and fbank
  void Init();

  // Compute the fbank feature for window, and return as feature. It also needs
  // a buffer_vector, which would be used in srfft computation. The size of
  // buffer_vector should be no less than frame_length_padded
  //
  // NOTE: window will be changed during computation
  //
  void ComputeFrame(
      Vector<float> *window,
      Vector<float> *feature,
      Vector<float> *buffer);

  // 3 -> 4, 6 -> 8, 30 -> 32, this function just copied from
  // src\base\kaldi-math.cc in Kaldi.
  int32_t RoundUpToNearestPowerOfTwo(int32_t n);

  // Calculate the frame number of given audio wave, the same as NumFrames in
  // Kaldi when --snip-edges=true
  int CalcNumFrames(const VectorBase<float> &wave);

  // Does all the windowing steps after actually extracting the windowed signal
  void ProcessWindow(
      const VectorBase<float> &window_function,
      VectorBase<float> *sub_window);

  // Extract a frame from wave and store into vector window. The frame should be
  // initialized and have more dimensions than FRAME_LENGTH (Usually the power of
  // 2).
  void ExtractWindow(
      const VectorBase<float> &wave,
      int frame_idx,
      const VectorBase<float> &window_function,
      Vector<float> *window);

  // ComputePowerSpectrum converts a complex FFT (as produced by the FFT
  // functions in matrix/matrix-functions.h), and converts it into
  // a power spectrum.  If the complex FFT is a vector of size n (representing
  // half the complex FFT of a real signal of size n, as described there),
  // this function computes in the first (n/2) + 1 elements of it, the
  // energies of the fft bins from zero to the Nyquist frequency.  Contents of the
  // remaining (n/2) - 1 elements are undefined at output.
  void ComputePowerSpectrum(Vector<float> *window);

  // Initialize the hamming window according to FRAME_LENGTH
  void HammingWindowInit(Vector<float> *window);
};

}  // namespace pocketkaldi

#endif  // POCKETKALDI_FBANK_H_
