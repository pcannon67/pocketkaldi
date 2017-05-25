// 2017-02-03

#include "fbank.h"

#include <assert.h>
#include <math.h>
#include <string>
#include <fstream>
#include <vector>
#include "pcm_reader.h"
#include "util.h"

using pocketkaldi::Fbank;

void TestFbank() {
  std::string wav_file = TESTDIR "data/en-us-hello.wav";
  std::string kaldi_featdump = TESTDIR "data/fbankmat_en-us-hello.wav.txt";

  // Read wave file
  pk_status_t status;
  pk_vector_t pcm_data;

  pk_status_init(&status);
  pk_vector_init(&pcm_data, 0, NAN);
  pk_16kpcm_read(wav_file.c_str(), &pcm_data, &status);
  assert(status.ok);

  // Calculate fbank feature
  Fbank fbank;
  pk_matrix_t fbank_feat;
  pk_matrix_init(&fbank_feat, 0, 0);

  fbank.Compute(&pcm_data, &fbank_feat);
  std::vector<float> fbank_featvec;
  for (int i = 0; i < fbank_feat.ncol; ++i) {
    pk_vector_t col = pk_matrix_getcol(&fbank_feat, i);
    for (int j = 0; j < col.dim; ++j) {
      fbank_featvec.push_back(col.data[j]);
    }
  }

  // Check the fbank_feat
  std::ifstream is(kaldi_featdump);
  assert(is.is_open());

  float val;
  int line_count = 0;
  for (int i = 0; is >> val; ++i) {
    assert(abs(val - fbank_featvec[i]) < 1e-5);
    ++line_count;
  }
  assert(line_count == 1880);

  pk_vector_destroy(&pcm_data);
  pk_matrix_destroy(&fbank_feat);
}


int main() {
  TestFbank();
  return 0;
}