// 2017-03-08

#include "cmvn.h"

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include "matrix.h"
#include "fbank.h"

using pocketkaldi::Fbank;
using pocketkaldi::CMVN;

// Read a matrix from a text file and store them into mat. nrows is the rows of
// mat
std::vector<float> ReadArray(const std::string &filename) {
  FILE *fd = fopen(filename.c_str(), "r");
  assert(fd != NULL);

  float val;
  std::vector<float> matrix_data;
  while (fscanf(fd, "%f", &val) != EOF) {
    matrix_data.push_back(val);
  }
  fclose(fd);
  assert(matrix_data.size() % 40 == 0);

  return matrix_data;
}

void TestOnlineCmvn() {
  std::string global_stats_path = TESTDIR "data/cmvn_stats.bin";
  std::string wav_file = TESTDIR "data/en-us-hello.wav";
  std::vector<float> corr = ReadArray(
      TESTDIR "data/fbankcmvnmat_en-us-hello.wav.txt");

  pk_status_t status;
  pk_status_init(&status);

  // Read audio feats
  pk_vector_t pcm_data;
  pk_status_init(&status);
  pk_vector_init(&pcm_data, 0, NAN);
  pk_16kpcm_read(wav_file.c_str(), &pcm_data, &status);
  assert(status.ok);

  // Computes fbank feature
  Fbank fbank;
  pk_matrix_t fbank_feat;
  pk_matrix_init(&fbank_feat, 0, 0);
  fbank.Compute(&pcm_data, &fbank_feat);

  // Read global stats from file
  pk_readable_t *fd = pk_readable_open(global_stats_path.c_str(), &status);
  assert(status.ok);

  pk_vector_t global_stats;
  pk_vector_init(&global_stats, 0, NAN);
  pk_vector_read(&global_stats, fd, &status);
  assert(status.ok);
  pk_readable_close(fd);

  // Initialize cmvn and cmvn instance
  CMVN cmvn(&global_stats, &fbank_feat);

  // Apply CMVN
  pk_vector_t feats;
  pk_vector_init(&feats, fbank_feat.nrow, NAN);
  for (int i = 0; i < fbank_feat.ncol; ++i) {
    cmvn.GetFrame(i, &feats);
    for (int d = 0; d < feats.dim; ++d) {
      assert(feats.data[d] - corr[i * feats.dim + d] < 0.00001);
    }
  }

  pk_vector_destroy(&pcm_data);
  pk_vector_destroy(&global_stats);
  pk_vector_destroy(&feats);
  pk_matrix_destroy(&fbank_feat);
}

int main() {
  TestOnlineCmvn();
  return 0;
}
