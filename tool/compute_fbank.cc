// featbin/compute-fbank-feats.cc

// Copyright 2009-2012  Microsoft Corporation
//                      Johns Hopkins University (author: Daniel Povey)

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

#define HAVE_ATLAS 1

#include <stdio.h>
#include "src/fbank.h"
#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "feat/wave-reader.h"

int main(int argc, char **argv) {
  using namespace kaldi;

  const char *usage =
      "Create Mel-filter bank (FBANK) feature files.\n"
      "Usage:  compute-fbank-feats [options...] <wav-rspecifier> <feats-wspecifier>\n";

  // construct all the global objects
  ParseOptions po(usage);
  int32 channel = -1;

  // Register the options
  po.Register("channel", &channel, "Channel to extract (-1 -> expect mono, 0 -> left, 1 -> right)");
  po.Read(argc, argv);
  
  std::string wav_rspecifier = po.GetArg(1);
  std::string output_wspecifier = po.GetArg(2);

  // Initialize the fbank computer
  pk_fbank_t fbank;
  pk_fbank_init(&fbank);

  SequentialTableReader<WaveHolder> reader(wav_rspecifier);
  BaseFloatMatrixWriter kaldi_writer;  // typedef to TableWriter<something>.

  if (!kaldi_writer.Open(output_wspecifier)) {
    KALDI_ERR << "Could not initialize output with wspecifier "
              << output_wspecifier;
  }

  int32 num_utts = 0, num_success = 0;
  for (; !reader.Done(); reader.Next()) {
    num_utts++;
    std::string utt = reader.Key();
    const WaveData &wave_data = reader.Value();
    int32 num_chan = wave_data.Data().NumRows(), this_chan = channel;
    
    // This block works out the channel (0=left, 1=right...)
    KALDI_ASSERT(num_chan > 0);  // should have been caught in
    // reading code if no channels.
    if (channel == -1) {
      this_chan = 0;
      if (num_chan != 1)
        KALDI_WARN << "Channel not specified but you have data with "
                   << num_chan  << " channels; defaulting to zero";
    } else {
      if (this_chan >= num_chan) {
        KALDI_WARN << "File with id " << utt << " has "
                   << num_chan << " channels but you specified channel "
                   << channel << ", producing no output.";
        continue;
      }
    }

    if (PK_SAMPLERATE != wave_data.SampFreq()) {
      KALDI_ERR << "Sample frequency mismatch: PK_SAMPLERATE is "
                << PK_SAMPLERATE << " but data has "
                << wave_data.SampFreq() << ". Utterance is " << utt;
    }

    SubVector<BaseFloat> waveform(wave_data.Data(), this_chan);

    // Store the waveform into pk_vector_t type
    pk_vector_t wave;
    pk_vector_init(&wave, waveform.Dim(), NAN);
    for (int i = 0; i < waveform.Dim(); ++i) {
      wave.data[i] = waveform(i);
    }

    pk_matrix_t feat;
    pk_matrix_init(&feat, 0, 0);
    
    // Compute the fbank feature
    pk_fbank_compute(&fbank, &wave, &feat);

    // Store result from feat to features. pk_matrix_t in pocketkaldi is
    // stored by column while Matrix in kaldi is stored by row
    Matrix<BaseFloat> features;
    features.Resize(feat.ncol, feat.nrow);
    for (int col = 0; col < feat.ncol; ++col) {
      pk_vector_t column = pk_matrix_getcol(&feat, col);
      for (int row = 0; row < column.dim; ++row) {
        features(col, row) = column.data[row];
      }
    }

    pk_vector_destroy(&wave);
    pk_matrix_destroy(&feat);

    kaldi_writer.Write(utt, features);

    if (num_utts % 10 == 0)
      KALDI_LOG << "Processed " << num_utts << " utterances";
    KALDI_VLOG(2) << "Processed features for key " << utt;
    num_success++;
  }
  KALDI_LOG << " Done " << num_success << " out of " << num_utts
            << " utterances.";

  pk_fbank_destroy(&fbank);

  return 0;
}

