// gmmbin/gmm-decode-simple.cc

// Copyright 2009-2011  Microsoft Corporation
//                2013  Johns Hopkins University

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

#include <assert.h>
#include <math.h>
#include "itf/decodable-itf.h"
#include "base/kaldi-common.h"
#include "hmm/transition-model.h"
#include "util/common-utils.h"
#include "nnet2/decodable-am-nnet.h"
#include "tree/context-dep.h"
#include "hmm/transition-model.h"
#include "fstext/lattice-utils.h"
#include "lat/kaldi-lattice.h"
#include "base/timer.h"
#include "fst.h"
#include "util.h"
#include "decoder.h"
#include "matrix.h"
#include "cmvn.h"
#include "am.h"
#include "util.h"
#include "transition.h"


namespace kaldi {
namespace nnet2 {

/// DecodableAmNnet is a decodable object that decodes
/// with a neural net acoustic model of type AmNnet.
class PkDecodableAmNnet: public DecodableInterface {
 public:
  PkDecodableAmNnet(const std::string &id2pdf_filename,
                    const pk_am_t *am,
                    pk_matrix_t *feats,
                    bool pad_input = true, // if !pad_input, the NumIndices()
                                           // will be < feats.NumRows().
                    BaseFloat prob_scale = 1.0) {
    pk_status_t status;
    pk_status_init(&status);
    pk_readable_t *fd = pk_readable_open(id2pdf_filename.c_str(), &status);
    pk_transition_init(&trans_model_, fd, &status);
    pk_readable_close(fd);

    pk_vector_t y;
    pk_vector_t log_prior;
    pk_vector_init(&y, 0, NAN);
    pk_vector_init(&log_prior, 0, NAN);

    pk_matrix_init(&log_prob_, trans_model_.num_pdfs, feats->ncol);
    for (int frame = 0; frame < feats->ncol; ++frame) {
      pk_am_compute(am, feats, frame, &y);
      pk_vector_scale(&y, prob_scale);
      for (int d = 0; d < y.dim; ++d) {
        log_prob_.data[frame * log_prob_.nrow + d] = y.data[d];
      }
    }
  }

  // Note, frames are numbered from zero.  But transition_id is numbered
  // from one (this routine is called by FSTs).
  virtual BaseFloat LogLikelihood(int32 frame, int32 transition_id) {
    int pdf_id = pk_transition_tid2pdf(&trans_model_, transition_id);
    float val = log_prob_.data[frame * log_prob_.nrow + pdf_id];
    return val;
  }

  virtual int32 NumFramesReady() const { return log_prob_.ncol; }
  
  // Indices are one-based!  This is for compatibility with OpenFst.
  virtual int32 NumIndices() const { return trans_model_.num_transition_ids; }
  
  virtual bool IsLastFrame(int32 frame) const {
    KALDI_ASSERT(frame < NumFramesReady());
    return (frame == NumFramesReady() - 1);
  }

 protected:
  pk_transition_t trans_model_;
  pk_matrix_t log_prob_; 

  KALDI_DISALLOW_COPY_AND_ASSIGN(PkDecodableAmNnet);
};

}
}

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    typedef kaldi::int32 int32;
    using fst::SymbolTable;
    using fst::VectorFst;
    using fst::StdArc;
    using fst::ReadFstKaldi;
    using namespace kaldi::nnet2;

    const char *usage =
        "Decode features using GMM-based model.\n"
        "Viterbi decoding, Only produces linear sequence; any lattice\n"
        "produced is linear\n"
        "\n"
        "Usage:   gmm-decode-simple [options] <model-in> <fst-in> "
        "<features-rspecifier> <words-wspecifier> [<alignments-wspecifier>] "
        "[<lattice-wspecifier>]";
    ParseOptions po(usage);
    Timer timer;
    bool allow_partial = true; 
    BaseFloat acoustic_scale = 0.1;

    std::string word_syms_filename;
    std::string cmvn_stats_filename;
    BaseFloat beam = 16.0;
    po.Register("beam", &beam, "Decoding log-likelihood beam");
    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods");
    po.Register("word-symbol-table", &word_syms_filename,
                "Symbol table for words [for debug output]");
    po.Register("cmvn-stats", &cmvn_stats_filename,
                "Vecotr of cmvn stats");
    po.Read(argc, argv);

    if (po.NumArgs() < 4 || po.NumArgs() > 6) {
      po.PrintUsage();
      exit(1);
    }

    std::string am_filename = po.GetArg(1),
        fst_in_filename = po.GetArg(2),
        feature_rspecifier = po.GetArg(3),
        words_wspecifier = po.GetArg(4),
        transition_in_filename = po.GetOptArg(5);

    pk_fst_t fst;
    pk_status_t status;

    pk_am_t am;
    pk_am_init(&am);
    pk_readable_t *fd_am = pk_readable_open(am_filename.c_str(), &status);
    assert(status.ok);
    pk_am_read(&am, fd_am, &status);
    assert(status.ok);

    pk_status_init(&status);
    pk_fst_read(&fst, fst_in_filename.c_str(), &status);
    if (!status.ok) {
      puts(status.message);
      exit(1);
    }

    Int32VectorWriter words_writer(words_wspecifier);
    fst::SymbolTable *word_syms = NULL;
    if (word_syms_filename != "") 
      if (!(word_syms = fst::SymbolTable::ReadText(word_syms_filename)))
        KALDI_ERR << "Could not read symbol table from file "
                   << word_syms_filename;

    SequentialBaseFloatMatrixReader feature_reader(feature_rspecifier);

    BaseFloat tot_like = 0.0;
    kaldi::int64 frame_count = 0;
    int num_success = 0, num_fail = 0;
    PkSimpleDecoder decoder(&fst, beam);

    // Read CMVN stats
    pk_readable_t *fdcmvn = pk_readable_open(
        cmvn_stats_filename.c_str(),
        &status);
    assert(status.ok);

    pk_vector_t cmvn_stats;
    pk_vector_init(&cmvn_stats, 0, NAN);
    pk_vector_read(&cmvn_stats, fdcmvn, &status);
    assert(status.ok);

    pk_readable_close(fdcmvn);

    for (; !feature_reader.Done(); feature_reader.Next()) {
      std::string utt = feature_reader.Key();
      Matrix<BaseFloat> feature_in(feature_reader.Value());

      // Copy matrix feature to raw_feats
      int nrow = feature_in.NumCols();
      int ncol = feature_in.NumRows();
      pk_matrix_t raw_feats;
      pk_matrix_init(&raw_feats, nrow, ncol);
      for (int frame = 0; frame < ncol; ++frame) {
        pk_vector_t col = pk_matrix_getcol(&raw_feats, frame);
        for (int d = 0; d < col.dim; ++d) {
          col.data[d] = feature_in(frame, d);
        }
      }

      // Apply CMVN on raw_feats
      pk_cmvn_t cmvn;
      pk_cmvn_init(&cmvn, &cmvn_stats, &raw_feats);
      pk_vector_t frame_feats;
      pk_vector_init(&frame_feats, 0, NAN);
      pk_matrix_t feats;
      pk_matrix_init(&feats, raw_feats.nrow, raw_feats.ncol);
      for (int frame = 0; frame < raw_feats.ncol; ++frame) {
        pk_cmvn_getframe(&cmvn, frame, &frame_feats);
        for (int d = 0; d < frame_feats.dim; ++d) {
          feats.data[frame * feats.nrow + d] = frame_feats.data[d];
        }
      }

      pk_matrix_destroy(&raw_feats);
      pk_vector_destroy(&frame_feats);
      pk_cmvn_destroy(&cmvn);

      feature_reader.FreeCurrent();
      bool pad_input = true;
      PkDecodableAmNnet nnet_decodable(transition_in_filename,
                                       &am,
                                       &feats,
                                       pad_input,
                                       acoustic_scale);
      decoder.Decode(&nnet_decodable);

      VectorFst<LatticeArc> decoded;  // linear FST.
      pk_decoder_result_t best_path;

      if ( (allow_partial || decoder.ReachedFinal())
           && decoder.GetBestPath(&best_path) ) {
        if (!decoder.ReachedFinal())
          KALDI_WARN << "Decoder did not reach end-state, "
                     << "outputting partial traceback since --allow-partial=true";
        num_success++;

        std::vector<int32> alignment(
            best_path.alignment,
            best_path.alignment + best_path.alignment_size);
        std::vector<int32> words(
            best_path.words,
            best_path.words + best_path.size);
        float weight = best_path.weight;
        pk_decoder_result_destroy(&best_path);

        frame_count += feats.ncol;
        words_writer.Write(utt, words);
        if (word_syms != NULL) {
          std::cerr << utt << ' ';
          for (size_t i = 0; i < words.size(); i++) {
            std::string s = word_syms->Find(words[i]);
            if (s == "")
              KALDI_ERR << "Word-id " << words[i] <<" not in symbol table.";
            std::cerr << s << ' ';
          }
          std::cerr << '\n';
        }
        BaseFloat like = -weight;
        tot_like += like;
        KALDI_LOG << "Log-like per frame for utterance " << utt << " is "
                  << (like / feats.ncol) << " over "
                  << feats.ncol << " frames.";
      } else {
        num_fail++;
        KALDI_WARN << "Did not successfully decode utterance " << utt
                   << ", len = " << feats.ncol;
      }
    }

    double elapsed = timer.Elapsed();
    KALDI_LOG << "Time taken "<< elapsed
              << "s: real-time factor assuming 100 frames/sec is "
              << (elapsed*100.0/frame_count);
    KALDI_LOG << "Done " << num_success << " utterances, failed for "
              << num_fail;
    KALDI_LOG << "Overall log-likelihood per frame is " << (tot_like/frame_count) << " over "
              << frame_count<<" frames.";

    delete word_syms;
    pk_fst_destroy(&fst);
    if (num_success != 0) return 0;
    else return 1;
  } catch(const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}


