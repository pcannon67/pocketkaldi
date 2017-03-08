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

#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "gmm/am-diag-gmm.h"
#include "tree/context-dep.h"
#include "hmm/transition-model.h"
#include "gmm/decodable-am-diag-gmm.h"
#include "fstext/lattice-utils.h"
#include "lat/kaldi-lattice.h"
#include "base/timer.h"
#include "fst.h"
#include "util.h"
#include "decoder.h"


int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    typedef kaldi::int32 int32;
    using fst::SymbolTable;
    using fst::VectorFst;
    using fst::StdArc;
    using fst::ReadFstKaldi;

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
    BaseFloat beam = 16.0;
    po.Register("beam", &beam, "Decoding log-likelihood beam");
    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods");
    po.Register("word-symbol-table", &word_syms_filename,
                "Symbol table for words [for debug output]");
    po.Register("allow-partial", &allow_partial,
                "Produce output even when final state was not reached");
    po.Read(argc, argv);

    if (po.NumArgs() < 4 || po.NumArgs() > 6) {
      po.PrintUsage();
      exit(1);
    }

    std::string model_in_filename = po.GetArg(1),
        fst_in_filename = po.GetArg(2),
        feature_rspecifier = po.GetArg(3),
        words_wspecifier = po.GetArg(4),
        alignment_wspecifier = po.GetOptArg(5),
        lattice_wspecifier = po.GetOptArg(6);

    TransitionModel trans_model;
    AmDiagGmm am_gmm;
    {
      bool binary;
      Input ki(model_in_filename, &binary);
      trans_model.Read(ki.Stream(), binary);
      am_gmm.Read(ki.Stream(), binary);
    }

    pk_fst_t fst;
    pk_status_t status;

    pk_status_init(&status);
    pk_fst_read(&fst, fst_in_filename.c_str(), &status);
    if (!status.ok) {
      puts(status.message);
      exit(1);
    }

    Int32VectorWriter words_writer(words_wspecifier);

    Int32VectorWriter alignment_writer(alignment_wspecifier);

    CompactLatticeWriter clat_writer(lattice_wspecifier);

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

    for (; !feature_reader.Done(); feature_reader.Next()) {
      std::string utt = feature_reader.Key();
      Matrix<BaseFloat> features (feature_reader.Value());
      feature_reader.FreeCurrent();
      if (features.NumRows() == 0) {
        KALDI_WARN << "Zero-length utterance: " << utt;
        num_fail++;
        continue;
      }

      DecodableAmDiagGmmScaled gmm_decodable(am_gmm, trans_model, features,
                                             acoustic_scale);
      decoder.Decode(&gmm_decodable);

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

        frame_count += features.NumRows();
        words_writer.Write(utt, words);
        if (alignment_wspecifier != "")
          alignment_writer.Write(utt, alignment);
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
                  << (like / features.NumRows()) << " over "
                  << features.NumRows() << " frames.";
      } else {
        num_fail++;
        KALDI_WARN << "Did not successfully decode utterance " << utt
                   << ", len = " << features.NumRows();
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


