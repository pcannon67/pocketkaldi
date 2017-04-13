// Created at 2017-03-29

#include "pocketkaldi.h"

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "am.h"
#include "cmvn.h"
#include "decodable.h"
#include "decoder.h"
#include "fbank.h"
#include "fst.h"
#include "list.h"
#include "nnet.h"
#include "symbol_table.h"
#include "transition.h"
#include "pcm_reader.h"

PKLIST_DEFINE(char, byte_list)

// The internal version of an utterance. It stores the intermediate state in
// decoding.
typedef struct pk_utterance_internal_t {
  pk_vector_t raw_wave;
} pk_utterance_internal_t;

void pk_init(pk_t *self) {
  self->fst = NULL;
  self->am = NULL;
  self->cmvn_global_stats = NULL;
  self->trans_model = NULL;
  self->symbol_table = NULL;
  self->fbank = NULL;
}

void pk_destroy(pk_t *self) {
  if (self->fst) {
    pk_fst_destroy(self->fst);
    free(self->fst);
    self->fst = NULL;
  }

  if (self->am) {
    pk_am_destroy(self->am);
    free(self->am);
    self->am = NULL;
  }

  if (self->cmvn_global_stats) {
    pk_vector_destroy(self->cmvn_global_stats);
    free(self->cmvn_global_stats);
    self->cmvn_global_stats = NULL;
  }

  if (self->trans_model) {
    pk_transition_destroy(self->trans_model);
    free(self->trans_model);
    self->trans_model = NULL;
  }

  if (self->symbol_table) {
    pk_symboltable_destroy(self->symbol_table);
    free(self->symbol_table);
    self->symbol_table = NULL;
  }

  if (self->fbank) {
    pk_fbank_destroy(self->fbank);
    free(self->fbank);
    self->fbank = NULL;
  }
}

// Pocketkaldi model struct
//   FST
//   CMVN
//   TRANS_MODEL
//   AM
//   SYMBOL_TABLE
void pk_load(pk_t *self, const char *filename, pk_status_t *status) {
  pk_readable_t *fd = pk_readable_open(filename, status);

  // FST
  self->fst = (pk_fst_t *)malloc(sizeof(pk_fst_t));
  pk_fst_init(self->fst);
  pk_fst_read(self->fst, fd, status);
  if (!status->ok) goto pk_load_failed;

  // CMVN
  self->cmvn_global_stats = (pk_vector_t *)malloc(sizeof(pk_vector_t));
  pk_vector_init(self->cmvn_global_stats, 0, NAN);
  pk_vector_read(self->cmvn_global_stats, fd, status);
  if (!status->ok) goto pk_load_failed;

  // TRANS_MODEL
  self->trans_model = (pk_transition_t *)malloc(sizeof(pk_transition_t));
  pk_transition_init(self->trans_model);
  pk_transition_read(self->trans_model, fd, status);
  if (!status->ok) goto pk_load_failed;

  // AM
  self->am = (pk_am_t *)malloc(sizeof(pk_am_t));
  pk_am_init(self->am);
  pk_am_read(self->am, fd, status);
  if (!status->ok) goto pk_load_failed;

  // SYMBOL TABLE
  self->symbol_table = (pk_symboltable_t *)malloc(sizeof(pk_symboltable_t));
  pk_symboltable_init(self->symbol_table);
  pk_symboltable_read(self->symbol_table, fd, status);
  if (!status->ok) goto pk_load_failed;

  // Initialize fbank feature extractor
  self->fbank = (pk_fbank_t *)malloc(sizeof(pk_fbank_t));
  pk_fbank_init(self->fbank);

  if (false) {
pk_load_failed:
    pk_destroy(self);
  }
  pk_readable_close(fd);
}

void pk_utterance_init(pk_utterance_t *utt) {
  utt->internal = (pk_utterance_internal_t *)malloc(
      sizeof(pk_utterance_internal_t));
  pk_vector_init(&utt->internal->raw_wave, 0, NAN);
  utt->hyp = NULL;
  utt->loglikelihood_per_frame = 0.0f;
}

void pk_utterance_destroy(pk_utterance_t *utt) {
  if (utt->internal) {
    pk_vector_destroy(&utt->internal->raw_wave);

    free(utt->internal);
    utt->internal = NULL;
  }

  free(utt->hyp);
  utt->hyp = NULL;
  utt->loglikelihood_per_frame = 0.0f;
}

void pk_read_audio(
    pk_utterance_t *utt,
    const char *filename,
    pk_status_t *status) {
  assert(utt->internal && "pk_read_audio: utterance is not initialized");

  pk_16kpcm_read(filename, &utt->internal->raw_wave, status);
}

void pk_process(pk_t *recognizer, pk_utterance_t *utt) {
  assert(utt->hyp == NULL && "utt->hyp == NULL expected");

  // Handle empty utterance
  if (utt->internal->raw_wave.dim == 0) {
    utt->hyp = (char *)malloc(sizeof(char));
    *(utt->hyp) = '\0';
    return;
  }

  clock_t t;

  // Extract fbank feats from raw_wave
  t = clock();
  pk_matrix_t raw_feats;
  pk_matrix_init(&raw_feats, 0, 0);
  pk_fbank_compute(recognizer->fbank, &utt->internal->raw_wave, &raw_feats);
  t = clock() - t;
  fprintf(stderr, "Fbank: %lfms\n", ((float)t) / CLOCKS_PER_SEC  * 1000);

  // Apply CMVN to raw_wave
  t = clock();
  pk_cmvn_t cmvn;
  pk_cmvn_init(&cmvn, recognizer->cmvn_global_stats, &raw_feats);
  pk_matrix_t feats;
  pk_matrix_init(&feats, raw_feats.nrow, raw_feats.ncol);
  for (int frame = 0; frame < raw_feats.ncol; ++frame) {
    pk_vector_t frame_col = pk_matrix_getcol(&feats, frame);
    pk_cmvn_getframe(&cmvn, frame, &frame_col);
  }
  t = clock() - t;
  fprintf(stderr, "CMVN: %lfms\n", ((float)t) / CLOCKS_PER_SEC  * 1000);

  // Start to decode
  pk_decoder_t decoder;
  pk_decoder_init(&decoder, recognizer->fst);
  pk_decodable_t decodable;
  t = clock();
  pk_decodable_init(
      &decodable,
      recognizer->am,
      recognizer->trans_model,
      0.1,
      &feats);
  t = clock() - t;
  fprintf(stderr, "NNET: %lfms\n", ((float)t) / CLOCKS_PER_SEC  * 1000);
  
  pk_decoder_decode(&decoder, &decodable);

  // Get final result
  pk_decoder_result_t best_path;
  byte_list_t hyp;
  byte_list_init(&hyp);
  pk_decoder_result_init(&best_path);
  if (pk_decoder_reachedfinal(&decoder) &&
      pk_decoder_bestpath(&decoder, &best_path) &&
      best_path.size > 0) {
    for (int idx = 0; idx < best_path.size; ++idx) {
      int word_id = best_path.words[idx];

      // Append the word into hyp
      const char *word = pk_symboltable_get(recognizer->symbol_table, word_id);
      while (*word) {
        byte_list_push_back(&hyp, *word);
        ++word;
      }
      byte_list_push_back(&hyp, ' ');
    }

    // Change last whitespace into '\0'
    hyp.data[hyp.size - 1] = '\0';

    // Copy hyp to utt->hyp
    utt->hyp = (char *)malloc(sizeof(char) * hyp.size);
    pk_strlcpy(utt->hyp, hyp.data, hyp.size);
    utt->loglikelihood_per_frame = best_path.weight / feats.ncol;
  } else {
    utt->hyp = (char *)malloc(sizeof(char));
    *(utt->hyp) = '\0';
  }

  pk_matrix_destroy(&raw_feats);
  pk_cmvn_destroy(&cmvn);
  pk_matrix_destroy(&feats);
  pk_decoder_destroy(&decoder);
  pk_decodable_destroy(&decodable);
  pk_decoder_result_destroy(&best_path);
  byte_list_destroy(&hyp);
}
