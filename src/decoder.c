// Created at 2017-04-12

#include "decoder.h"

#include <time.h>
#include "list.h"
#include "hashtable.h"

#define OLABEL_BEGINIDX -1
#define CUTOFF_SAMPLES 200
#define CUTOFF_SAMPLESEED 0x322

// olabel_t is the struct records the list of output labels of a tok in beam. 
// prev_idx pointes to the previous olabel_t like a link list. And the prev_idx
// of root node is OLABEL_BEGINIDX
typedef struct olabel_t {
  int prev_idx;
  int olabel;
} olabel_t;

// token_t represents a state in the viterbi lattice. olabel_idx is the index
// of its corresponded outpu label link-list in the list impl->olabels
typedef struct token_t {
  int olabel_idx;
  int state;
  float cost;
} tok_t;

PKLIST_DEFINE(tok_t, toklist);
PKLIST_DEFINE(olabel_t, olabellist);
PKLIST_DEFINE(float, floatlist);
PKLIST_DEFINE(int, intlist);

// A frame in the viterbi lattice. toks is a list of states and tok_idx records
// the mapping of state_id of WFST to index in toks
typedef struct beam_t {
  toklist_t toks;
  pk_hashtable_t tok_idx;
} beam_t;

typedef struct pk_decoder_impl_t {
  float original_adaptive_beam;
  olabellist_t olabels;

  // Only used in get_cutoff()
  floatlist_t costs;

  const pk_fst_t *fst;
  int32_t num_frames_decoded;
  beam_t beams[2];
  beam_t *beam;
  beam_t *prev_beam;
} pk_decoder_impl_t;


// Compare two float numbers pointered by a and b
static int float_cmp(const void *a, const void *b, void *thunk) {
  PK_UNUSED(thunk);

  float float_a = *((float *)a);
  float float_b = *((float *)b);
  if (float_a > float_b) {
    return 1;
  } else if (float_a < float_b) {
    return -1;
  } else {
    return 0;
  }
}

// Gets the weight cutoff. We won't go throuth all the toks in beam here to
// calculate cutoff, Since it takes a long time. Instead, we randomly sample N
// costs from the beam and GUESS the cutoff value.
static double get_cutoff(
    pk_decoder_impl_t *self,
    toklist_t *toks,
    float *adaptive_beam,
    tok_t **best_tok) {
  double best_cost = INFINITY;
  floatlist_clear(&self->costs);
  uint64_t next_random = CUTOFF_SAMPLESEED;

  // Probability of sample a cost into self->costs
  float sample_prob = CUTOFF_SAMPLES / (float)toks->size;
  
  for (int i = 0; i < toks->size; ++i) {
    tok_t *tok = toks->data + i;

    // Random sample costs from beam. To be consistant even in multi-thread
    // environment, do don't use rand() here
    next_random = next_random * (uint64_t)25214903917 + 11;
    float random_f = (next_random & 0xffff) / (float)65535;
    if (random_f < sample_prob) {
      floatlist_push_back(&self->costs, tok->cost);
    }

    if (tok->cost < best_cost) {
      best_cost = tok->cost;
      if (best_tok) *best_tok = tok;
    }
  }

  double beam_cutoff = best_cost + self->original_adaptive_beam;
  double max_active_cutoff = NAN;
  
  if (toks->size > PK_DECODER_BEAMSIZE) {
    int cutoff_idx = self->costs.size * PK_DECODER_BEAMSIZE / toks->size;

    // The same as std::nth_element
    pk_introselect_r(
        self->costs.data,
        self->costs.size,
        sizeof(float),
        cutoff_idx,
        float_cmp,
        NULL);
    max_active_cutoff = self->costs.data[cutoff_idx];
  }

  if (max_active_cutoff < beam_cutoff) { // max_active is tighter than beam.
    *adaptive_beam = max_active_cutoff - best_cost + PK_DECODER_BEAMDELTA;
    beam_cutoff = max_active_cutoff;
  } else {
    *adaptive_beam = self->original_adaptive_beam;
  }

  return beam_cutoff;
}

// Initialize the beam
static void beam_init(beam_t *self) {
  toklist_init(&self->toks);
  toklist_reserve(&self->toks, PK_DECODER_BEAMSIZE * 2);
  pk_hashtable_init(&self->tok_idx, PK_DECODER_BEAMSIZE * 4);
}

// Clear all elements in the beam
static void beam_clear(beam_t *self) {
  toklist_clear(&self->toks);
  pk_hashtable_clear(&self->tok_idx);
}

// Destroy the beam
static void beam_destroy(beam_t *self) {
  toklist_destroy(&self->toks);
  pk_hashtable_destroy(&self->tok_idx);
}

// Insert tok into self->beam. The tok is created according to arc. And it will
// either insert a new token or update existing token in the beam.
// from_tok is the previous token of it in viterbi
// total_cost is the cost of new token
// Return true if successfully inserted. Otherwise, when the cost of existing
// tok is less than new one, return false
static bool insert_tok(
    pk_decoder_impl_t *self,
    const pk_fst_arc_t *arc,
    const tok_t *from_tok,
    float total_cost) {
  int next_state = arc->next_state;
  int tok_idx = pk_hashtable_find(&self->beam->tok_idx, next_state, -1);

  // We must get the value of from_tok->olabel_idx here. Since the memory
  // from_tok points to may be realloc-ed by the following
  // toklist_emplace_back() call
  int from_olabelidx = from_tok ? from_tok->olabel_idx : OLABEL_BEGINIDX;
  
  // Create the olabel for next tok when the output_label of arc
  // is not 0 (epsilon)
  int olabel_idx = OLABEL_BEGINIDX;
  if (arc->output_label != 0) {
    olabel_idx = self->olabels.size;
    olabel_t *olabel = olabellist_emplace_back(&self->olabels);
    olabel->prev_idx = from_tok->olabel_idx;
    olabel->olabel = arc->output_label;
  }

  // Insert new or update existing token in the beam
  tok_t *next_tok = NULL;
  if (tok_idx < 0) {
    int num_toks = self->beam->toks.size;
    next_tok = toklist_emplace_back(&self->beam->toks);
    pk_hashtable_upsert(&self->beam->tok_idx, next_state, num_toks);
  } else {
    next_tok = self->beam->toks.data + tok_idx;

    // If the cost of existing token is less than the new one, just discard
    // inserting and return false
    if (next_tok->cost < total_cost) return false;
  }

  next_tok->cost = total_cost;
  next_tok->state = next_state;
  if (olabel_idx != OLABEL_BEGINIDX) {
    next_tok->olabel_idx = olabel_idx;
  } else {
    next_tok->olabel_idx = from_tok ? from_olabelidx : OLABEL_BEGINIDX;
  }

  return true;
}

// Process the emitting (non-epsilon) arcs of each states in the beam
static float process_emitting(
    pk_decoder_impl_t *self,
    pk_decodable_t *decodable) {
  // Clear the prev_beam
  beam_clear(self->prev_beam);

  // Move beam to prev_beam by swapping them
  beam_t *t = self->prev_beam;
  self->prev_beam = self->beam;
  self->beam = t;

  // Calculate beam_cutoff of beam
  float adaptive_beam = INFINITY;
  tok_t *best_tok = NULL;
  float weight_cutoff = get_cutoff(
      self,
      &self->prev_beam->toks,
      &adaptive_beam,
      &best_tok);

  // This is the cutoff we use after adding in the log-likes (i.e.
  // for the next frame).  This is a bound on the cutoff we will use
  // on the next frame.
  double next_weight_cutoff = INFINITY;

  // First process the best token to get a hopefully
  // reasonably tight bound on the next cutoff.
  if (best_tok) {
    float state = best_tok->state;
    pk_fst_iter_t arc_iter;
    pk_fst_iterate_arc(self->fst, state, &arc_iter);
    const pk_fst_arc_t *arc = NULL;
    while ((arc = pk_fst_iter_next(&arc_iter)) != NULL) {
      if (arc->input_label != 0) {
        float acoustic_cost = -pk_decodable_loglikelihood(
            decodable,
            self->num_frames_decoded,
            arc->input_label);
        double total_cost = best_tok->cost + arc->weight + acoustic_cost;
        if (total_cost + adaptive_beam < next_weight_cutoff) {
          next_weight_cutoff = total_cost + adaptive_beam;
        }
      }
    }
  }
  assert(next_weight_cutoff != INFINITY);
  
  int processed_tok = 0;
  toklist_t *prev_toks = &self->prev_beam->toks;
  for (int from_tokidx = 0; from_tokidx < prev_toks->size; ++from_tokidx) {
    tok_t from_tok = prev_toks->data[from_tokidx];
    int state = from_tok.state;

    // weight_cutoff is computed according to beam size
    // So there are only top beam_size toks less than weight_cutoff
    if (from_tok.cost > weight_cutoff) continue;
    ++processed_tok;

    pk_fst_iter_t arc_iter;
    pk_fst_iterate_arc(self->fst, state, &arc_iter);
    const pk_fst_arc_t *arc = NULL;
    while ((arc = pk_fst_iter_next(&arc_iter)) != NULL) {
      if (arc->input_label == 0) continue;

      float ac_cost = -pk_decodable_loglikelihood(
          decodable,
          self->num_frames_decoded,
          arc->input_label);
      double total_cost = from_tok.cost + arc->weight + ac_cost;
      
      // Prune the toks whose cost is too high
      if (total_cost > next_weight_cutoff) continue;
      if (total_cost + adaptive_beam < next_weight_cutoff) {
        next_weight_cutoff = total_cost + adaptive_beam;
      }

      // Create and insert the tok into self->beam
      insert_tok(self, arc, &from_tok, total_cost);
    }
  }
  self->num_frames_decoded++;
  return next_weight_cutoff;
}

// Processes nonemitting arcs for one frame.  Propagates within cur_toks_.
static void process_nonemitting(pk_decoder_impl_t *self, double cutoff) {
  intlist_t queue;
  intlist_init(&queue);

  toklist_t *cur_toks = &self->beam->toks;
  for (int i = 0; i < cur_toks->size; ++i) {
    tok_t *tok = &cur_toks->data[i];
    intlist_push_back(&queue, tok->state);
  }
  
  // Loop until no state in beam have out-going epsilon arc
  while (!intlist_empty(&queue)) {
    int state = intlist_back(&queue);
    intlist_pop_back(&queue);

    // Get tok by state
    int tok_idx = pk_hashtable_find(&self->beam->tok_idx, state, -1);
    assert(tok_idx >= 0);

    pk_fst_iter_t arc_iter;
    pk_fst_iterate_arc(self->fst, state, &arc_iter);
    const pk_fst_arc_t *arc = NULL;
    while ((arc = pk_fst_iter_next(&arc_iter)) != NULL) {
      // propagate nonemitting only...
      if (arc->input_label != 0) continue;

      const float ac_cost = 0.0;
      const tok_t *from_tok = &cur_toks->data[tok_idx];
      double total_cost = from_tok->cost + arc->weight + ac_cost;
      if (total_cost > cutoff) {
        continue;
      }

      // Create and insert tok into beam
      // If the token successfully inserted or updated in the beam, 
      // inserted will be true and then we will push the new state into queue
      bool inserted = insert_tok(self, arc, from_tok, total_cost);
      if (inserted) intlist_push_back(&queue, arc->next_state);
    }
  }

  intlist_destroy(&queue);
}

// Initializes the decoding
static void init_decoding(pk_decoder_impl_t *self) {
  // Prepare beams
  self->prev_beam = self->beams;
  self->beam = self->beams + 1;
  beam_clear(self->prev_beam);
  beam_clear(self->beam);

  // Initialize decoding:
  int start_state = pk_fst_startstate(self->fst);
  assert(start_state >= 0);
  pk_fst_arc_t dummy_arc;
  dummy_arc.input_label = 0;
  dummy_arc.output_label = 0;
  dummy_arc.next_state = start_state;
  dummy_arc.weight = 0;

  insert_tok(self, &dummy_arc, NULL, 0.0f);
  self->num_frames_decoded = 0;
  process_nonemitting(self, INFINITY);
}

void pk_decoder_init(pk_decoder_t *self, const pk_fst_t *fst) {
  pk_decoder_impl_t *impl = (pk_decoder_impl_t *)malloc(
      sizeof(pk_decoder_impl_t));
  impl->original_adaptive_beam = 16.0;
  olabellist_init(&impl->olabels);
  floatlist_init(&impl->costs);
  impl->fst = fst;
  impl->num_frames_decoded = 0;
  beam_init(impl->beams);
  beam_init(impl->beams + 1);
  impl->beam = NULL;
  impl->prev_beam = NULL;
  self->impl = impl;
}

void pk_decoder_destroy(pk_decoder_t *self) {
  pk_decoder_impl_t *impl = (pk_decoder_impl_t *)self->impl;
  olabellist_destroy(&impl->olabels);
  floatlist_destroy(&impl->costs);
  beam_destroy(impl->beams);
  beam_destroy(impl->beams + 1);
  impl->beam = NULL;
  impl->prev_beam = NULL;
  impl->original_adaptive_beam = 0.0f;
  impl->fst = NULL;
  impl->num_frames_decoded = 0;
  free(impl);
  self->impl = NULL;
}

bool pk_decoder_decode(pk_decoder_t *self, pk_decodable_t *decodable) {
  pk_decoder_impl_t *impl = (pk_decoder_impl_t *)self->impl;
  clock_t t;

  // Extract fbank feats from raw_wave
  clock_t t_all = clock();
  clock_t t_emitting = 0;
  clock_t t_nonemitting = 0;
  
  init_decoding(impl);
  while(!pk_decodable_islastframe(decodable, impl->num_frames_decoded - 1)) {
    t = clock();
    double cutoff = process_emitting(impl, decodable);
    t_emitting += clock() - t;
    
    t = clock();
    process_nonemitting(impl, cutoff);
    t_nonemitting = clock() - t;
  }

  t_all = clock() - t_all;
  fprintf(stderr, "decode: %lfms\n", ((float)t_all) / CLOCKS_PER_SEC  * 1000);
  fprintf(stderr, "  process_emitting: %lfms\n", ((float)t_emitting) / CLOCKS_PER_SEC  * 1000);
  fprintf(stderr, "  process_nonemitting: %lfms\n", ((float)t_nonemitting) / CLOCKS_PER_SEC  * 1000);
  return impl->beam->toks.size > 0;
}

bool pk_decoder_reachedfinal(pk_decoder_t *self) {
  return true;
}

// Outputs best_path corresponding to the single best path through the lattice.
bool pk_decoder_bestpath(
    pk_decoder_t *self,
    pk_decoder_result_t *best_path) {
  pk_decoder_impl_t *impl = (pk_decoder_impl_t *)self->impl;
  
  // Find the best token
  tok_t *best_tok = NULL;
  double best_cost = INFINITY;

  toklist_t *cnt_toks = &impl->beam->toks;
  for (int i = 0; i < cnt_toks->size; ++i) {
    tok_t *tok = &cnt_toks->data[i];
    int state = tok->state;
    double this_cost = tok->cost + pk_fst_final(impl->fst, state);
    if (this_cost != INFINITY && this_cost < best_cost) {
      best_cost = this_cost;
      best_tok = tok;
    }
  }
  if (best_tok == NULL) return false;  // No output.

  // Get all output labels from best_tok
  int best_olabelidx = best_tok->olabel_idx;
  int olabel_idx = best_olabelidx;
  intlist_t words;
  intlist_init(&words);
  while (olabel_idx != OLABEL_BEGINIDX) {
    const olabel_t *olabel = &impl->olabels.data[olabel_idx];
    intlist_push_back(&words, olabel->olabel);
    olabel_idx = olabel->prev_idx;
  }

  best_path->words = (int32_t *)malloc(words.size * sizeof(int32_t));
  for (int idx = 0; idx < words.size; ++idx) {
    best_path->words[words.size - idx - 1] = words.data[idx];
  }
  best_path->size = words.size;
  best_path->weight = best_cost;
  best_path->weight += pk_fst_final(impl->fst, best_tok->state);

  intlist_destroy(&words);
  return true;
}

void pk_decoder_result_init(pk_decoder_result_t *self) {
  self->words = NULL;
  self->size = 0;
  self->weight = 0.0f;
}

void pk_decoder_result_destroy(pk_decoder_result_t *best_path) {
  best_path->size = 0;
  free(best_path->words);
  best_path->words = NULL;

  best_path->weight = 0;
}
