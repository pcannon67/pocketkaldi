// Create at 2017-02-23

#include "decodable.h"

#include <math.h>
#include "am.h"

void pk_decodable_init(
    pk_decodable_t *self,
    pk_am_t *am,
    pk_transition_t *trans_model,
    float prob_scale,
    const pk_matrix_t *feats) {
  pk_matrix_init(&self->log_prob, trans_model->num_pdfs, feats->ncol);
  pk_am_compute(am, feats, &self->log_prob);
  pk_matrix_scale(&self->log_prob, prob_scale);

  self->trans_model = trans_model;
}

void pk_decodable_destroy(pk_decodable_t *self) {
  pk_matrix_destroy(&self->log_prob);
  self->trans_model = NULL;
}

float pk_decodable_loglikelihood(
    pk_decodable_t *self,
    int frame,
    int trans_id) {
  int pdf_id = pk_transition_tid2pdf(self->trans_model, trans_id);
  pk_vector_t frame_prob = pk_matrix_getcol(&self->log_prob, frame);
  return frame_prob.data[pdf_id];
}

bool pk_decodable_islastframe(pk_decodable_t *self, int frame) {
  assert(frame < self->log_prob.ncol);
  return frame == self->log_prob.ncol - 1;
}