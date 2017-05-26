// Create at 2017-02-23

#include "decodable.h"

#include <math.h>
#include "am.h"

void pk_decodable_init(
    pk_decodable_t *self,
    AcousticModel *am,
    float prob_scale,
    const pk_matrix_t *feats) {
  pk_matrix_init(&self->log_prob, am->num_pdfs(), feats->ncol);
  am->Compute(feats, &self->log_prob);
  pk_matrix_scale(&self->log_prob, prob_scale);
  self->am = am;
}

void pk_decodable_destroy(pk_decodable_t *self) {
  pk_matrix_destroy(&self->log_prob);
  self->am = NULL;
}

float pk_decodable_loglikelihood(
    pk_decodable_t *self,
    int frame,
    int trans_id) {
  int pdf_id = self->am->TransitionIdToPdfId(trans_id);
  pk_vector_t frame_prob = pk_matrix_getcol(&self->log_prob, frame);
  return frame_prob.data[pdf_id];
}

bool pk_decodable_islastframe(pk_decodable_t *self, int frame) {
  assert(frame < self->log_prob.ncol);
  return frame == self->log_prob.ncol - 1;
}