// Created at 2017-03-09

#include <assert.h>
#include "transition.h"

void pk_transition_init(
    pk_transition_t *self,
    pk_readable_t *fd,
    pk_status_t *status) {
  int section_size = pk_readable_readsectionhead(
      fd,
      PK_TRANSITION_SECTION,
      status);
  pk_bytebuffer_t bytebuffer;
  pk_bytebuffer_init(&bytebuffer);
  pk_bytebuffer_reset(&bytebuffer, section_size);
  if (status->ok) {
    if (section_size < 2 * sizeof(int32_t)) {
      PK_STATUS_CORRUPTED(status, "%s", fd->filename);
    }
  }
  if (status->ok) pk_readable_readbuffer(fd, &bytebuffer, status);

  int32_t size = 0;
  int32_t num_transtion = 0;
  int32_t num_pdf = 0;
  if (status->ok) {
    num_pdf = pk_bytebuffer_readint32(&bytebuffer);
    num_transtion = pk_bytebuffer_readint32(&bytebuffer);
    int expected_sectionsize = (num_transtion + 3) * sizeof(int32_t);
    if (section_size != expected_sectionsize) {
      PK_STATUS_CORRUPTED(status, "%s", fd->filename);
    }
  }
  if (status->ok) {
    self->num_transition_ids = num_transtion;
    self->num_pdfs = num_pdf;

    // Index of transition_id start from 1
    self->id2pdf = (int32_t *)malloc(sizeof(int32_t) * num_transtion + 1);
    for (int i = 0; i <= num_transtion; ++i) {
      self->id2pdf[i] = pk_bytebuffer_readint32(&bytebuffer);
    }
  }
  pk_bytebuffer_destroy(&bytebuffer);
}

void pk_transition_destroy(pk_transition_t *self) {
  self->num_transition_ids = 0;
  self->num_pdfs = 0;
  free(self->id2pdf);
  self->id2pdf = NULL;
}

int pk_transition_tid2pdf(const pk_transition_t *self, int transtion_id) {
  assert(transtion_id <= self->num_transition_ids);
  return self->id2pdf[transtion_id];
}
