// Create at 2017-03-27

#include "symbol_table.h"

#include <assert.h>
#include <stdlib.h>
#include "util.h"

void pk_symboltable_init(pk_symboltable_t *self) {
  self->size = 0;
  self->buffer = NULL;
  self->buffer_index = NULL;
}

void pk_symboltable_destroy(pk_symboltable_t *self) {
  self->size = 0;
  free(self->buffer);
  self->buffer = NULL;
  free(self->buffer_index);
  self->buffer_index = NULL;
}

void pk_symboltable_read(
    pk_symboltable_t *self,
    pk_readable_t *fd,
    pk_status_t *status) {
  int buffer_size;
  int expected_section_size;

  pk_bytebuffer_t bytebuffer;
  pk_bytebuffer_init(&bytebuffer);

  int section_size = pk_readable_readsectionhead(
      fd,
      PK_SYMBOLTABLE_SECTION,
      status);
  if (!status->ok) goto pk_symboltable_read_failed;
  pk_bytebuffer_reset(&bytebuffer, section_size);

  // Read size and check section size
  pk_readable_readbuffer(fd, &bytebuffer, status);
  if (!status->ok) goto pk_symboltable_read_failed;
  self->size = pk_bytebuffer_readint32(&bytebuffer);
  buffer_size = pk_bytebuffer_readint32(&bytebuffer);
  expected_section_size = 8 + self->size * sizeof(int) + buffer_size;
  if (section_size != expected_section_size) {
    PK_STATUS_CORRUPTED(
        status,
        "pk_symboltable_read: section_size = %d expected, but %d found (%s)",
        expected_section_size,
        section_size,
        fd->filename);
    goto pk_symboltable_read_failed;
  }

  // Read index
  self->buffer_index = (int *)malloc(self->size * sizeof(int));
  for (int i = 0; i < self->size; ++i) {
    int symbol_idx = pk_bytebuffer_readint32(&bytebuffer);
    self->buffer_index[i] = symbol_idx;
  }

  // Read symbol buffer
  self->buffer = (char *)malloc(buffer_size * sizeof(char));
  pk_bytebuffer_readbytes(&bytebuffer, self->buffer, buffer_size);

  if (false) {
pk_symboltable_read_failed:
    pk_symboltable_destroy(self);
  }

  pk_bytebuffer_destroy(&bytebuffer);
}

const char *pk_symboltable_get(const pk_symboltable_t *self, int symbol_id) {
  assert(symbol_id < self->size && "symbol_id out of boundary");
  int idx = self->buffer_index[symbol_id];
  return &self->buffer[idx];
}
