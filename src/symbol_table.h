// Create at 2017-03-27

#ifndef POCKETKALDI_SYMBOL_TABLE_H_
#define POCKETKALDI_SYMBOL_TABLE_H_

#include "util.h"

#define PK_SYMBOLTABLE_SECTION "SYM0"

// Store a list of symbols. And the symbol string could be got by
//   pk_symboltable_get(symboltable, symbol_id)
typedef struct pk_symboltable_t {
  int size;
  char *buffer;
  int *buffer_index;
} pk_symboltable_t;

// Initialize the symbol table
POCKETKALDI_EXPORT
void pk_symboltable_init(pk_symboltable_t *self);

// Destroy the symbol table
POCKETKALDI_EXPORT
void pk_symboltable_destroy(pk_symboltable_t *self);

// Read the symbol table from fd. If failed, status->ok will be set to false
POCKETKALDI_EXPORT
void pk_symboltable_read(
    pk_symboltable_t *self,
    pk_readable_t *fd,
    pk_status_t *status);

// Get the symbol by id, return the string of symbol.
POCKETKALDI_EXPORT
const char *pk_symboltable_get(const pk_symboltable_t *self, int symbol_id); 

#endif  // POCKETKALDI_SYMBOL_TABLE_H_
