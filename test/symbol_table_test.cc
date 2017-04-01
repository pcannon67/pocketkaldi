// Created at 2017-04-01

#include "symbol_table.h"

#include <assert.h>
#include <string.h>
#include "pocketkaldi.h"
#include "util.h"

void TestSymbolTable() {
  pk_status_t status;
  pk_status_init(&status);

  pk_readable_t *fd = pk_readable_open(
      TESTDIR "data/symboltable_test.bin",
      &status);
  assert(status.ok);

  pk_symboltable_t symbol_table;
  pk_symboltable_init(&symbol_table);
  pk_symboltable_read(&symbol_table, fd, &status);
  assert(status.ok);

  assert(strcmp(pk_symboltable_get(&symbol_table, 0), "hello") == 0);
  assert(strcmp(pk_symboltable_get(&symbol_table, 1), "world") == 0);
  assert(strcmp(pk_symboltable_get(&symbol_table, 2), "cat") == 0);
  assert(strcmp(pk_symboltable_get(&symbol_table, 3), "milk") == 0);

  pk_symboltable_destroy(&symbol_table);
  pk_readable_close(fd);
}

int main() {
  TestSymbolTable();
  return 0;
}
