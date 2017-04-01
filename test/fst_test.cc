// Created at 2016-11-24

#include <assert.h>
#include <stdio.h>
#include "fst.h"
#include "util.h"

void TestFst() {
  pk_fst_t fst;
  pk_status_t status;

  pk_status_init(&status);
  pk_readable_t *fd = pk_readable_open(TESTDIR "data/testinput.fst", &status);
  assert(status.ok);

  pk_fst_read(&fst, fd, &status);
  puts(TESTDIR "data/testinput.fst");
  assert(status.ok);

  assert(pk_fst_startstate(&fst) == 0);
  assert(pk_fst_final(&fst, 0) == 0);
  assert(pk_fst_final(&fst, 1) == 0);
  assert(pk_fst_final(&fst, 2) == 3.5f);

  pk_fst_iter_t it;
  const pk_fst_arc_t *arc;

  pk_fst_iterate_arc(&fst, 0, &it);
  arc = pk_fst_iter_next(&it);
  assert(arc);
  assert(arc->next_state == 1);
  assert(arc->input_label == 1);
  assert(arc->output_label == 1);
  assert(arc->weight == 0.5f);
  arc = pk_fst_iter_next(&it);
  assert(arc);
  assert(arc->next_state == 1);
  assert(arc->input_label == 2);
  assert(arc->output_label == 2);
  assert(arc->weight == 1.5f);
  arc = pk_fst_iter_next(&it);
  printf("%d\n", it.total);
  assert(arc == NULL);

  pk_fst_iterate_arc(&fst, 1, &it);
  arc = pk_fst_iter_next(&it);
  assert(arc);
  assert(arc->next_state == 2);
  assert(arc->input_label == 3);
  assert(arc->output_label == 3);
  assert(arc->weight == 2.5f);
  arc = pk_fst_iter_next(&it);
  assert(arc == NULL);

  pk_fst_iterate_arc(&fst, 2, &it);
  arc = pk_fst_iter_next(&it);
  assert(arc == NULL);

  pk_readable_close(fd);
}

int main() {
  TestFst();
  return 0;
}
