// Created at 2016-11-24

#include <assert.h>
#include <stdio.h>
#include <limits>
#include "fst.h"
#include "util.h"

using pocketkaldi::Fst;
using pocketkaldi::util::ReadableFile;
using pocketkaldi::Status;

void TestFst() {
  Status status;
  ReadableFile fd;
  status = fd.Open(TESTDIR "data/testinput.fst");
  assert(status.ok());


  Fst fst;
  status = fst.Read(&fd);
  assert(status.ok());

  assert(fst.start_state() == 0);
  assert(fst.final(0) == std::numeric_limits<double>::infinity());
  assert(fst.final(1) == std::numeric_limits<double>::infinity());
  assert(fst.final(2) == 3.5f);

  
  const Fst::Arc *arc = nullptr;

  Fst::ArcIterator arc_iter = fst.IterateArcs(0);
  arc = arc_iter.Next();
  assert(arc);
  assert(arc->next_state == 1);
  assert(arc->input_label == 1);
  assert(arc->output_label == 1);
  assert(arc->weight == 0.5f);
  arc = arc_iter.Next();
  assert(arc);
  assert(arc->next_state == 1);
  assert(arc->input_label == 2);
  assert(arc->output_label == 2);
  assert(arc->weight == 1.5f);
  arc = arc_iter.Next();
  assert(arc == nullptr);

  arc_iter = fst.IterateArcs(1);
  arc = arc_iter.Next();
  assert(arc);
  assert(arc->next_state == 2);
  assert(arc->input_label == 3);
  assert(arc->output_label == 3);
  assert(arc->weight == 2.5f);
  arc = arc_iter.Next();
  assert(arc == nullptr);

  arc_iter = fst.IterateArcs(2);
  arc = arc_iter.Next();
  assert(arc == nullptr);
}

int main() {
  TestFst();
  return 0;
}
