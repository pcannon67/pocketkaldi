// Created at 2016-11-24

#ifndef POCKETKALDI_FST_H_
#define POCKETKALDI_FST_H_

#include <stdint.h>
#include "util.h"
#include "status.h"
#include "pocketkaldi.h"

namespace pocketkaldi {

class Fst {
public:
  // An out-going arc in Fst
  struct Arc {
    int32_t next_state;
    int32_t input_label;
    int32_t output_label;
    float weight;
  };

  // Iterators of out-going arcs for a state
  class ArcIterator {
   public:
    ArcIterator(int base, int total, const Arc *arcs);
    ~ArcIterator();

    // If next arc exists, retrun it and move the iterator forward, else return
    // nullptr
    const Arc *Next();

   private:
    int base_;
    int cnt_pos_;
    int total_;
    const Arc *arcs_;
  };

  // Section name in binary file
  static constexpr const char *kSectionName = "pk::fst_0";

  Fst();

  // Read fst from binary file.
  Status Read(util::ReadableFile *fd);

  // Start state of this Fst
  int start_state() const {
    return start_state_;
  }

  // Get the final score of state. If the state is non-terminal, returns 0
  float final(int state_id) const {
    assert(state_id < final_.size());
    return final_[state_id];
  }

  // Iterate out-going arcs for a state
  ArcIterator IterateArcs(int state) const;

 private:
  int start_state_;
  std::vector<Arc> arcs_;
  std::vector<int32_t> state_idx_;
  std::vector<float> final_;

  // Calcuate the number of outcoming arcs for state
  int CountArcs(int state) const;
};

}  // namespace pocketkaldi

#endif  // POCKETKALDI_FST_H_