// pool.h --- Created at 2016-05-17

#ifndef POCKETKALDI_POOL_H_
#define POCKETKALDI_POOL_H_

#include <stdio.h>
#include <vector>
#include <memory>
#include <type_traits>
#include "util.h"

namespace pocketkaldi {

// Pool is a class to optimize allocating large number of a class T.
template<typename T, int BLOCK_SIZE = 4096>
class Pool {
 public:
  Pool() : current_block_(0), current_pos_(0) {
    T *block = reinterpret_cast<T *>(::operator new(sizeof(T) * BLOCK_SIZE));
    blocks_.push_back(block);
  }

  ~Pool() {
    Clear();
    for (T *block : blocks_) {
      ::operator delete(block);
    }
    current_block_ = 0;
    current_pos_ = 0;
  }

  // Allocate a class T with constructor parameter args 
  template<typename... Args>
  T *Alloc(Args&&... args) {
    T *memory;
    if (free_.empty()) {
      assert(current_pos_ <= BLOCK_SIZE);
      if (current_pos_ == BLOCK_SIZE) {
        if (current_block_ == blocks_.size() - 1) {
          T *block = reinterpret_cast<T *>(
              ::operator new(sizeof(T) * BLOCK_SIZE));
          blocks_.push_back(block);
        }
        ++current_block_;
        current_pos_ = 0;
      }
      memory = blocks_[current_block_] + current_pos_;
      ++current_pos_;
    } else {
      memory = free_.back();
      free_.pop_back();
    }

    new (memory) T(std::forward<Args>(args)...);
    return memory;
  }

  // Allocate a class T with constructor parameter args 
  template<typename... Args>
  void Dealloc(T *pointer) {
    if (!std::is_trivially_destructible<T>::value) {
      pointer->~T();
    }
    free_.push_back(pointer);
  }

  // Clear all allocated class T
  void Clear() {
    // Call the dtor
    if (!std::is_trivially_destructible<T>::value) {
      for (int block_idx = 0; block_idx < blocks_.size() - 1; ++block_idx) {
        T *block = blocks_[block_idx];
        for (int idx = 0; idx < BLOCK_SIZE; ++idx) {
          block[idx].~T();
        }
      }

      // The last block
      T *last_block = blocks_.back();
      for (int idx = 0; idx < current_pos_; ++idx) {
        last_block[idx].~T();
      }
    }

    current_block_ = 0;
    current_pos_ = 0;
    free_.clear();
  }

 private:
  std::vector<T *> blocks_;
  std::vector<T *> free_;
  int current_block_;
  int current_pos_;
};

}  // namespace pocketkaldi

#endif  // POCKETKALDI_POOL_H_