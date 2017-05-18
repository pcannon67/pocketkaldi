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

  virtual ~Pool() {
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

  // Free and allocated nodes in this pool
  int free_nodes() const {
    return free_.size() +
           (blocks_.size() - current_block_) * BLOCK_SIZE -
           current_pos_;
  }
  int allocated_nodes() const {
    return current_block_ * BLOCK_SIZE + current_pos_ - free_.size();
  }

 protected:
  std::vector<T *> blocks_;
  std::vector<T *> free_;
  int current_block_;
  int current_pos_;
};

// Base class for any class T could be used in GCPool
class Collectable {
 public:
  enum {
    kMarked,
    kUnmarked,
    kFreed,
    kUsing
  };
  Collectable() : state_(kFreed) {}

  // Indicates this node is freed or under using
  void Free() { state_ = kFreed; }
  void Use() { state_ = kUsing; }

  // Mark and ummark in GC
  void Mark() { 
    assert(state_ != kFreed);
    state_ = kMarked; 
  }
  void Unmark() { 
    assert(state_ != kFreed);
    state_ = kUnmarked;
  }

  bool is_marked() const { return state_ == kMarked; }
  bool is_freed() const { return state_ == kFreed; }

 private:
  int8_t state_;
};

// GC pool add gabbage collection into Pool. Using mark-and-sweep algorithm. 
// And the class T linked with each othher with 'previous()' pointer.
// So the type T should have at least 4 members
//   T *previous()
//   void mark()
//   void unmark()
//   bool marked()
template<typename T, int BLOCK_SIZE = 4096>
class GCPool : public Pool<T, BLOCK_SIZE> {
 public:
  GCPool() {
  }

  template<typename... Args>
  T *Alloc(Args&&... args) {
    T *node = Pool<T, BLOCK_SIZE>::Alloc(std::forward<Args>(args)...);
    node->Use();
    return node;
  }

  void Dealloc(T *pointer) {
    pointer->Free();
    Pool<T, BLOCK_SIZE>::Dealloc(pointer);
  }

  // Start gabbage collection. Root set for GC is in root_nodes
  void GC(const std::vector<T *> root_nodes) {
    PK_DEBUG("=== GC Start ===");
    int free_nodes = this->free_nodes();
    int allocated_nodes = this->allocated_nodes();

    // Unmark all elements
    for (int block = 0; block <= this->current_block_; ++block) {
      int max_pos = block == this->current_block_ ?
          this->current_pos_ - 1 :
          BLOCK_SIZE;
      for (int pos = 0; pos < max_pos; ++pos) {
        // Bypass the nodes already freed
        if (this->blocks_[block][pos].is_freed()) continue;
        this->blocks_[block][pos].Unmark();
      }
    }

    // Push the root set into queue. Then traverse all nodes linked directly or
    // indirectly by root set
    std::vector<T *> queue(root_nodes.begin(), root_nodes.end());
    while (!queue.empty()) {
      T *node = queue.back();
      queue.pop_back();

      // Ignore nodes that already marked
      if (node->is_marked()) continue;
      
      // Freed node shouldn't be pointed by existing nodes
      assert(!node->is_freed());

      node->Mark();
      if (node->previous()) {
        queue.push_back(node->previous());
      }
    }

    // Free all unmarked nodes
    int free_count = 0;
    for (int block = 0; block <= this->current_block_; ++block) {
      int max_pos = block == this->current_block_ ?
          this->current_pos_ - 1 :
          BLOCK_SIZE;
      for (int pos = 0; pos < max_pos; ++pos) {
        // Bypass the nodes already freed
        T *node = this->blocks_[block] + pos;
        if (node->is_freed()) continue;
        if (node->is_marked()) continue;
        Dealloc(node);
        ++free_count;
      }
    }
    PK_DEBUG(util::Format(
        "Freed {} nodes, Allocated = {}",
        free_count,
        this->allocated_nodes()));
    PK_DEBUG("=== GC Finsihed ===");
  }
};

}  // namespace pocketkaldi

#endif  // POCKETKALDI_POOL_H_