// hashlist.h --- Created at 2016-11-08
// hashtable.h ---  Created at 2017-04-10

#ifndef POCKETKALDI_HASHTABLE_H_
#define POCKETKALDI_HASHTABLE_H_

#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <functional>
#include "pocketkaldi.h"

namespace pocketkaldi {

extern int gPrimeNumbers[32];

// Hash fucntions
inline int32_t hash(int32_t x) {
  return x;
}

// A open addressing implementation of hashtable
// It assumes both K and V are prime types. So value copy won't cost too much
// time. Usually, K is int and V is int or pointer
template<typename K, typename V>
class HashTable {
 public:
  static constexpr int kDefaultCap = 1024;
  static constexpr double kLoadFactor = 0.5;

  // Initialize the hashtable. And the bucket size will be the nearest prime
  // number of initial_cap.
  HashTable(): HashTable(kDefaultCap) {}
  HashTable(int initial_cap) {
    prime_idx_ = NearestPrimeIdx(initial_cap);
    bucket_size_ = gPrimeNumbers[prime_idx_];

    CreateBuckets(&buckets_, &empty_, bucket_size_);
    size_ = 0;
  }

  // Destroy the hashtable.
  ~HashTable() {
    size_ = 0;
    prime_idx_ = 0;
    bucket_size_ = 0;

    delete[] empty_;
    empty_ = NULL;

    delete[] buckets_;
    buckets_ = NULL;
  }

  // Insert or update an element in the hash table
  void Insert(K key, V value) {
    double size = size_;
    if (size / bucket_size_ > kLoadFactor) {
      ExtendBuckets();
    }

    int idx = Lookup(buckets_, empty_, bucket_size_, key);
    if (empty_[idx]) {
      // Insert
      buckets_[idx].key = key;
      buckets_[idx].value = value;
      empty_[idx] = false;
      ++size_;
    } else {
      // Update
      assert(buckets_[idx].key == key);
      buckets_[idx].value = value;
    }
  }

  // Find element in the hash table. If the key doesn't exist, return
  // default_val, else return the value
  V Find(K key, V default_val) {
    int idx = Lookup(buckets_, empty_, bucket_size_, key);
    if (empty_[idx]) {
      return default_val;
    } else {
      assert(buckets_[idx].key == key);
      return buckets_[idx].value;
    }
  }

  // Clear all elements in hashtable
  void Clear() {
    size_ = 0;
    for (int i = 0; i < bucket_size_; ++i) {
      empty_[i] = true;
    }
  }

 private:
  struct Item {
    K key;
    V value;
  };

  Item *buckets_;
  bool *empty_;  // empty_[i] == true means buckets_[i] is empty

  // bucket_size always equal to prime_numbers[prime_idx]
  int bucket_size_;
  int prime_idx_;

  int size_;

  // Get index of the nearest prime number of num in gPrimeNumbers. And the
  // prime number should greater or equal than num
  int NearestPrimeIdx(int num) const {
    assert(num <= 2147483647);  // num < 2^31 - 1

    int num_primes = sizeof(gPrimeNumbers) / sizeof(int);
    for (int prime_idx = 0; prime_idx < num_primes; ++prime_idx) {
      int prime_number = gPrimeNumbers[prime_idx];
      if (prime_number >= num) return prime_idx;
    }

    assert(false);
    return -1;
  }

  // Create and initialize buckets and is_empty array
  void CreateBuckets(Item **buckets, bool **empty, int bucket_size) const {
    *buckets = new Item[bucket_size];
    *empty = new bool[bucket_size];
    for (int i = 0; i < bucket_size; ++i) {
      (*empty)[i] = true;
    }
  }

  // Find the position that the key will be put. If the key already exists in
  // buckets, returns the position to it. If the key doesn't exists, returns
  // the position to put.
  int Lookup(Item *buckets, bool *empty, int buckets_size, K key) const {
    int hashval = hash(key) % buckets_size;
    for (int i = 1; ; ++i) {
      assert(i < 65536);
      for (int sign = -1; sign <= 1; sign += 2) {
        int squared = i * i;
        int idx = hashval + (sign == -1 ? -squared : squared);

        // Boundary
        idx = idx % buckets_size;
        idx = idx < 0 ? idx + buckets_size : idx;

        if (empty[idx] == true) return idx;
        if (buckets[idx].key == key) return idx;
      }
    }
  }

  // Extend the hashtable, and put the elements into the new bucket
  void ExtendBuckets() {
    int new_size = gPrimeNumbers[prime_idx_ + 1];

    // Allocate and initialize the new buckets
    Item *new_buckets = NULL;
    bool *new_empty = NULL;
    CreateBuckets(&new_buckets, &new_empty, new_size);

    // Put each element into new buckets
    for (int i = 0; i < bucket_size_; ++i) {
      if (empty_[i]) continue;
      K key = buckets_[i].key;
      V value = buckets_[i].value;

      int new_idx = Lookup(new_buckets, new_empty, new_size, key);
      assert(new_empty[new_idx] && "ExtendBuckets: element already exists");

      new_buckets[new_idx].key = key;
      new_buckets[new_idx].value = value;
      new_empty[new_idx] = false;
    }

    // Free the original buckets and is_empty array
    delete[] buckets_;
    delete[] empty_;

    // Replace the buckets and is_empty array
    buckets_ = new_buckets;
    empty_ = new_empty;
    bucket_size_ = new_size;
    prime_idx_ = prime_idx_ + 1;
  }
};

}  // namespace pocketkaldi

#endif  // POCKETKALDI_HASHTABLE_H_
