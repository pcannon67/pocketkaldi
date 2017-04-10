// Created at 2017-04-10

#include "hashtable.h"

#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>

#define hash(x) x

static int prime_numbers[] = {
  2,
  2,
  5,
  11,
  17,
  37,
  67,
  131,
  257,
  521,
  1031,
  2053,
  4099,
  8209,
  16411,
  32771,
  65537,
  131101,
  262147,
  524309,
  1048583,
  2097169,
  4194319,
  8388617,
  16777259,
  33554467,
  67108879,
  134217757,
  268435459,
  536870923,
  1073741827,
  2147483647
};

// Get index of the nearest prime number of num in prime_numbers. And the
// prime number should greater or equal than num
static int nearest_prime_idx(int num) {
  assert(num <= 2147483647);  // num < 2^31 - 1

  int num_primes = sizeof(prime_numbers) / sizeof(int);
  for (int prime_idx = 0; prime_idx < num_primes; ++prime_idx) {
    int prime_number = prime_numbers[prime_idx];
    if (prime_number >= num) return prime_idx;
  }

  assert(false);
  return -1;
}

// Initialize buckets and is_empty array
static void init_buckets(
    pk_hashtable_bucket_t *buckets,
    bool *is_empty,
    int bucket_size) {
  memset(buckets, 0, sizeof(pk_hashtable_bucket_t) * bucket_size);
  for (int i = 0; i < bucket_size; ++i) {
    is_empty[i] = true;
  }
}

// Create and initialize buckets and is_empty array
static void create_buckets(
    pk_hashtable_bucket_t **buckets,
    bool **is_empty,
    int bucket_size) {
  *buckets = (pk_hashtable_bucket_t *)malloc(
      sizeof(pk_hashtable_bucket_t) * bucket_size);
  *is_empty = (bool *)malloc(sizeof(bool) * bucket_size);
  init_buckets(*buckets, *is_empty, bucket_size);
}

// Find the position that the key will be put. If the key already exists in
// buckets, returns the position to it. If the key doesn't exists, returns
// the position to put.
static inline int lookup(
    pk_hashtable_bucket_t *buckets,
    bool *is_empty,
    int buckets_size,
    int32_t key) {
  int hashval = hash(key) % buckets_size;
  for (int i = 1; ; ++i) {
    assert(i < 65536);
    for (int sign = -1; sign <= 1; sign += 2) {
      int squared = i * i;
      int idx = hashval + (sign == -1 ? -squared : squared);

      // Boundary
      idx = (idx + buckets_size) % buckets_size;

      if (is_empty[idx] == true) return idx;
      if (buckets[idx].key == key) return idx;
    }
  }
}

// Extend the hashtable, and put the elements into the new bucket
static void extend_buckets(pk_hashtable_t *self) {
  int new_size = prime_numbers[self->prime_idx + 1];

  // Allocate and initialize the new buckets
  pk_hashtable_bucket_t *new_buckets = NULL;
  bool *is_empty = NULL;
  create_buckets(&new_buckets, &is_empty, new_size);

  // Put each element into new buckets
  for (int i = 0; i < self->bucket_size; ++i) {
    if (self->is_empty[i]) continue;
    int new_idx = lookup(new_buckets, is_empty, new_size, self->buckets[i].key);
    assert(is_empty[new_idx] && "extend_buckets: element already exists");
    new_buckets[new_idx].key = self->buckets[i].key;
    new_buckets[new_idx].value = self->buckets[i].value;
    is_empty[new_idx] = false;
  }

  // Free the original buckets and is_empty array
  free(self->buckets);
  free(self->is_empty);

  // Replace the buckets and is_empty array
  self->buckets = new_buckets;
  self->is_empty = is_empty;
  self->bucket_size = new_size;
  self->prime_idx = self->prime_idx + 1;
}

void pk_hashtable_init(pk_hashtable_t *self, int initial_cap) {
  self->prime_idx = nearest_prime_idx(initial_cap);
  self->bucket_size = prime_numbers[self->prime_idx];

  // Initialize self->buckets and self->is_empty
  create_buckets(&self->buckets, &self->is_empty, self->bucket_size);

  self->size = 0;
}

// Destroy the hashtable.
void pk_hashtable_destroy(pk_hashtable_t *self) {
  self->size = 0;
  self->prime_idx = 0;
  self->bucket_size = 0;

  free(self->is_empty);
  self->is_empty = NULL;

  free(self->buckets);
  self->buckets = NULL;
}

void pk_hashtable_upsert(pk_hashtable_t *self, int32_t key, int32_t value) {
  double size = self->size;
  if (size / self->bucket_size > PK_HASHTABLE_LOADFACTOR) {
    extend_buckets(self);
  }

  int idx = lookup(
      self->buckets,
      self->is_empty,
      self->bucket_size,
      key);
  if (self->is_empty[idx]) {
    // Insert
    self->buckets[idx].key = key;
    self->buckets[idx].value = value;
    self->is_empty[idx] = false;
  } else {
    // Update
    assert(self->buckets[idx].key == key);
    self->buckets[idx].value = value;
  }
}

// Find element in the hash table. If the key doesn't exist, return default_val.
// Else return the value
int32_t pk_hashtable_find(
    pk_hashtable_t *self,
    int32_t key,
    int32_t default_val) {
  int idx = lookup(
      self->buckets,
      self->is_empty,
      self->bucket_size,
      key);
  if (self->is_empty[idx]) {
    return default_val;
  } else {
    assert(self->buckets[idx].key == key);
    return self->buckets[idx].value;
  }
}

// Clear all elements in hashtable
void pk_hashtable_clear(pk_hashtable_t *self) {
  self->size = 0;
  init_buckets(self->buckets, self->is_empty, self->bucket_size);
}
