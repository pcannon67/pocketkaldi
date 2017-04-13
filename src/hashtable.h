// hashlist.c --- Created at 2016-11-08
// hashtable.c ---  Created at 2017-04-10

#ifndef POCKETKALDI_HASHTABLE_H_
#define POCKETKALDI_HASHTABLE_H_

#include <stdint.h>
#include <stdbool.h>
#include "pocketkaldi.h"

#define PK_HASHTABLE_DEFAULTCAP 1024
#define PK_HASHTABLE_LOADFACTOR 0.5

typedef struct pk_hashtable_bucket_t {
  int32_t key;
  int32_t value;
} pk_hashtable_bucket_t;

// pk_hashtable_t is a struct to map int32_t to int32_t using hash table
typedef struct pk_hashtable_t {
  pk_hashtable_bucket_t *buckets;
  bool *is_empty;

  // bucket_size always equal to prime_numbers[prime_idx]
  int bucket_size;
  int prime_idx;

  int size;
} pk_hashtable_t;

// Initialize the hashtable. And the bucket size will be the nearest prime
// number of initial_cap. And initial_cap coule be PK_HASHTABLE_DEFAULTCAP
POCKETKALDI_EXPORT
void pk_hashtable_init(pk_hashtable_t *self, int initial_cap);

// Destroy the hashtable.
POCKETKALDI_EXPORT
void pk_hashtable_destroy(pk_hashtable_t *self);

// Insert or update an element in the hash table
POCKETKALDI_EXPORT
void pk_hashtable_upsert(pk_hashtable_t *self, int32_t key, int32_t value);

// Find element in the hash table. If the key doesn't exist, return default_val.
// Else return the value
POCKETKALDI_EXPORT
int32_t pk_hashtable_find(
    pk_hashtable_t *self,
    int32_t key,
    int32_t default_val);

// Clear all elements in hashtable
POCKETKALDI_EXPORT
void pk_hashtable_clear(pk_hashtable_t *self);

#endif  // POCKETKALDI_HASHTABLE_H_