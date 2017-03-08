// Created at 2016-11-08

#include <stdbool.h>
#include <stdint.h>
#include "pocketkaldi.h"

#ifndef POCKETKALDI_HASHLIST_H_
#define POCKETKALDI_HASHLIST_H_

typedef int32_t pk_hashlist_key_t;
typedef void * pk_hashlist_value_t;

struct pk_hashlist_elem_t {
  pk_hashlist_key_t key;
  pk_hashlist_value_t value;
  struct pk_hashlist_elem_t *next_bucket_elem;
  struct pk_hashlist_elem_t *next;
};
typedef struct pk_hashlist_elem_t pk_hashlist_elem_t;

struct pk_hashlist_t {
  pk_hashlist_elem_t **buckets;
  int bucket_size;
  int size;
  pk_hashlist_elem_t *head;
  pk_hashlist_elem_t *tail;  // Points to the last element of list
  pk_hashlist_elem_t *empty_head;
};
typedef struct pk_hashlist_t pk_hashlist_t;

// Initialize the hashlist
POCKETKALDI_EXPORT
void pk_hashlist_init(pk_hashlist_t *hashlist);

// Destory the hashlist, free all elements it allocates. But the memory that
// each value points to will not be freed
POCKETKALDI_EXPORT
void pk_hashlist_destroy(pk_hashlist_t *hashlist);

// Insert a key-value pair into hashlist. If the key already exists, just update
// the value corresponded to it.
POCKETKALDI_EXPORT
void pk_hashlist_insert(
    pk_hashlist_t *hashlist,
    pk_hashlist_key_t key,
    pk_hashlist_value_t value);

// Find the key-value pair (element) from the hashlist by its key. If the key
// not exist, returns NULL
POCKETKALDI_EXPORT
pk_hashlist_elem_t *pk_hashlist_find(
    pk_hashlist_t *hashlist,
    pk_hashlist_key_t key);

// Clear all elements in hashlist
POCKETKALDI_EXPORT
void pk_hashlist_clear(pk_hashlist_t *hashlist);

// Swap two hashlists. They should have the same allocator
POCKETKALDI_EXPORT
void pk_hashlist_swap(pk_hashlist_t *hashlist1, pk_hashlist_t *hashlist2);

#endif  // POCKETKALDI_HASHLIST_H_
