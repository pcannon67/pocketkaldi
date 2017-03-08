// Created at 2016-11-08


#include "hashlist.h"

#include <assert.h>
#include "pocketkaldi.h"

#define PK_HASHLIST_INITIALSIZE 16
#define MAX_ELEMBUCKET_RATIO 1.0
#define INCREASE_RATE 1.5
#define hash(x) x

// Insert one element into bucket. Only changes buckets amd the elememt, other
// parts like link list or sizes will not be changed
static void insert_to_bucket(
    pk_hashlist_t *hashlist,
    pk_hashlist_elem_t *elem) {
  int hashval = hash(elem->key) % hashlist->bucket_size;
  pk_hashlist_elem_t **first_elem = &(hashlist->buckets[hashval]);

  // Put elem at the head of bucket
  elem->next_bucket_elem = *first_elem;
  *first_elem = elem;
}

// Extend bucket number to `new_size`. And put each element into its new
// position
static void extend_buckets(pk_hashlist_t *hashlist, int new_size) {
  pk_hashlist_elem_t **new_buckets = (pk_hashlist_elem_t **)pk_alloc(
      hashlist->alloc,
      sizeof(pk_hashlist_elem_t *) * new_size);
  for (int i = 0; i < new_size; ++i) {
    new_buckets[i] = NULL;
  }
  
  pk_free(hashlist->alloc, hashlist->buckets);
  hashlist->buckets = new_buckets;
  hashlist->bucket_size = new_size;

  // Put each element into new bucket
  pk_hashlist_elem_t *elem = hashlist->head;
  while (elem) {
    insert_to_bucket(hashlist, elem);
    elem = elem->next;
  }
}

void pk_hashlist_init(pk_hashlist_t *hashlist, pk_alloc_t *alloc) {
  hashlist->alloc = alloc;
  hashlist->buckets = (pk_hashlist_elem_t **)pk_alloc(
      hashlist->alloc,
      sizeof(pk_hashlist_elem_t *) * PK_HASHLIST_INITIALSIZE);
  for (int i = 0; i < PK_HASHLIST_INITIALSIZE; ++i) {
    hashlist->buckets[i] = NULL;
  }
  hashlist->bucket_size = PK_HASHLIST_INITIALSIZE;
  hashlist->size = 0;
  hashlist->head = NULL;
  hashlist->empty_head = NULL;
}

void pk_hashlist_destroy(pk_hashlist_t *hashlist) {
  pk_free(hashlist->alloc, hashlist->buckets);
  hashlist->buckets = NULL;
  hashlist->bucket_size = 0;
  hashlist->size = 0;
  hashlist->alloc = NULL;

  // Free all elements
  pk_hashlist_elem_t *elem = hashlist->head;
  while (elem) {
    pk_hashlist_elem_t *next = elem->next;
    pk_free(hashlist->alloc, elem);
    elem = next;
  }

  elem = hashlist->empty_head;
  while (elem) {
    pk_hashlist_elem_t *next = elem->next;
    pk_free(hashlist->alloc, elem);
    elem = next;
  }

  hashlist->head = NULL;
  hashlist->empty_head = NULL;
}

void pk_hashlist_insert(
    pk_hashlist_t *hashlist,
    pk_hashlist_key_t key,
    pk_hashlist_value_t value) {
  // If there is no enough spaces, extend buckets
  double float_size = (double)hashlist->size;
  if (float_size / hashlist->bucket_size > MAX_ELEMBUCKET_RATIO) {
    int new_size = (int)(hashlist->bucket_size * INCREASE_RATE);
    extend_buckets(hashlist, new_size);
  }

  // Check if it is already exists in hashlist
  int hashval = hash(key) % hashlist->bucket_size;
  pk_hashlist_elem_t *bucket = hashlist->buckets[hashval];
  pk_hashlist_elem_t *start = bucket;
  bool has_updated = false;
  while (start) {
    if (start->key == key) {
      start->value = value;
      has_updated = true;
      break;
    }
    start = start->next_bucket_elem;
  }

  // If key is not exist, allocate the element and put into hashlist
  if (!has_updated) {
    pk_hashlist_elem_t *elem = NULL;

    // If empty list have unused elements, just pick one from it
    if (hashlist->empty_head) {
      elem = hashlist->empty_head;
      hashlist->empty_head = elem->next;
    } else {
      elem = (pk_hashlist_elem_t *)pk_alloc(
          hashlist->alloc,
          sizeof(pk_hashlist_elem_t));
    }
    elem->key = key;
    elem->value = value;
    elem->next = hashlist->head;
    elem->next_bucket_elem = NULL;  // It will be changed soon
    hashlist->head = elem;

    // Insert into buckets
    insert_to_bucket(hashlist, elem);
    ++hashlist->size;
  }

}

pk_hashlist_elem_t *pk_hashlist_find(
    pk_hashlist_t *hashlist,
    pk_hashlist_key_t key) {
  int hashval = hash(key) % hashlist->bucket_size;
  pk_hashlist_elem_t *bucket = hashlist->buckets[hashval];
  pk_hashlist_elem_t *p = bucket;
  while (p) {
    if (p->key == key) return p;
    p = p->next_bucket_elem;
  }

  return NULL; 
}

void pk_hashlist_clear(pk_hashlist_t *hashlist) {
  for (int i = 0; i < hashlist->bucket_size; ++i) {
    hashlist->buckets[i] = NULL;
  }
  hashlist->size = 0;
  hashlist->empty_head = hashlist->head;
  hashlist->head = NULL;
}

void pk_hashlist_swap(pk_hashlist_t *hashlist1, pk_hashlist_t *hashlist2) {
  // We must assure they are using the same allocator
  assert(hashlist1->alloc == hashlist2->alloc);

  pk_hashlist_elem_t **buckets = hashlist1->buckets;
  hashlist1->buckets = hashlist2->buckets;
  hashlist2->buckets = buckets;

  int size = hashlist1->bucket_size;
  hashlist1->bucket_size = hashlist2->bucket_size;
  hashlist2->bucket_size = size;
  size = hashlist1->size;
  hashlist1->size = hashlist2->size;
  hashlist2->size = size;

  pk_hashlist_elem_t *head = hashlist1->head;
  hashlist1->head = hashlist2->head;
  hashlist2->head = head;
  head = hashlist1->empty_head;
  hashlist1->empty_head = hashlist2->empty_head;
  hashlist2->empty_head = head;
}
