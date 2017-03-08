// Created at 2016-11-08


#include "pocketkaldi.h"
#include "hashlist.h"

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
  pk_hashlist_elem_t *elem = hashlist->list_head;
  while (elem) {
    insert_to_bucket(hashlist, elem);
    elem = elem->next_list_elem;
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
  hashlist->elem_size = 0;
  hashlist->list_head = NULL;
}

void pk_hashlist_destroy(pk_hashlist_t *hashlist) {
  pk_free(hashlist->alloc, hashlist->buckets);
  hashlist->buckets = NULL;
  hashlist->bucket_size = 0;
  hashlist->elem_size = 0;
  hashlist->alloc = NULL;

  // Free all elements
  pk_hashlist_elem_t *elem = hashlist->list_head;
  while (elem) {
    pk_hashlist_elem_t *next = elem->next_list_elem;
    pk_free(hashlist->alloc, elem);
    elem = next;
  }
  hashlist->list_head = NULL;
}

void pk_hashlist_insert(
    pk_hashlist_t *hashlist,
    pk_hashlist_key_t key,
    pk_hashlist_value_t value) {
  // If there is no enough spaces, extend buckets
  double float_size = (double)hashlist->elem_size;
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
    pk_hashlist_elem_t *elem = (pk_hashlist_elem_t *)pk_alloc(
      hashlist->alloc,
      sizeof(pk_hashlist_elem_t));
    elem->key = key;
    elem->value = value;
    elem->next_list_elem = hashlist->list_head;
    hashlist->list_head = elem;

    // Insert into buckets
    insert_to_bucket(hashlist, elem);
    ++hashlist->elem_size;
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
  }

  return NULL; 
}

