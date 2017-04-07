// Created at 2016-11-09

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unordered_map>
#include "hashlist.h"

void TestHashList(pk_hashlist_t *hashlist) {
  std::unordered_map<int, int> hashmap;
  char *null = NULL;

  for (int i = 0; i < 3000; ++i) {
    int random_key = rand() % 3000;
    int random_value = rand() % 200000;

    pk_hashlist_insert(hashlist, random_key, null + random_value);
    hashmap[random_key] = random_value;
  }

  assert(hashlist->size == hashmap.size());
  pk_hashlist_elem_t *elem = hashlist->head;
  while (elem) {
    int key = elem->key;
    int value = (char *)elem->value - null;
    assert(hashmap[key] == value);
    elem = elem->next;
  }
} 

void SimpleTest() {
  pk_hashlist_t hashlist;
  pk_hashlist_init(&hashlist);

  TestHashList(&hashlist);

  pk_hashlist_destroy(&hashlist);
}

void TestClear() {
  pk_hashlist_t hashlist;
  pk_hashlist_init(&hashlist);

  TestHashList(&hashlist);
  pk_hashlist_clear(&hashlist);
  pk_hashlist_clear(&hashlist);
  TestHashList(&hashlist);
  pk_hashlist_clear(&hashlist);
  TestHashList(&hashlist);

  pk_hashlist_destroy(&hashlist);
}

void TestSwap() {
  char *NULL_CPTR = NULL;
  pk_hashlist_t hashlist1, hashlist2;
  pk_hashlist_init(&hashlist1);
  pk_hashlist_init(&hashlist2);

  pk_hashlist_insert(&hashlist1, 23, NULL_CPTR + 23);
  pk_hashlist_insert(&hashlist1, 29, NULL_CPTR + 29);
  pk_hashlist_insert(&hashlist1, 31, NULL_CPTR + 31);
  pk_hashlist_insert(&hashlist1, 37, NULL_CPTR + 37);

  pk_hashlist_insert(&hashlist2, 8, NULL_CPTR + 8);
  pk_hashlist_insert(&hashlist2, 16, NULL_CPTR + 16);

  pk_hashlist_swap(&hashlist1, &hashlist2);
  pk_hashlist_elem_t *elem = pk_hashlist_find(&hashlist1, 16);
  assert(elem && elem->value == NULL_CPTR + 16);
  elem = pk_hashlist_find(&hashlist1, 23);
  assert(elem == NULL);
  elem = pk_hashlist_find(&hashlist2, 8);
  assert(elem == NULL);
  elem = pk_hashlist_find(&hashlist2, 31);
  assert(elem && elem->value == NULL_CPTR + 31);

  pk_hashlist_destroy(&hashlist1);
  pk_hashlist_destroy(&hashlist2);
}

void Benchmark() {
  pk_hashlist_t hashlist;
  std::unordered_map<int, int> hashmap;
  char *null = NULL;
  srand(1024);

  int n_samples = 10000000;

  pk_hashlist_init(&hashlist);
  int sum = 0;
  clock_t start = clock();
  for (int i = 0; i < n_samples; ++i) {
    int random_key = rand() % 3000;
    int random_value = rand() % 200000;

    pk_hashlist_insert(&hashlist, random_key, null + random_value);
  }
  pk_hashlist_elem_t *elem = hashlist.head;
  while (elem) {
    int key = elem->key;
    int value = (char *)elem->value - null;
    sum += value;
    elem = elem->next;
  }
  clock_t end = clock();
  printf("pk_hashlist: sum = %d, time = %ld\n", sum, end - start);
  pk_hashlist_destroy(&hashlist);

  srand(1024);
  sum = 0;
  start = clock();
  for (int i = 0; i < n_samples; ++i) {
    int random_key = rand() % 3000;
    int random_value = rand() % 200000;

    hashmap[random_key] = random_value;
  }
  for (const std::pair<int, int> &item : hashmap) {
    sum += item.second;
  }
  end = clock();
  printf("std::unordered_map: sum = %d, time = %ld\n", sum, end - start);  
}

int main(int argc, char **argv) {
  SimpleTest();
  TestClear();

  if (argc == 2 && strcmp(argv[1], "benchmark") == 0) {
    Benchmark();
  }
  
  return 0;
}

