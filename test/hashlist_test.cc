// Created at 2016-11-09

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unordered_map>
#include "hashlist.h"

void SimpleTest() {
  pk_hashlist_t hashlist;
  pk_alloc_t alloc;
  std::unordered_map<int, int> hashmap;
  char *null = NULL;

  pk_alloc_init(&alloc);
  pk_hashlist_init(&hashlist, &alloc);
  for (int i = 0; i < 100000; ++i) {
    int random_key = rand() % 3000;
    int random_value = rand() % 200000;

    pk_hashlist_insert(&hashlist, random_key, null + random_value);
    hashmap[random_key] = random_value;
  }

  printf(
      "hashlist.elem_size = %d, hashmap.size() = %d\n",
      hashlist.elem_size,
      (int)hashmap.size());
  assert(hashlist.elem_size == hashmap.size());
  pk_hashlist_elem_t *elem = hashlist.list_head;
  while (elem) {
    int key = elem->key;
    int value = (char *)elem->value - null;
    assert(hashmap[key] == value);
    elem = elem->next_list_elem;
  }

  pk_hashlist_destroy(&hashlist);
}

void Benchmark() {
  pk_hashlist_t hashlist;
  pk_alloc_t alloc;
  std::unordered_map<int, int> hashmap;
  char *null = NULL;
  srand(1024);

  int n_samples = 10000000;

  pk_alloc_init(&alloc);
  pk_hashlist_init(&hashlist, &alloc);
  int sum = 0;
  clock_t start = clock();
  for (int i = 0; i < n_samples; ++i) {
    int random_key = rand() % 3000;
    int random_value = rand() % 200000;

    pk_hashlist_insert(&hashlist, random_key, null + random_value);
  }
  pk_hashlist_elem_t *elem = hashlist.list_head;
  while (elem) {
    int key = elem->key;
    int value = (char *)elem->value - null;
    sum += value;
    elem = elem->next_list_elem;
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

  if (argc == 2 && strcmp(argv[1], "benchmark") == 0) {
    Benchmark();
  }
  
  return 0;
}

