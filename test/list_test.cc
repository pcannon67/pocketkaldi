// Created at 2016-11-22

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "alloc.h"
#include "list.h"

PKLIST_DEFINE(int, int_list)

void TestIntList() {
  int_list_t pklist_int;
  pk_alloc_t alloc;

  pk_alloc_init(&alloc);
  int_list_init(&pklist_int, &alloc);

  assert(int_list_empty(&pklist_int));

  std::vector<int> stdvector_int;
  int N = 1000000;
  for (int i = 0; i < N; ++i) {
    int random_num = rand();
    stdvector_int.push_back(random_num);
    int_list_push_back(&pklist_int, random_num);
  }

  assert(!int_list_empty(&pklist_int));

  for (int i = 0; i < N; ++i) {
    assert(stdvector_int[i] == pklist_int.data[i]);
  }
  assert(pklist_int.size == stdvector_int.size());
  assert(int_list_back(&pklist_int) == stdvector_int.back());

  // Test pop_back()
  int_list_pop_back(&pklist_int);
  stdvector_int.pop_back();
  assert(int_list_back(&pklist_int) == stdvector_int.back());
  assert(pklist_int.size == stdvector_int.size());

  int_list_destroy(&pklist_int);
}

int main() {
  TestIntList();
  return 0;
}
