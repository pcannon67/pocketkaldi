// Created at 2016-11-22

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "alloc.h"
#include "vector.h"

PKVECTOR_DEFINE(int, int_vector)

void TestIntVector() {
  int_vector_t pkvector_int;
  pk_alloc_t alloc;

  pk_alloc_init(&alloc);
  int_vector_init(&pkvector_int, &alloc);

  assert(int_vector_empty(&pkvector_int));

  std::vector<int> stdvector_int;
  int N = 1000000;
  for (int i = 0; i < N; ++i) {
    int random_num = rand();
    stdvector_int.push_back(random_num);
    int_vector_push_back(&pkvector_int, random_num);
  }

  assert(!int_vector_empty(&pkvector_int));

  for (int i = 0; i < N; ++i) {
    assert(stdvector_int[i] == pkvector_int.data[i]);
  }
  assert(pkvector_int.size == stdvector_int.size());
  assert(int_vector_back(&pkvector_int) == stdvector_int.back());

  // Test pop_back()
  int_vector_pop_back(&pkvector_int);
  stdvector_int.pop_back();
  assert(int_vector_back(&pkvector_int) == stdvector_int.back());
  assert(pkvector_int.size == stdvector_int.size());

  int_vector_destroy(&pkvector_int);
}

int main() {
  TestIntVector();
  return 0;
}
