// Created at 2016-11-22

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "alloc.h"
#include "vector.h"

PKVECTOR_DEFINE(int)

void TestIntVector() {
  pk_vector_int_t pkvector_int;
  pk_alloc_t alloc;

  pk_alloc_init(&alloc);
  pk_vector_int_init(&pkvector_int, &alloc);

  std::vector<int> stdvector_int;
  int N = 1000000;
  for (int i = 0; i < N; ++i) {
    int random_num = rand();
    stdvector_int.push_back(random_num);
    pk_vector_int_push_back(&pkvector_int, random_num);
  }

  for (int i = 0; i < N; ++i) {
    assert(stdvector_int[i] == pkvector_int.data[i]);
  }

  pk_vector_int_destroy(&pkvector_int);
}

int main() {
  TestIntVector();
  return 0;
}
