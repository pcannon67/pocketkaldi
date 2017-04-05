/* Copyright (c) 2017-2017 Weng Xuetian
 * Copyright (c) 2007-2014 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

// Created at 2017-04-04

#include "util.h"

#include <assert.h>
#include <stdlib.h>
#include <array>
#include <algorithm>

int intcmp(const void *a, const void *b, void *thunk) {
  return (*((int*)a)) - (*((int*)b));
}

#define ARRAY_SIZE(ARRAY) (sizeof(ARRAY) / sizeof(ARRAY[0]))

int main() {
  srand(time(NULL));

  for (int i = 0; i < 100; ++i) {
    std::array<int, 1024> qsort_array;
    std::array<int, 1024> array;
    for (int j = 0; j < qsort_array.size(); j++) {
      qsort_array[j] = array[j] = rand();
    }

    int idx = rand() % 1024;

    // Here we test the boundary
    switch (i) {
    case 0:
      idx = 0;
    case 1:
      idx = qsort_array.size() - 1;
    }

    pk_introselect_r(
        array.data(),
        array.size(),
        sizeof(int),
        idx,
        intcmp,
        NULL);
    std::sort(qsort_array.begin(), qsort_array.end());

    assert(qsort_array[idx] == array[idx]);
  }

  return 0;
}