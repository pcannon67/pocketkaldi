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

#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>

/* swap size bytes between a_ and b_ */
static inline void swap(void *a_, void *b_, size_t size) {
  if (a_ == b_) return;
  {
    size_t i, nlong = size / sizeof(long);
    long *a = (long *) a_, *b = (long *) b_;
    for (i = 0; i < nlong; ++i) {
      long c = a[i];
      a[i] = b[i];
      b[i] = c;
    }
    a_ = (void*)(a + nlong);
    b_ = (void*)(b + nlong);
  }
  {
    size_t i;
    char *a = (char *) a_, *b = (char *) b_;
    size = size % sizeof(long);
    for (i = 0; i < size; ++i) {
      char c = a[i];
      a[i] = b[i];
      b[i] = c;
    }
  }
}

static inline void insertion_sort(
    void *base_,
    size_t nmemb,
    size_t size,
    int (*compar)(const void *, const void *, void *),
    void *thunk) {
  char *base = (char *) base_;
  size_t i, j;
  for (i = 0; i < nmemb; ++i) {
    for (j = i ; j-- > 0;) {
      if (compar(base + j * size, base + (j + 1) * size, thunk) <= 0) {
        break;
      }
      swap(base + j * size, base + (j + 1) * size, size);
    }
  }
}

static size_t log2_helper(size_t s) {
  size_t i = 0;
  while(s) {
    s >>= 1;
    i++;
  }
  return i;
}

static void sift_down_r(
    void *base_,
    size_t start,
    size_t nmemb,
    size_t size,
    int (*compar)(const void *, const void *, void *),
    void *thunk) {
  char *base = (char *) base_;
  while (start < nmemb / 2) {
    size_t child = start * 2 + 1;
    // choose larger child
    if (child < nmemb - 1 &&
        compar(base + child * size, base + (child + 1) * size, thunk) < 0) {
      child++;
    }
    if (compar(base + child * size, base + start * size, thunk) < 0) {
      return;
    }
    swap(base + start * size, base + child * size, size);
    start = child;
  }
}

static void make_heap_r(
    void *base_,
    size_t nmemb,
    size_t size,
    int (*compar)(const void *, const void *, void *),
    void *thunk) {
  char *base = (char *) base_;
  for (size_t start = nmemb / 2; start-- > 0 ; ) {
    sift_down_r(base_, start, nmemb, size, compar, thunk);
  }
}

static void pop_and_push_r(
    void *base_, size_t nmemb,
    size_t size, size_t nth,
    int (*compar)(const void *, const void *, void *),
    void *thunk) {
  char *base = (char *) base_;
  swap(base_, base+ nth * size, size);
  sift_down_r(base_, 0, nmemb, size, compar, thunk);
}

void pk_introselect_r(
    void *base_ptr,
    size_t nmemb,
    size_t size,
    size_t nth,
    int (*compar)(const void *, const void *, void *),
    void *thunk) {
  char *base = (char *) base_ptr;
  if (nmemb <= 1) {
    return;
  }

  size_t depth = log2_helper(nmemb * 2);
  while( nmemb > 3) {
    if (depth == 0) {
      make_heap_r(base, nth + 1, size, compar, thunk);
      size_t i;
      for (i = nth + 1; i < nmemb; i++) {
        if (compar(base + i * size, base, thunk) < 0) {
          pop_and_push_r(base, nth + 1, size, i, compar, thunk);
        }
      }
      swap(base, base + nth * size, size);
      return;
    }

    size_t i, pivot, npart;
    /* pick median of first/middle/last elements as pivot */
    {
      const char *a = base, *b = base + (nmemb / 2) * size,
            *c = base + (nmemb - 1) * size;
      pivot = compar(a, b, thunk) < 0
          ? (compar(b, c, thunk) < 0 ? nmemb / 2 :
             (compar(a, c, thunk) < 0 ? nmemb - 1 : 0))
            : (compar(a, c, thunk) < 0 ? 0 :
               (compar(b, c, thunk) < 0 ? nmemb - 1 : nmemb / 2));
    }

    /* partition array */
    swap(base + pivot * size, base + (nmemb - 1) * size, size);
    pivot = (nmemb - 1) * size;
    for (i = npart = 0; i < nmemb - 1; ++i)
      if (compar(base + i * size, base + pivot, thunk) <= 0)
        swap(base + i * size, base + (npart++)*size, size);
    swap(base + npart * size, base + pivot, size);

    if (nth == npart) {
      return;
    } else if (nth < npart) {
      nmemb = npart;
    } else { // nth > npart
      npart += 1;
      base = base + npart * size;
      nmemb -= npart;
      nth -= npart;
    }

    depth--;
  }

  insertion_sort(base, nmemb, size, compar, thunk);
}
