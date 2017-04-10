// Created at 2017-04-10

#include "hashtable.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>
#include <unordered_map>

#define PRIME_NUM 16777259

std::vector<std::pair<int, int>> GenerateTestCases(int n, int max_val) {
  std::vector<std::pair<int, int>> test_cases;
  int key = 801;
  int value = 702;
  for (int i = 0; i < n; ++i) {
    key = abs((PRIME_NUM * key) % max_val);
    value = abs((PRIME_NUM * value) % max_val);
    test_cases.emplace_back(key, value);
  }

  return test_cases;
}

void TestHashTable() {
  std::vector<std::pair<int, int>> test_cases = GenerateTestCases(32767, 32767);
  pk_hashtable_t hashtable;
  pk_hashtable_init(&hashtable, PK_HASHTABLE_DEFAULTCAP);

  std::unordered_map<int, int> hash_map;

  for (const std::pair<int, int> &test_case : test_cases) {
    pk_hashtable_upsert(&hashtable, test_case.first, test_case.second);
    hash_map[test_case.first] = test_case.second;
  }
  for (int i = 0; i < test_cases.size(); ++i) {
    int val = pk_hashtable_find(&hashtable, i, -1);
    if (val == -1) {
      assert(hash_map.find(i) == hash_map.end());
    } else {
      assert(hash_map[i] == val);
    }
  }
}

void Benchmark() {
  int N = 100000;
  int M = 10000;
  std::vector<std::pair<int, int>>
  test_cases = GenerateTestCases(100000, 100000);

  pk_hashtable_t hashtable;
  pk_hashtable_init(&hashtable, PK_HASHTABLE_DEFAULTCAP);

  std::unordered_map<int, int> hash_map;

  puts("upsert:");
  clock_t t = clock();
  for (int i = 0; i < M; ++i) {
    for (const std::pair<int, int> &test_case : test_cases) {
      pk_hashtable_upsert(&hashtable, test_case.first, test_case.second);
    }
    pk_hashtable_clear(&hashtable);
  }
  t = clock() - t;
  printf("  pk_hashtable_t: %lfms\n", ((float)t) / CLOCKS_PER_SEC  * 1000);

  t = clock();
  for (int i = 0; i < M; ++i) {
    for (const std::pair<int, int> &test_case : test_cases) {
      hash_map[test_case.first] = test_case.second;
    }
    hash_map.clear();
  }
  t = clock() - t;
  printf("  std::unordered_map: %lfms\n", ((float)t) / CLOCKS_PER_SEC  * 1000);

  puts("lookup:");
  int sum = 0;
  t = clock();
  for (int i = 0; i < M; ++i) {
    for (int i = 0; i < test_cases.size(); ++i) {
      int val = pk_hashtable_find(&hashtable, i, -1);
      sum += val;
    }
  }
  t = clock() - t;
  printf("  pk_hashtable_t: %lfms\n", ((float)t) / CLOCKS_PER_SEC  * 1000);

  t = clock();
  for (int i = 0; i < M; ++i) {
    for (int i = 0; i < test_cases.size(); ++i) {
      std::unordered_map<int, int>::iterator it = hash_map.find(i);
      if (it != hash_map.end()) {
        sum += it->second;
      }
    }
  }
  t = clock() - t;
  printf("  std::unordered_map: %lfms\n", ((float)t) / CLOCKS_PER_SEC  * 1000);
  assert(sum != 0);
}

int main(int argc, char **argv) {
  TestHashTable();
  if (argc == 2 && strcmp(argv[1], "benchmark") == 0) {
    Benchmark();
  }
  return 0;
}
