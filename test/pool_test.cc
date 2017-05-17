// pool_test.cc --- Created at 2016-05-17

#include "pool.h"

#include <assert.h>

struct A {
  A(int a) : a_(a) {}
  int a() const { return a_; }
  int a_;
};

struct B {
  B(int b) {
    b_ = new int(b);
  }
  ~B() {
    delete b_;
  }
  int a() const { return *b_; }
  int *b_;
};

template<typename T>
int TestType() {
  pocketkaldi::Pool<T, 256> pool;
  for (int i = 0; i < 1000; ++i) {
    T *p = pool.Alloc(i);
    assert(p->a() == i);
  }
  pool.Clear();
  std::vector<T *> pointer_list;
  for (int i = 1000; i < 4000; ++i) {
    T *p = pool.Alloc(i);
    pointer_list.push_back(p);
  }
  for (typename std::vector<T *>::iterator
       it = pointer_list.begin(); it != pointer_list.end(); ++it) {
    assert((*it)->a() == it - pointer_list.begin() + 1000);
  }
  return 0;
}

int main() {
  TestType<A>();
  TestType<B>();

  return 0;
}
