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

class C : public pocketkaldi::Collectable {
 public:
  C *referent;
  C(C *referent): referent(referent) {}
  C *previous() const {
    return referent;
  }
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

void TestGC() {
  pocketkaldi::GCPool<C> pool;
  C *a1 = pool.Alloc(nullptr);
  C *a2 = pool.Alloc(a1);
  C *b1 = pool.Alloc(nullptr);
  C *c1 = pool.Alloc(nullptr);
  C *c2 = pool.Alloc(c1);
  C *c3 = pool.Alloc(c2);
  C *d1 = pool.Alloc(c2);
  pool.Dealloc(a2);

  std::vector<C *> root_set{d1};
  pool.GC(root_set);

  assert(a1->is_freed());
  assert(a2->is_freed());
  assert(b1->is_freed());
  assert(c3->is_freed());

  assert(!c1->is_freed());
  assert(!c2->is_freed());
  assert(!d1->is_freed());

  pool.Dealloc(c2);
  assert(c2->is_freed());
}

int main() {
  TestType<A>();
  TestType<B>();
  TestGC();

  return 0;
}
