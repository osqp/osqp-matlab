#include <mex.h>

void* c_calloc(size_t num, size_t size) {
  void *m = mxCalloc(num, size);
  mexMakeMemoryPersistent(m);
  return m;
}

void* c_malloc(size_t size) {
  void *m = mxMalloc(size);
  mexMakeMemoryPersistent(m);
  return m;
}

void* c_realloc(void *ptr, size_t size) {
  void *m = mxRealloc(ptr, size);
  mexMakeMemoryPersistent(m);
  return m;
}