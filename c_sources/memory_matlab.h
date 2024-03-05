/* Memory managment for MATLAB */
#include "mex.h"


#ifdef __cplusplus
extern "C" {
#endif

  void* c_calloc(size_t num, size_t size);
  void* c_malloc(size_t size);
  void* c_realloc(void *ptr, size_t size);

#ifdef __cplusplus
}
#endif

#define c_free mxFree