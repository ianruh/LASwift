#ifndef CLIB_SWIFT_CBLAS
#define CLIB_SWIFT_CBLAS

#if __APPLE__
// This may be redundant because the same header is imported by CLAPACK,
// but just incase CBLAS is used without CLAPACK, we keep it.
#include <Accelerate/Accelerate.h>;
#elif __linux__
#include <cblas.h>
#else
#error "Unknown compiler"
#endif

#endif
