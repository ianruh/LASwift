#ifndef CLIB_SWIFT_CLAPACK
#define CLIB_SWIFT_CLAPACK

#if __APPLE__
    #include <Accelerate/Accelerate.h>;
#elif __linux__
    #include <lapacke.h>
#else
#   error "Unknown compiler"
#endif

#endif
