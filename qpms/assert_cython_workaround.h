// The purpose of this file is to enable assertions in cython modules.
// By default, cython includes -DNDEBUG argument when running gcc and 
// it seems this can not be disabled. Therefore, we force undefining
// NDEBUG in the code if DISABLE_NDEBUG is defined.
#ifndef ASSERT_CYTHON_H
#define ASSERT_CYTHON_H

#ifdef DISABLE_NDEBUG
#undef NDEBUG
#endif

#include <assert.h>

#endif // ASSERT_CYTHON_H
