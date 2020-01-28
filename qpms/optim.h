/** \file optim.h
 * \brief Macros for compiler optimisation.
 */
#ifndef QPMS_OPTIM_H
#define QPMS_OPTIM_H


#if ((defined __GNUC__) || (defined __clang__)) && !(defined QPMS_NO_BUILTIN_EXPECT)
/// Wrapper over gcc's and clang's __builtin_expect.
/** If expands to __builtin_expect if gcc or clang are used,
 * else expands only to the first argument.
 */
#define QPMS_EXPECT(exp, c) __builtin_expect(exp, c)
#define QPMS_LIKELY(x) __builtin_expect(!!(x), 1)
#define QPMS_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define QPMS_LIKELY(x) (x)
#define QPMS_UNLIKELY(x)
#define QPMS_EXPECT(exp,c) (exp)
#endif

#if (defined(__GNUC__) && __GNUC__ >= 3) || \
    (defined(__GNUC__) && defined(__GNUC_MINOR__) && __GNUC__ == 2 && __GNUC_MINOR__ >= 8) 
// TODO clang
#define QPMS_NORETURN __attribute__((noreturn))
#endif

#endif // OPTIM_H
