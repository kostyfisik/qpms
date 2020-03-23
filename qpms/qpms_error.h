/*! \file qpms_error.h
 *
 * \brief QPMS miscellanous internal error handling functions and macros.
 *
 */
#ifndef QPMS_ERROR_H
#define QPMS_ERROR_H

#include "optim.h"

/// Print error message and abort();
QPMS_NORETURN void qpms_pr_error(const char *fmt, ...);
//void qpms_error(const char *fmt, ...);

/// Print an error message, indicating source, function name and line number, and abort().
/** 
 * Usually not used directly, but rather via some of the macros
 * that fill the first arguments automatically.
 *
 * \see QPMS_PR_ERROR
 */
QPMS_NORETURN void qpms_pr_error_at_flf(const char *filename, unsigned int linenum,
		const char *func,
	       const char *fmt, ...);

/// Print a warning message to stderr and flush the buffer. Don't call this directly, use QPMS_WARN().
void qpms_warn_at_flf(const char *filename, unsigned int linenum,
		const char *func,
	       const char *fmt, ...);

/// Classification of debugging messages.
/**
 * \see qpms_dbgmsg_enabled, qpms_dbgmsg_enable(), qpms_dbgmsg_disable()
 */
typedef enum {
	QPMS_DBGMSG_MISC = 1, 
	QPMS_DBGMSG_THREADS = 2, ///< Multithreading-related debug messages.
	QPMS_DBGMSG_INTEGRATION = 4 ///< Quadrature-related debug messages.
} qpms_dbgmsg_flags;

/// Print a debugging message to stderr and flush the buffer. Don't call this directly, use QPMS_DEBUG().
void qpms_debug_at_flf(const char *filename, unsigned int linenum,
		const char *func,
		qpms_dbgmsg_flags type,
	        const char *fmt, ...);

/// Global variable determining which types of debug messages shall be printed with QPMS_DEBUG().
/**
 * Use qpms_dbgmsg_enable() and qpms_dbgmsg_disable() to manipulate
 * this variable.
 *
 * \see QPMS_DEBUG()
 */
extern qpms_dbgmsg_flags qpms_dbgmsg_enabled;

/// Enable debugging messages of given \a types.
/** \see qpms_dbgmsg_disable() */
qpms_dbgmsg_flags qpms_dbgmsg_disable(qpms_dbgmsg_flags types);

/// Disable debugging messages of given \a types.
/** \see qpms_dbgmsg_enable() */
qpms_dbgmsg_flags qpms_dbgmsg_enable(qpms_dbgmsg_flags types);

/// Print a warning to stderr and flush the buffer.
/**
 * The arguments are the same as in standard printf().
 */
#define QPMS_WARN(msg, ...) qpms_warn_at_flf(__FILE__,__LINE__,__func__,msg, ##__VA_ARGS__)

/// Print a debugging message to stderr and flush the buffer.
/**
 * The arguments after \a type are the same as in standard printf().
 *
 * The debugging message is printed only if the corresponding \a type flag
 * is set in qpms_dbgmsg_enabled.
 *
 * \see qpms_dbgmsg_enabled
 */
#define QPMS_DEBUG(type /**< Debugging message type flag, see qpms_dbgmsg_flags */ , msg, ...) qpms_debug_at_flf(__FILE__,__LINE__,__func__,type,msg, ##__VA_ARGS__)

/// Wrapper macro of standard malloc(), crashing on failure.
/** 
 * Normally corresponds to a `pointer = malloc(size)`
 * statement;  however, if NULL is returned, this prints an 
 * error message and abort()s the program.
 *
 * The arguments are expanded several times.
 *
 * Note that this macro expands to a code block, to be kept in mind when using
 * together with if/else etc.
 *
 * Assigned memory block is to be deallocated with standard free().
 *
 * \see QPMS_CRASHING_CALLOC.
 */
#define QPMS_CRASHING_MALLOC(pointer, size) {\
	(pointer) = malloc(size);\
       	if(QPMS_UNLIKELY(!(pointer) && (size)))\
       		qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,\
			"Allocation of %zd bytes for " #pointer " failed.",\
			(size_t) (size));\
}

/// Allocate and copy.
/**
 * Behaves as QPMS_CRASHING_MALLOC(dest, size), but 
 * additionaly copies a chunk of memory from \a src to \a dest.
 *
 * \see QPMS_CRASHING_MALLOC()
 */
#define QPMS_CRASHING_MALLOCPY(dest, src, size) {\
	(dest) = malloc(size);\
       	if(QPMS_UNLIKELY(!(dest) && (size)))\
       		qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,\
			"Allocation of %zd bytes for " #dest " failed.",\
			(size_t) (size));\
	memcpy((dest), (src), (size));\
}

/// Wrapper macro of standard calloc(), crashing on failure.
/** 
 * Normally corresponds to a `pointer = calloc(nmemb, size)`
 * statement;  however, if NULL is returned, this prints an 
 * error message and abort()s the program.
 *
 * The arguments are expanded several times.
 *
 * Note that this macro expands to a code block, to be kept in mind when using
 * together with if/else etc.
 *
 * Assigned memory block is to be deallocated with standard free().
 *
 * \see QPMS_CRASHING_MALLOC, QPMS_CRASHING_REALLOC.
 */
#define QPMS_CRASHING_CALLOC(pointer, nmemb, size) {\
	(pointer) = calloc((nmemb), (size));\
       	if(QPMS_UNLIKELY(!(pointer) && (nmemb) && (size)))\
       		qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,\
			"Allocation of %zd bytes for " #pointer " failed.",\
			(size_t)((nmemb)*(size)));\
}

/// Wrapper macro of standard realloc(), crashing on failure.
/** 
 * Normally corresponds to a `pointer = realloc(pointer, size)`
 * statement;  however, if NULL is returned, this prints an 
 * error message and abort()s the program.
 *
 * The arguments are expanded several times.
 *
 * Note that this macro expands to a code block, to be kept in mind when using
 * together with if/else etc.
 *
 * Assigned memory block is to be deallocated with standard free().
 *
 * \see QPMS_CRASHING_MALLOC.
 */
#define QPMS_CRASHING_REALLOC(pointer, size) {\
	(pointer) = realloc((pointer), size);\
       	if(QPMS_UNLIKELY(!(pointer) && (size)))\
       		qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,\
			"Rellocation of %zd bytes for " #pointer " failed.",\
			(size_t) (size));\
}

/// Prints an "unexpected error" message and aborts the program.
/** Usually only put to presumably unreachable places in the code and similar. */
#define QPMS_WTF qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,"Unexpected error.")

/// Aborts the program with "invalid enumerator" error message.
#define QPMS_INVALID_ENUM(x) qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,"Invalid enum value (" #x " == %d)", (int) (x))

/// Prints an "untested function/feature" warning once when reached in the code.
#define QPMS_UNTESTED {\
	static _Bool already_bitched = 0; \
	if (QPMS_UNLIKELY(!already_bitched)) {\
		qpms_warn_at_flf(__FILE__,__LINE__,__func__,"Warning: untested function/feature!");\
		already_bitched = 1;\
	}\
}

/// Prints a given error message and aborts the program.
/** The arguments are as in standard printf(). */
#define QPMS_PR_ERROR(msg, ...) qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,msg, ##__VA_ARGS__)

/// Raises an error if \a x is not zero.
#define QPMS_ENSURE_SUCCESS(x) { \
	int errorcode = (x); /* evaluate x only once */ \
	if(QPMS_UNLIKELY(errorcode)) \
		qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,"Unexpected error (%d)", errorcode); \
}

/// Raises an error if \a x is not zero, with custom error message.
#define QPMS_ENSURE_SUCCESS_M(x, msg, ...) { \
	int errorcode = (x); \
	if(QPMS_UNLIKELY(errorcode)) \
		qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,msg, ##__VA_ARGS__); \
}

/// Raises an error if \a x is not zero or one of the values listed in arguments.
#define QPMS_ENSURE_SUCCESS_OR(x, ...) { \
	int errorcode = (x); /* evaluate x only once */ \
	static const int allowed_errorcodes[] = {0, ##__VA_ARGS__};\
	static const int n_allowed = sizeof(allowed_errorcodes) / sizeof(int); \
	int i; \
	for(i = 0; i < n_allowed; ++i) \
		if (errorcode == allowed_errorcodes[i]) break; \
	if (QPMS_UNLIKELY(i >= n_allowed))  \
		qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,"Unexpected error (%d)", errorcode); \
}

/// Raises an error if \a x is not true, with custom error message.
#define QPMS_ENSURE(x, msg, ...) {if(QPMS_UNLIKELY(!(x))) qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,msg, ##__VA_ARGS__); }

/// Raises an error if \a x is false.
/** 
 * Currently, this is always expanded, ignoring the possible NDEBUG macro.
 * In places where the evaluation could have significant performance impact,
 * consider using QPMS_PARANOID_ASSERT() instead.
 */
#define QPMS_ASSERT(x) {\
	if(QPMS_UNLIKELY(!(x)))\
       		qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,\
			"Unexpected error. This is most certainly a bug.");\
}

#ifdef QPMS_EVALUATE_PARANOID_ASSERTS
  /** \brief Raises an error if \a x is false.
   *
   * Expanded only if QPMS_EVALUATE_PARANOID_ASSERTS macro is defined.
   */
  #define QPMS_PARANOID_ASSERT(x) {\
	if(QPMS_UNLIKELY(!(x)))\
       		qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,\
			"Unexpected error. This is most certainly a bug.");\
}
#else
  #define QPMS_PARANOID_ASSERT(x) {;}
#endif

/// Raises a "not implemented" error with additional custom message.
/** Serves also as a label/placeholder of not implemented parts of the code. */
#define QPMS_NOT_IMPLEMENTED(msg, ...) qpms_pr_error_at_flf(__FILE__,__LINE__,__func__, \
		"Not implemented:" msg, ##__VA_ARGS__)

/// Prints an "incomplete implementation" warning once with a custom message.
/** Serves mainly as a label/placeholder of incomplete parts of the code. */
#define QPMS_INCOMPLETE_IMPLEMENTATION(msg, ...) {\
	static _Bool already_bitched = 0; \
	if (QPMS_UNLIKELY(!already_bitched)) {\
		qpms_warn_at_flf(__FILE__,__LINE__,__func__,msg, ##__VA_ARGS__);\
		already_bitched = 1;\
	}\
}

#endif 
