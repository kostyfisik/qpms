#ifndef QPMS_ERROR_H
#define QPMS_ERROR_H

#include "optim.h"

/// Provisional error message with abort();
QPMS_NORETURN void qpms_pr_error(const char *fmt, ...);
//void qpms_error(const char *fmt, ...);

/// Provisional error message with abort(), indicating source and line number.
void qpms_pr_error_at_line(const char *filename, unsigned int linenum,
		const char *fmt, ...);

QPMS_NORETURN void qpms_pr_error_at_flf(const char *filename, unsigned int linenum,
		const char *func,
	       const char *fmt, ...);

/// Print a warning message to stderr and flush the buffer. Don't call this directly, use QPMS_WARN().
void qpms_warn_at_flf(const char *filename, unsigned int linenum,
		const char *func,
	       const char *fmt, ...);

void qpms_pr_debug_at_flf(const char *filename, unsigned int linenum,
		const char *func,
	       const char *fmt, ...);
//void qpms_error_at_line(const char *filename, unsigned int linenum,
//		const char *fmt, ...);


typedef enum {
	QPMS_DBGMSG_MISC = 1, 
	QPMS_DBGMSG_THREADS = 2, // Multithreading-related debug messages.
	QPMS_DBGMSG_INTEGRATION = 4 // Quadrature-related debug messages.
} qpms_dbgmsg_flags;

void qpms_debug_at_flf(const char *filename, unsigned int linenum,
		const char *func,
		qpms_dbgmsg_flags type,
	        const char *fmt, ...);

extern qpms_dbgmsg_flags qpms_dbgmsg_enabled;

qpms_dbgmsg_flags qpms_dbgmsg_disable(qpms_dbgmsg_flags types);
qpms_dbgmsg_flags qpms_dbgmsg_enable(qpms_dbgmsg_flags types);


#define QPMS_WARN(msg, ...) qpms_warn_at_flf(__FILE__,__LINE__,__func__,msg, ##__VA_ARGS__)

#define QPMS_DEBUG(type, msg, ...) qpms_debug_at_flf(__FILE__,__LINE__,__func__,type,msg, ##__VA_ARGS__)

#define QPMS_CRASHING_MALLOC(pointer, size) {\
	(pointer) = malloc(size);\
       	if(QPMS_UNLIKELY(!(pointer) && (size)))\
       		qpms_pr_debug_at_flf(__FILE__,__LINE__,__func__,\
			"Allocation of %zd bytes for " #pointer " failed.",\
			(size_t) (size));\
}

#define QPMS_CRASHING_MALLOCPY(dest, src, size) {\
	(dest) = malloc(size);\
       	if(QPMS_UNLIKELY(!(dest) && (size)))\
       		qpms_pr_debug_at_flf(__FILE__,__LINE__,__func__,\
			"Allocation of %zd bytes for " #dest " failed.",\
			(size_t) (size));\
	memcpy((dest), (src), (size));\
}

#define QPMS_CRASHING_CALLOC(pointer, nmemb, size) {\
	(pointer) = calloc((nmemb), (size));\
       	if(QPMS_UNLIKELY(!(pointer) && (nmemb) && (size)))\
       		qpms_pr_debug_at_flf(__FILE__,__LINE__,__func__,\
			"Allocation of %zd bytes for " #pointer " failed.",\
			(size_t)((nmemb)*(size)));\
}

#define QPMS_CRASHING_REALLOC(pointer, size) {\
	(pointer) = realloc((pointer), size);\
       	if(QPMS_UNLIKELY(!(pointer) && (size)))\
       		qpms_pr_debug_at_flf(__FILE__,__LINE__,__func__,\
			"Rellocation of %zd bytes for " #pointer " failed.",\
			(size_t) (size));\
}

#define QPMS_WTF qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,"Unexpected error.")

#define QPMS_INVALID_ENUM(x) qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,"Invalid enum value (" #x " == %d)", (int) (x))

#define QPMS_UNTESTED {\
	static _Bool already_bitched = 0; \
	if (QPMS_UNLIKELY(!already_bitched)) {\
		qpms_warn_at_flf(__FILE__,__LINE__,__func__,"Warning: untested function/feature!");\
		already_bitched = 1;\
	}\
}

#define QPMS_PR_ERROR(msg, ...) qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,msg, ##__VA_ARGS__)

/// Raises an error if \a x is not zero.
#define QPMS_ENSURE_SUCCESS(x) { \
	int errorcode = (x); /* evaluate x only once */ \
	if(QPMS_UNLIKELY(errorcode)) \
		qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,"Unexpected error (%d)", errorcode); \
}

// Same as previous, but with message.
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

#define QPMS_ENSURE(x, msg, ...) {if(QPMS_UNLIKELY(!(x))) qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,msg, ##__VA_ARGS__); }

#define QPMS_ASSERT(x) {\
	if(QPMS_UNLIKELY(!(x)))\
       		qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,\
			"Unexpected error. This is most certainly a bug.");\
}

#ifdef QPMS_EVALUATE_PARANOID_ASSERTS
  #define QPMS_PARANOID_ASSERT(x) {\
	if(QPMS_UNLIKELY(!(x)))\
       		qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,\
			"Unexpected error. This is most certainly a bug.");\
}
#else
  #define QPMS_PARANOID_ASSERT(x) {;}
#endif

#define QPMS_NOT_IMPLEMENTED(msg, ...) qpms_pr_error_at_flf(__FILE__,__LINE__,__func__, \
		"Not implemented:" msg, ##__VA_ARGS__)

#define QPMS_INCOMPLETE_IMPLEMENTATION(msg, ...) {\
	static _Bool already_bitched = 0; \
	if (QPMS_UNLIKELY(!already_bitched)) {\
		qpms_warn_at_flf(__FILE__,__LINE__,__func__,msg, ##__VA_ARGS__);\
		already_bitched = 1;\
	}\
}

#endif 
