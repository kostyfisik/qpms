#ifndef QPMS_ERROR_H
#define QPMS_ERROR_H

/// Provisional error message with abort();
void qpms_pr_error(const char *fmt, ...);
//void qpms_error(const char *fmt, ...);

/// Provisional error message with abort(), indicating source and line number.
void qpms_pr_error_at_line(const char *filename, unsigned int linenum,
		const char *fmt, ...);

void qpms_pr_error_at_flf(const char *filename, unsigned int linenum,
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


#define QPMS_WARN(msg, ...) qpms_warn_at_flf(__FILE__,__LINE__,__func__,msg, ##__VA_ARGS__)

#define QPMS_CRASHING_MALLOC(pointer, size) {(pointer) = malloc(size); if(!pointer && (size)) qpms_pr_debug_at_flf(__FILE__,__LINE__,__func__, "Allocation of %zd bytes for " #pointer " failed.", (size_t) (size));}

#define QPMS_CRASHING_REALLOC(pointer, size) {(pointer) = realloc(pointer, size); if(!pointer && (size)) qpms_pr_debug_at_flf(__FILE__,__LINE__,__func__, "Rellocation of %zd bytes for " #pointer " failed.", (size_t) (size));}

#define QPMS_WTF qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,"Unexpected error.")

#define QPMS_PR_ERROR(msg, ...) qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,msg, ##__VA_ARGS__)

#define QPMS_ENSURE_SUCCESS(x) {if(x) QPMS_WTF;}

#define QPMS_ENSURE(x, msg, ...) {if(!(x)) qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,msg, ##__VA_ARGS__); }

#define QPMS_NOT_IMPLEMENTED(msg, ...) qpms_pr_error_at_flf(__FILE__,__LINE__,__func__, \
		"Not implemented:" msg, ##__VA_ARGS__)

#endif 
