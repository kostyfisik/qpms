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

void qpms_pr_debug_at_flf(const char *filename, unsigned int linenum,
		const char *func,
	       const char *fmt, ...);
//void qpms_error_at_line(const char *filename, unsigned int linenum,
//		const char *fmt, ...);

#define QPMS_CRASHING_MALLOC(pointer, size) {(pointer) = malloc(size); if(!pointer) qpms_pr_debug_at_flf(__FILE__,__LINE__,__func__, "Allocation of %zd bytes for " #pointer " failed.", (size_t) (size));}

#define QPMS_WTF qpms_pr_error_at_flf(__FILE__,__LINE__,__func__,"Unexpected error.")

#define QPMS_ENSURE_SUCCESS(x) {if(x) QPMS_WTF;}

#endif 
