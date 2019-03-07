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

#endif 
