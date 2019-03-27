#include "qpms_error.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

qpms_dbgmsg_flags qpms_dbgmsg_enabled = QPMS_DBGMSG_MISC;

qpms_dbgmsg_flags qpms_dbgmsg_enable(qpms_dbgmsg_flags types) {
  return (qpms_dbgmsg_enabled |= types);
}

qpms_dbgmsg_flags qpms_dbgmsg_disable(qpms_dbgmsg_flags types) {
  return (qpms_dbgmsg_enabled &= ~types);
}

void qpms_pr_error(const char *fmt, ...) {
  fflush(stdout);
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  fflush(stderr);
  abort();
}


//void qpms_error(const char *fmt, ...);

void qpms_pr_error_at_line(const char *filename, unsigned int linenum,
		const char *fmt, ...) {
  fflush(stdout);
  fprintf(stderr, "%s:%u: ", filename, linenum);
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  fflush(stderr);
  abort();
}

void qpms_pr_error_at_flf(const char *filename, unsigned int linenum,
    const char *func,
		const char *fmt, ...) {
  fflush(stdout);
  fprintf(stderr, "%s:%u, %s: ", filename, linenum, func);
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  fflush(stderr);
  abort();
}

void qpms_warn_at_flf(const char *filename, unsigned int linenum,
    const char *func,
		const char *fmt, ...) {
  fprintf(stderr, "%s:%u, %s: ", filename, linenum, func);
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  fflush(stderr);
}

void qpms_pr_debug_at_flf(const char *filename, unsigned int linenum,
    const char *func,
		const char *fmt, ...) {
  fprintf(stderr, "%s:%u, %s: ", filename, linenum, func);
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  fflush(stderr);
}

void qpms_debug_at_flf(const char *filename, unsigned int linenum,
    const char *func, qpms_dbgmsg_flags type,
		const char *fmt, ...) {
  if (!(type & qpms_dbgmsg_enabled)) return;
  fprintf(stderr, "%s:%u, %s: ", filename, linenum, func);
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  fflush(stderr);
}

//void qpms_error_at_line(const char *filename, unsigned int linenum,
//		const char *fmt, ...);

