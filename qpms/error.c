#include "qpms_error.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

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
//void qpms_error_at_line(const char *filename, unsigned int linenum,
//		const char *fmt, ...);

