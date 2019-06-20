#include "parsing.h"
#include "qpms_error.h"
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

size_t qpms_parse_ndoubles(
  double *target,
  size_t n,
  const char *orig
) {
  QPMS_ENSURE(target, "The target parameter must not be NULL");
  char * const dup = strdup(orig);
  QPMS_ENSURE(dup, "Memory error in a strdup() call.");

  // Replace commas and semicolons with whitespaces
  for (char *c = dup; *c; ++c) 
    if (*c == ',' || *c == ';')
      *c = ' ';

  errno = 0;
  size_t i = 0;

  const char *beg = dup;
  while(*beg) {
    char *endptr;
    double parsed = strtod(beg, &endptr);
    if (endptr > beg) {
      if (i >= n) {
        errno = EOVERFLOW;
        if (i == n) QPMS_WARN("Supplied string contains additional numbers"
           " (expected only %zd numbers): %s\n", n, beg);
      }
      else
        target[i] = parsed;
      ++i;
      beg = endptr;
    } else {
      while (*beg) {
        if (!isspace(*beg)) {
          QPMS_WARN("Invalid character (expected a double), leaving the rest of the string unprocessed: %s\n", beg);
          errno = EILSEQ;
          goto qpms_parse_ndoubles_cleanup;
        }
        ++beg;
      }
    }
  }

qpms_parse_ndoubles_cleanup:
  free(dup);
  return i;
}


size_t qpms_parse_doubles(
  double **target,
  size_t start_index,
  const char *orig
) {
  QPMS_ENSURE(target, "The target parameter must not be NULL");
  char * const dup = strdup(orig);
  QPMS_ENSURE(dup, "Memory error in a strdup() call.");

  size_t capacity = start_index * 2;
  if (capacity < 128) capacity = 128;

  // Replace commas and semicolons with whitespaces
  for (char *c = dup; *c; ++c) 
    if (*c == ',' || *c == ';')
      *c = ' ';

  size_t i = start_index;
  errno = 0;

  const char *beg = dup;
  while(*beg) {
    char *endptr;
    errno = 0;
    double parsed = strtod(beg, &endptr);
    if (endptr > beg) {
      (*target)[i] = parsed;
      ++i;
      if (i >= capacity) {
        capacity *= 2;
        QPMS_CRASHING_REALLOC(*target, capacity * sizeof(double));
      }
      beg = endptr;
    } else {
      while (*beg) {
        if (!isspace(*beg)) {
          QPMS_WARN("Invalid character (expected a double), leaving the rest of the string unprocessed: %s\n", beg);
          errno = EILSEQ;
          goto qpms_parse_doubles_cleanup;
        }
        ++beg;
      }
    }
  }

qpms_parse_doubles_cleanup:
  free(dup);
  return i;
}


size_t qpms_parse_doubles_fromfile(
  double **target,
  size_t start_index, //< Starting index for writing the parsed values.
  const char *filepath //< File to read from, or NULL, "", "-" to read from stdin.
) {
  QPMS_ENSURE(target, "The target parameter must not be NULL");

  FILE *src;

  if (!filepath || !strcmp(filepath, "-") || filepath[0]=='\0')
    src = stdin;
  else 
    QPMS_ENSURE(src = fopen(filepath, "f"),
        "Could not open file %s: %s", filepath, strerror(errno));

  char buf[1024];
  int scanresult;
  while (1 == (scanresult = fscanf(src, "%1023s", buf)))
    start_index = qpms_parse_doubles(target, start_index, buf);

  if (errno) QPMS_WARN("Problem reading %s: %s", 
     (src==stdin) ? "stdin" : filepath, strerror(errno));

qpms_parse_doubles_files_cleanup:
  if (src != stdin)
    QPMS_ENSURE(!fclose(src),
        "Could not close file %s: %s", filepath, strerror(errno));

  return start_index;
}
