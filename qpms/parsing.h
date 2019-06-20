#ifndef QPMS_PARSING_H
#define QPMS_PARSING_H

#include <stddef.h>

/** Parse a given number of doubles from a string.
 *
 * The doubles can be separated by whitespaces, comma or semicolon.
 *
 * \return If the string included up to n doubles, number of parsed doubles.
 * If more, n+1.
 */
size_t qpms_parse_ndoubles(
  double *target,
  size_t n,
  const char *orig
);

/** Parse doubles from a string.
 *
 * The doubles can be separated by whitespaces, comma or semicolon.
 * The parsed numbers are saved into an array specified by *target
 * that has been preallocated with malloc() to contain at least start_index
 * members. If start_index is nonzero, the newly parsed numbers are 
 * saved to the positions starting from start_index.
 *
 * If *target is NULL, the function allocates the necessary space.
 *
 * \return Number of newly parsed doubles + start_index.
 */
size_t qpms_parse_doubles(
  double **target,
  size_t start_index,
  const char *orig
);

/** Parse doubles from a file.
 *
 * The doubles can be separated by whitespaces, comma or semicolon.
 * The parsed numbers are saved into an array specified by *target
 * that has been preallocated with malloc() to contain at least start_index
 * members. If start_index is nonzero, the newly parsed numbers are 
 * saved to the positions starting from start_index.
 *
 * If *target is NULL, the function allocates the necessary space.
 *
 * If filepath is NULL, "" or "-", read from stdin.
 *
 * \return Number of newly parsed doubles + start_index.
 */
size_t qpms_parse_doubles_fromfile(
  double **target,
  size_t start_index, //< Starting index for writing the parsed values.
  const char *filepath //< File to read from, or NULL, "", "-" to read from stdin.
);

#endif // QPMS_PARSING_H
