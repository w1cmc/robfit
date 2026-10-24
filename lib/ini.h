/* inih -- simple .INI file parser

SPDX-License-Identifier: BSD-3-Clause

Copyright (C) 2009-2025, Ben Hoyt

inih is released under the New BSD license (see LICENSE.txt). Go to the project
home page for more info:

https://github.com/benhoyt/inih

*/

#ifndef INI_H
#define INI_H

/* Make this header file easier to include in C++ code */
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

/* Nonzero if ini_handler callback should accept lineno parameter. */
#ifndef INI_HANDLER_LINENO
#define INI_HANDLER_LINENO 0
#endif

/* Visibility symbols, required for Windows DLLs */
#ifndef INI_API
#define INI_API
#endif

/* Typedef for prototype of handler function. */
typedef int (*ini_handler)(void* user, const char* section,
                           const char* name, const char* value);

/* Typedef for prototype of fgets-style reader function. */
typedef char* (*ini_reader)(char* str, int num, void* stream);

/* Parse given INI-style file. May have [section]s, name=value pairs
   (whitespace stripped), and comments starting with ';' (semicolon). Section
   is "" if name=value pair parsed before any section heading.

   For each name=value pair parsed, call handler function with given user
   pointer as well as section, name, and value (data only valid for duration
   of handler call). Handler should return nonzero on success, zero on error.

   Returns 0 on success, line number of first error on parse error (doesn't
   stop on first error), -1 on file open error, or -2 on memory allocation
   error (only when INI_USE_STACK is zero).
*/
INI_API int ini_parse(const char* filename, ini_handler handler, void* user);

/* Same as ini_parse(), but takes a FILE* instead of filename. */
INI_API int ini_parse_file(FILE* file, ini_handler handler, void* user);

/* Same as ini_parse(), but takes an ini_reader function pointer instead of
   filename. Used for implementing custom or string-based I/O. */
INI_API int ini_parse_stream(ini_reader reader, void* stream, ini_handler handler,
                             void* user);

/* Same as ini_parse(), but takes a zero-terminated string with the INI data
   instead of a file. */
INI_API int ini_parse_string(const char* string, ini_handler handler, void* user);

/* Configuration options */
#ifndef INI_ALLOW_MULTILINE
#define INI_ALLOW_MULTILINE 1
#endif

#ifndef INI_ALLOW_BOM
#define INI_ALLOW_BOM 1
#endif

#ifndef INI_START_COMMENT_PREFIXES
#define INI_START_COMMENT_PREFIXES ";#"
#endif

#ifndef INI_ALLOW_INLINE_COMMENTS
#define INI_ALLOW_INLINE_COMMENTS 1
#endif

#ifndef INI_INLINE_COMMENT_PREFIXES
#define INI_INLINE_COMMENT_PREFIXES ";"
#endif

#ifndef INI_USE_STACK
#define INI_USE_STACK 1
#endif

#ifndef INI_MAX_LINE
#define INI_MAX_LINE 512
#endif

#ifndef INI_MAX_SECTION
#define INI_MAX_SECTION 50
#endif

#ifndef INI_MAX_NAME
#define INI_MAX_NAME 50
#endif

#ifndef INI_STOP_ON_FIRST_ERROR
#define INI_STOP_ON_FIRST_ERROR 0
#endif

#ifndef INI_CALL_HANDLER_ON_NEW_SECTION
#define INI_CALL_HANDLER_ON_NEW_SECTION 0
#endif

#ifndef INI_ALLOW_NO_VALUE
#define INI_ALLOW_NO_VALUE 0
#endif

#ifdef __cplusplus
}
#endif

#endif /* INI_H */
