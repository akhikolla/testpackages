#ifndef _ERR_DUMMY_H_
#define _ERR_DUMMY_H_

#include <stdlib.h>

#define warnx(...) do { \
        Rprintf ( __VA_ARGS__); \
        Rprintf ( "\n"); \
} while (0)

#define errx(code, ...) do { \
        Rprintf ( __VA_ARGS__); \
        Rprintf ( "\n"); \
        stop ("Error occured", code); \
} while (0)

#endif /* !_ERR_H_ */
