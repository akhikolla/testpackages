#include "double_comp.h"

#include <math.h> /* fabs, fmin */

// Comparison function used by quicksort implementation in C
int doublecomp(const void *elem1, const void *elem2) {
    if (*(const double *)elem1 < *(const double *)elem2)
        return -1;
    return *(const double *)elem1 > *(const double *)elem2;
}

// Equality test with relative error precision
bool doubleeq(double a, double b, double reltol) {
    return (fabs(a-b) <= reltol*fmin(fabs(a),fabs(b)));
}
