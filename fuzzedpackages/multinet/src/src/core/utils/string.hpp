/**
 * This header defines:
 * basic mathematical/statistical functions (mean, standard deviation, intersection, union).
 *
 * History:
 * - 2018.01.01 file imported from version 1.0 of the multinet library
 */

#ifndef UU_CORE_UTILS_STRING_H_
#define UU_CORE_UTILS_STRING_H_

#include <string>

namespace uu {
namespace core {


void
to_upper_case(
    std::string& s
);

void
format(
    std::string& in
);

}
}

#endif
