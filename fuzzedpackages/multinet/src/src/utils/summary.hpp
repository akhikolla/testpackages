#ifndef UU_UTILS_SUMMARY_H_
#define UU_UTILS_SUMMARY_H_

#include <memory>
#include <string>

namespace uu {
namespace net {

/**
 * Returns a string summarizing the input network that can be printed on a single line.
 */
template<typename G>
std::string
summary_short(
    const G* g
);

}
}

#include "summary.ipp"

#endif
