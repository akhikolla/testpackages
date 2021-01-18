#ifndef UU_OPERATIONS_ANONYMIZE_H_
#define UU_OPERATIONS_ANONYMIZE_H_

#include <vector>
#include "core/exceptions/assert_not_null.hpp"
#include "core/utils/math.hpp"
#include "objects/Vertex.hpp"
#include "objects/EdgeMode.hpp"
#include "measures/degree.hpp"

namespace uu {
namespace net {


/**
 * @todo TO BE UPDATED
 * Replace actor names with random values
 * @param mnet A multilayer network.
 */
MLNetworkSharedPtr
anonymize_actors(const MLNetworkSharedPtr& mnet, const string& name);

}
}

#endif
