#ifndef UU_COMMUNITY_INFOMAP_H_
#define UU_COMMUNITY_INFOMAP_H_

#ifndef NS_INFOMAP
#define NS_INFOMAP
#endif

//#include "infomap.hpp"
#include <cstddef> // to prevent '::max_align_t' has not been declared error
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <string>
#include "io/Config.h"
#include "infomap/InfomapContext.h"
#include "core/exceptions/ExternalLibException.hpp"
#include "community/_impl/infomap_utils.hpp"
#include "community/VertexLayerCommunity.hpp"
#include "community/CommunityStructure.hpp"

namespace uu {
namespace net {

/*  */
template <typename M>
std::unique_ptr<CommunityStructure<VertexLayerCommunity<const typename M::layer_type>>>
infomap(const M* net,
        bool overlapping=false,
        bool directed=false,
        bool include_self_links=true
       );


}
}

#include "infomap.ipp"

#endif
