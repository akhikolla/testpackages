#ifndef UU_COMMUNITY_ABACUS_H_
#define UU_COMMUNITY_ABACUS_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdio>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "core/exceptions/FileNotFoundException.hpp"
#include "core/exceptions/ExternalLibException.hpp"
#include "community/CommunityStructure.hpp"
#include "community/label_propagation.hpp"
#include "community/PillarCommunity.hpp"
#include "community/_impl/abacus_utils.hpp"




namespace uu {
namespace net {


/**
 * Finds communities using the abacus algorithm, itself using a label propagation
 * single-layer clustering algorithm as in the original paper.
 * @param mnet input multilayer network
 * @param min_actors minimum number of actors in a community
 * @param min_layers minimum number of layers in a community
 * @return a set of actor communities
 */
template <typename M, typename L>
std::unique_ptr<CommunityStructure<PillarCommunity<L>>>
abacus(
    const M* mnet,
    int min_actors,
    int min_layers
);

}
}

#include "abacus.ipp"

#endif
