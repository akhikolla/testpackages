#ifndef UU_NETWORKS_WEIGHTEDNETWORK_H_
#define UU_NETWORKS_WEIGHTEDNETWORK_H_

#include "networks/Network.hpp"

namespace uu {
namespace net {

/**
 * A weighted network is a Network with a default edge attribute to store weights
 * and methods to set/get values on this attribute.
 */
class WeightedNetwork
    : public Network
{

  private:

    std::string kWEIGHT_ATTR_NAME = "weight";

    typedef Network super;

  public:

    /**
     * Creates a WeightedNetwork with directed or undirected simple edges and with or without loops.
     */
    WeightedNetwork(
        const std::string& name,
        EdgeDir dir = EdgeDir::UNDIRECTED,
        bool allow_loops = false
    );


    /**
     * Checks if the network is weighted.
     * Always returns true.
     */
    bool
    is_weighted(
    ) const override;


    /**
     * Sets the weight of an edge.
     */
    void
    set_weight(
        const Edge* e,
        double w
    );


    /**
     * Gets the weight of an edge.
     */
    core::Value<double>
    get_weight(
        const Edge* e
    ) const;

};

}
}

#endif
