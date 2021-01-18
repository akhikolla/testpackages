#ifndef UU_NETWORKS_TEMPORALNETWORK_H_
#define UU_NETWORKS_TEMPORALNETWORK_H_

#include "networks/MultiNetwork.hpp"

namespace uu {
namespace net {

/**
 * A TemporalNetwork is a MultiNetwork with a default edge attribute to store times
 * and methods to set/get values on this attribute.
 */

class TemporalNetwork
    : public MultiNetwork
{

  private:

    typedef MultiNetwork super;

  public:

    const std::string kTIME_ATTR_NAME = "t";

    /**
     * Creates a TemporalNetwork with directed or undirected multiedges and with or without loops.
     */
    TemporalNetwork(
        const std::string& name,
        EdgeDir dir = EdgeDir::UNDIRECTED,
        bool allow_loops = true
    );


    /**
     * Checks if the network is temporal.
     * Always returns true.
     */
    bool
    is_temporal(
    ) const override;


    /**
     * Sets the time of an edge.
     */
    void
    set_time(
        const Edge* e,
        core::Time t
    );


    /**
     * Gets the time of an edge.
     */
    core::Value<core::Time>
    get_time(
        const Edge* e
    ) const;


    /**
     * Gets the smallest edge time in the network.
     */
    core::Value<core::Time>
    get_min_time(
    ) const;


    /**
     * Gets the highest edge time in the network.
     */
    core::Value<core::Time>
    get_max_time(
    ) const;


};

}
}

#endif
