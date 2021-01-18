#ifndef UU_NETWORKS_PROBABILISTICNETWORK_H_
#define UU_NETWORKS_PROBABILISTICNETWORK_H_

#include "networks/Network.hpp"

namespace uu {
namespace net {

/**
 * A ProbabilisticNetwork is a Network with a default edge attribute to store probabilities
 * and methods to set/get values on this attribute.
 */
class ProbabilisticNetwork
    : public Network
{

  private:

    std::string kPROB_ATTR_NAME = "p";

    typedef Network super;

  public:

    /**
     * Creates a ProbabilisticNetwork with directed or undirected simple edges and
     * with or without loops.
     */
    ProbabilisticNetwork(
        const std::string& name,
        EdgeDir dir = EdgeDir::UNDIRECTED,
        bool allow_loops = false
    );


    /**
     * Checks if the network is probabilistic.
     * Always returns true.
     */
    bool
    is_probabilistic(
    ) const override;


    /**
     * Sets the probability of an edge.
     * @throw uu::core::WrongParameterException if p < 0 or p > 1.
     */
    void
    set_prob(
        const Edge* e,
        double p
    );


    /**
     * Gets the probability of an edge.
     */
    core::Value<double>
    get_prob(
        const Edge* e
    ) const;

};

}
}

#endif
