#ifndef UU_OBJECTS_MULTIPLEXCLIQUE_H_
#define UU_OBJECTS_MULTIPLEXCLIQUE_H_

#include <utility>
#include <string>
#include <memory>
#include <unordered_set>

namespace uu {
namespace net {

/**
 * A multiplex clique.
 */
template <typename M>
class MultiplexClique
{
  public:

    /** Constructors */
    MultiplexClique(
        const std::unordered_set<const Vertex*>& actors,
        const std::unordered_set<const typename M::layer_type*>& layers
    );

    MultiplexClique();


    /** Comparison operator: equality, based on the presence of the same actors and the same layers. */
    bool
    operator==(
        const MultiplexClique&
    ) const;

    /** Comparison operator: difference, based on the presence of the same actors and the same layers. */
    bool
    operator!=(
        const MultiplexClique&
    ) const;

    /** Comparison operator: less then, based on the < operator of the actors and (if needed) the layers in the MultiplexClique. */
    bool
    operator<(
        const MultiplexClique&
    ) const;

    /** Comparison operator: greater than, based on the > operator of the actors and (if needed) the layers in the MultiplexClique. */
    bool
    operator>(
        const MultiplexClique&
    ) const;

    /** Output function, presenting a complete description of the MultiplexClique */
    std::string
    to_string(
    );

    /**
     * Creates an empty MultiplexClique and returns a pointer to it.
     * @return a pointer to the new multilayer network

    static std::shared_ptr<MultiplexClique<M>>
                                            create(
                                            );

    static std::shared_ptr<MultiplexClique<M>>
                                            create(
                                                    const std::unordered_set<const Vertex*>& actors,
                                                    const std::unordered_set<const typename M::layer_type*>& layers
                                            );
     */

    std::set<const Vertex*> actors;

    typename std::set<const typename M::layer_type*> layers;
};


}
}

#include "MultiplexClique.ipp"

#endif
