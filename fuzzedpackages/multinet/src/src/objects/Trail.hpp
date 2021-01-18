#ifndef UU_OBJECTS_TRAIL_
#define UU_OBJECTS_TRAIL_

#include <unordered_set>
#include "objects/Walk.hpp"

namespace uu {
namespace net {

/**
 * A trail is a walk without repeated edges.
 */
class
    Trail :
    public Walk
{
  private:
    typedef Walk super;

  public:

    /**
     * Creates a path with a starting vertex.
     *
     * A trail is by definition not empty, therefore a starting vertex is needed.
     */
    Trail(
        const Vertex* v0
    );

    /**
     * Adds an edge (and the corresponding vertex) at the end of the trail.
     * @return the new end-vertex of the path.
     * @throw WrongParameterException if the edge does not start from the last vertex in the trail
     * @throw WrongParameterException if the edge is already part of the trail (except
     * the first)
     */
    const Vertex*
    extend(
        const Edge* v0
    );

    /**
     * Returns the number of edges in the trail.
     */
    size_t
    length(
    ) const;

    /**
     * Returns the vertices in the trail.
     */
    const std::list<const Vertex*>&
    vertices(
    ) const;

    /**
     * Returns the edges in the trail.
     */
    const std::list<const Edge*>&
    edges(
    ) const;

  private:

    /** Vertices */
    std::unordered_set<const Edge*> edge_set_;

};

}
}

#endif
