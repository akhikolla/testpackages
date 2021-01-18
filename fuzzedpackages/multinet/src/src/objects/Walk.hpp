#ifndef UU_OBJECTS_WALK_
#define UU_OBJECTS_WALK_

#include <list>
#include <iostream>
#include "objects/Edge.hpp"
#include "objects/Vertex.hpp"

namespace uu {
namespace net {

/**
 * A walk is a non-empty list (v0, e1, v1, ..., en, vn) where each edge has his preceding and
 * following vertices as end vertices. For directed edges, edge directionality must allow to
 * traverse the walk from left to right.
 */
class Walk
{

  public:

    /**
     * Creates a walk with a starting vertex.
     *
     * A walk is by definition not empty, therefore a starting vertex is needed.
     */
    Walk(
        const Vertex* v0
    );

    /**
     * Adds an edge (and the corresponding vertex) at the end of the walk.
     * @return the new end-vertex of the walk.
     * @throw WrongParameterException if the edge does not start from the last vertex in the walk
     */
    const Vertex*
    extend(
        const Edge* v0
    );

    /**
     * Returns the number of edges in the walk.
     */
    size_t
    length(
    ) const;

    /**
     * Returns the vertices in the walk.
     */
    const std::list<const Vertex*>&
    vertices(
    ) const;

    /**
     * Returns the edges in the walk.
     */
    const std::list<const Edge*>&
    edges(
    ) const;

    /** Output function, presenting a complete description of the walk. */
    std::string
    to_string(
    ) const;

  protected:

    /** Vertices */
    std::list<const Vertex*> vertices_;

    /** Edges */
    std::list<const Edge*> edges_;

};


std::ostream&
operator<<(std::ostream& os, const Walk& w);

}
}

#endif
