#ifndef UU_OBJECTS_PATH_
#define UU_OBJECTS_PATH_

#include <unordered_set>
#include "objects/Walk.hpp"

namespace uu {
namespace net {

/**
 * A path is a walk without repeated vertices (except the first and last, in which case
 * the path is a cycle).
 */
class
    Path :
    public Walk
{
  private:
    typedef Walk super;

  public:

    /**
     * Creates a path with a starting vertex.
     *
     * A path is by definition not empty, therefore a starting vertex is needed.
     */
    Path(
        const Vertex* v0
    );

    /**
     * Adds an edge (and the corresponding vertex) at the end of the path.
     * @return the new end-vertex of the path.
     * @throw WrongParameterException if the edge does not start from the last vertex in the walk
     * @throw WrongParameterException if the edge contains vertices already in the path (except
     * the first)
     */
    const Vertex*
    extend(
        const Edge* v0
    );

    /**
     * Returns the number of edges in the path.
     */
    size_t
    length(
    ) const;

    /**
     * Returns the vertices in the path.
     */
    const std::list<const Vertex*>&
    vertices(
    ) const;

    /**
     * Returns the edges in the path.
     */
    const std::list<const Edge*>&
    edges(
    ) const;

    /**
     * Returns true if the first and last vertices in the path are the same.
     */
    bool
    is_cycle(
    ) const;

  private:

    /** Vertices */
    std::unordered_set<const Vertex*> vertex_set_;

};

}
}

#endif
