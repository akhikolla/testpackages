/**

 */

#ifndef UU_NET_DATASTRUCTURES_STORES_MULTIEDGESTORE_H_
#define UU_NET_DATASTRUCTURES_STORES_MULTIEDGESTORE_H_

#include <unordered_set>
#include "networks/_impl/stores/EdgeStore.hpp"

namespace uu {
namespace net {


class MultiEdgeStore :
    public EdgeStore
{

  private:
    typedef EdgeStore super;

  public:

    MultiEdgeStore(
        EdgeDir dir
    );


    /**
     * Adds a new edge.
     * Multiple edges between the same pair of vertices are allowed.
     * @param node1 a pointer to the "from" vertex if directed, or to one end of
     * the edge if undirected.
     * @param node2 a pointer to the "to" vertex if directed, or one end of the
     * edge if undirected.
     * @return a pointer to the new edge, or nullptr if the edge already exists.
     **/
    virtual
    const Edge *
    add(
        const Vertex* vertex1,
        const Vertex* vertex2
    );

    /**
     * Adds a new edge.
     * Multiple edges between the same pair of vertices are allowed.
     * @param node1 a pointer to the "from" vertex if directed, or to one end of
     * the edge if undirected.
     * @param node2 a pointer to the "to" vertex if directed, or one end of the
     * edge if undirected.
     * @return a pointer to the new edge, or nullptr if the edge already exists.
     **/
    virtual
    const Edge *
    add(
        std::shared_ptr<const Edge>  e
    ) override;


    /**
     * Returns all edges between two vertices.
     * This function can also be used to check if an edge is present.
     * @param node1 a pointer to the "from" vertex if directed, or to one end
     * of the edge if undirected.
     * @param node2 a pointer to the "to" vertex if directed, or one end
     * of the edge if undirected.
     * @return a pointer to the requested edge, or NULL if it does not exist.
     **/
    core::SortedRandomSet<const Edge*>
    get(
        const Vertex* vertex1,
        const Vertex* vertex2
    ) const;


    using EdgeStore::neighbors;
    using EdgeStore::is_directed;
    using EdgeStore::attach;
    using EdgeStore::summary;

    using EdgeStore::size;
    using EdgeStore::edge_directionality;
    using EdgeStore::sidx_neighbors_out;
    using EdgeStore::sidx_neighbors_in;
    using EdgeStore::sidx_neighbors_all;
    using EdgeStore::observers;

    /**
     * @brief Deletes an existing edge.
     * Attribute values associated to this edge are also deleted.
     * @param edge a pointer to the edge to be deleted
     * @return true if the object has been erased, false if it was not present.
     **/
    bool
    erase(
        const Edge * const e
    ) override;

    void
    erase(
        const Vertex* const v
    ) override;

    std::string
    summary(
    ) const override;

  protected:

    // Indexes to objects (Component IDX):
    std::unordered_map<const Vertex*, std::unordered_map<const Vertex*, std::unordered_set<const Edge*> > >
    cidx_edges_by_vertices;

};

}
}

#endif
