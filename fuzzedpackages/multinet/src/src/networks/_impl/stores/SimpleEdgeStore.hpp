/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_DATASTRUCTURES_STORES_SIMPLEEDGESTORE2_H_
#define UU_NET_DATASTRUCTURES_STORES_SIMPLEEDGESTORE2_H_

#include <unordered_map>
#include "networks/_impl/stores/EdgeStore.hpp"

namespace uu {
namespace net {

class SimpleEdgeStore :
    public EdgeStore
{

  private:
    typedef EdgeStore super;

  public:

    /**
     * Constructor.
     */
    SimpleEdgeStore(
        EdgeDir dir
    ) :
        super(dir)
    {}

    /**
     * Adds a new edge.
     * @param e edge to be added.
     * @return a pointer to the new edge, or nullptr if the edge already exists.
     **/
    virtual
    const Edge*
    add(
        std::shared_ptr<const Edge>  e
    ) override;

    /**
     * Returns an edge.
     * This function can also be used to check if an edge is present.
     * @param vertex1 a pointer to the "from" actor if directed, or to one end
     * of the edge if undirected.
     * @param vertex2 a pointer to the "to" actor if directed, or one end
     * of the edge if undirected.
     * @return a pointer to the requested edge, or nullptr if it does not exist.
     **/
    const Edge*
    get(
        const Vertex* vertex1,
        const Vertex* vertex2
    ) const;

    const Edge*
    add(
        const Vertex* v1,
        const Vertex* v2
    );

    using super::add;
    using super::neighbors;
    using super::is_directed;
    using super::attach;
    using super::size;
    using super::summary;

    using super::edge_directionality;
    using super::sidx_neighbors_out;
    using super::sidx_neighbors_in;
    using super::sidx_neighbors_all;
    using super::sidx_incident_out;
    using super::sidx_incident_in;
    using super::sidx_incident_all;
    using super::observers;

    virtual
    bool
    erase(
        const Edge* const e
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
    std::unordered_map<const Vertex*, std::unordered_map<const Vertex*, const Edge*> > cidx_edge_by_vertexes;

};

}
}

#endif
