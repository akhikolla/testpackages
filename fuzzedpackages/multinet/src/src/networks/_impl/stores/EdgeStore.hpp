/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_DATASTRUCTURES_STORES_EDGESTORE2_H_
#define UU_NET_DATASTRUCTURES_STORES_EDGESTORE2_H_

#include <unordered_set>
#include <unordered_map>
#include "core/datastructures/containers/SharedPtrSortedRandomSet.hpp"
#include "core/datastructures/observers/Subject.hpp"
#include "objects/Edge.hpp"
#include "objects/Vertex.hpp"
#include "objects/EdgeDir.hpp"
#include "objects/EdgeMode.hpp"
#include "networks/_impl/containers/GenericObjectList.hpp"
#include "core/exceptions/OperationNotSupportedException.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace net {

/**
 * An Edge Store is a basic class to store edges. It provides functionality shared among
 * different types of edge stores.
 *
 * Once an edge store has been created with a given directionality (DIRECTED or UNDIRECTED),
 * only edges with the same directionality can be added to the store.
 *
 * E must be a class with members:
 * const Vertex* v1, const Vertex* v2, EdgeDir dir
 */
class EdgeStore :
    public core::SharedPtrSortedRandomSet<const Edge>,
    public core::Subject<const Edge>
{

  private:
    typedef core::SharedPtrSortedRandomSet<const Edge> super;

  protected:

    EdgeStore(
        EdgeDir dir
    );

  public:

    using super::add;
    using super::erase;
    using super::size;

    virtual
    const Edge*
    add(
        std::shared_ptr<const Edge> e
    ) override;

    virtual
    bool
    erase(
        const Edge* e
    ) override = 0;

    /**
     * Deletes all edges incident to a vertex.
     **/
    virtual
    void
    erase(
        const Vertex* v
    ) = 0;


    /**
     * @brief Returns the nodes with an edge from/to the input vertex.
     * @param node pointer to the node.
     * @param mode IN, OUT or INOUT.
     * @return the list of neighbors.
     **/
    const
    GenericObjectList<Vertex>*
    neighbors(
        const Vertex* vertex,
        EdgeMode mode = EdgeMode::INOUT
    ) const;

    /**
     * @brief Returns the nodes with an edge from/to the input vertex.
     * @param node pointer to the node.
     * @param mode IN, OUT or INOUT.
     * @return the list of neighbors.
     **/
    const
    GenericObjectList<Edge>*
    incident(
        const Vertex* vertex,
        EdgeMode mode = EdgeMode::INOUT
    ) const;


    bool
    is_directed(
    );

    virtual
    std::string
    summary(
    ) const;


  protected:

    /** Edge directionality */
    EdgeDir edge_directionality;

    // Indexes to sets of objects (Set IDX):
    std::unordered_map<const Vertex*, std::unique_ptr<GenericObjectList<Vertex>>> sidx_neighbors_out;
    std::unordered_map<const Vertex*, std::unique_ptr<GenericObjectList<Vertex>>> sidx_neighbors_in;
    std::unordered_map<const Vertex*, std::unique_ptr<GenericObjectList<Vertex>>> sidx_neighbors_all;

    std::unordered_map<const Vertex*, std::unique_ptr<GenericObjectList<Edge>>> sidx_incident_out;
    std::unordered_map<const Vertex*, std::unique_ptr<GenericObjectList<Edge>>> sidx_incident_in;
    std::unordered_map<const Vertex*, std::unique_ptr<GenericObjectList<Edge>>> sidx_incident_all;

};


}
}

#endif
