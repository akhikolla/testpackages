/**
 * History:
 * - 2019.07.21 File created
 */

#ifndef UU_MNET_OLAP_ECUBE_H_
#define UU_MNET_OLAP_ECUBE_H_

#include <memory>
#include <string>
#include <unordered_map>
#include "core/olap/datastructures/CCube.hpp"
#include "networks/_impl/stores/MDSimpleEdgeStore.hpp"
#include "networks/_impl/olap/VCube.hpp"
#include "objects/Vertex.hpp"
#include "objects/MLEdge.hpp"

namespace uu {
namespace net {

/**
 * Similar to an edge store, but updates happen at cell level
 */
class
    ECube
{

  private:

    std::unique_ptr<core::CCube<MDSimpleEdgeStore<VCube>>> cube_;

    std::string name_;

    typedef MLEdge<Vertex, VCube> IEdge;

  public:

    // ECube(const std::vector<size_t>& dim);

    ECube(
        const std::string& name,
        VCube* vc1,
        VCube* vc2,
        EdgeDir dir,
        const std::vector<std::string>& dim,
        const std::vector<std::vector<std::string>>& members
    );

    /** Returns a const iterator to the first object in the cube
     @todo should this be const_iterator? */
    core::SortedRandomBag<const IEdge *>::iterator
    begin(
    ) const;

    /** Returns a const iterator after the last object in the cube
     @todo should this be const_iterator? */
    core::SortedRandomBag<const IEdge *>::iterator
    end(
    ) const;

    /**
     * Returns the number of vertices in the cube
     */
    size_t
    size(
    ) const;

    /*const core::UnionSortedRandomSet<typename CONTAINER_TYPE::value_type>*
    elements(
    ) const;*/

    /** Returns true if an object with the input id is present in the collection */
    bool
    contains(
        const IEdge* v
    ) const;

    const IEdge*
    get(
        const Vertex* v1,
        const VCube* l1,
        const Vertex* v2,
        const VCube* l2
    ) const;

    /** Returns the object at the given position in the collection.
     * @throw ElementNotFoundException if the index is outside the bounds on the set
     */
    const IEdge*
    at(
        size_t pos
    ) const;

    /** Returns a random object, uniform probability */
    const IEdge*
    get_at_random(
    ) const;


    /** Returns the position of the input value in the collection, or -1 */
    int
    index_of(
        const IEdge* v
    ) const;


    core::AttributeStore<IEdge>*
    attr(
    );


    const core::AttributeStore<IEdge>*
    attr(
    ) const;


    /**
     * Returns the order (number of dimensions) of this cube.
     */
    size_t
    order(
    ) const;


    /**
     * Returns the dimensions of this cube.
     */
    const std::vector<std::string>&
    dim(
    ) const;


    /**
     * Returns the members of a dimension.
     */
    const std::vector<std::string>&
    members(
        const std::string& dim
    ) const;


    /**
     * Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    MDSimpleEdgeStore<VCube>*
    operator[](
        const std::vector<size_t>& index
    );

    /**
     * Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    MDSimpleEdgeStore<VCube>*
    operator[](
        const std::vector<std::string>& index
    );

    /**
     * Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    const MDSimpleEdgeStore<VCube>*
    operator[](
        const std::vector<size_t>& index
    ) const;

    /**
     * Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    const MDSimpleEdgeStore<VCube>*
    operator[](
        const std::vector<std::string>& index
    ) const;


    /** Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    MDSimpleEdgeStore<VCube>*
    at(
        const std::vector<size_t>& index
    );

    /** Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    MDSimpleEdgeStore<VCube>*
    at(
        const std::vector<std::string>& index
    );

    /** Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    const MDSimpleEdgeStore<VCube>*
    at(
        const std::vector<size_t>& index
    ) const;

    /** Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    const MDSimpleEdgeStore<VCube>*
    at(
        const std::vector<std::string>& index
    ) const;

    /*
    friend
    std::unique_ptr<CCube<CONTAINER_TYPE>>
    vslice(
           CCube<CONTAINER_TYPE>* cube,
           const std::vector<std::vector<size_t>>& indexes
           );
    */

    std::string
    to_string(
    ) const;

    void
    attach(
        core::Observer<const IEdge>* obs
    );

    /*
    static
    std::unique_ptr<ECube>
    create(
        const std::string& name,
        const std::vector<std::string>& dim,
        const std::vector<std::vector<std::string>>& members
    );*/


    bool
    is_directed(
    );


  private:

    /** Edge directionality */
    EdgeDir edge_directionality;

    std::vector<size_t>
    index_of(
        const std::vector<std::string>& pos
    ) const;

};



}
}

#endif
