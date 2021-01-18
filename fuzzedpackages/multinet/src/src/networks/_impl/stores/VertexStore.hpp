/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_STORES_VERTEXSTORE_H_
#define UU_NET_STORES_VERTEXSTORE_H_

#include <memory>
#include "core/datastructures/observers/Observer.hpp"
#include "objects/Vertex.hpp"
#include "core/stores/ObjectStore.hpp"

namespace uu {
namespace net {

/**
 * A VertexStore allows to create, store, retrieve and erase a set of vertices.
 */
class
    VertexStore :
    public std::enable_shared_from_this<VertexStore>
{

  private:

    std::unique_ptr<core::ObjectStore<Vertex>> store_;

  public:

    typedef Vertex value_type;
    typedef std::string key_type;

    VertexStore(
    );


    /** Returns an iterator to the first object in the collection */
    core::PtrSortedRandomSet<const Vertex,std::shared_ptr<const Vertex>,core::SharedPtrLT<const Vertex>,core::SharedPtrEQ<const Vertex>>::iterator
            begin(
            ) const;

    /** Returns an iterator after the last object in the collection */
    core::PtrSortedRandomSet<const Vertex,std::shared_ptr<const Vertex>,core::SharedPtrLT<const Vertex>,core::SharedPtrEQ<const Vertex>>::iterator
            end(
            ) const;

    /** Returns the number of objects in the collection */
    size_t
    size(
    ) const;

    /**
     * Inserts a new object in the collection.
     * @ret
     */
    const Vertex*
    add(
        std::shared_ptr<const Vertex> v
    );

    /**
     * Inserts a new object in the collection.
     * @ret
     */
    const Vertex*
    add(
        const Vertex* v
    );


    /**
     * Creates a new vertex and adds it to the store.
     */
    const Vertex*
    add(
        const std::string& key
    );

    /** Returns true if an object with the input id is present in the collection */
    bool
    contains(
        const Vertex* v
    ) const;

    const Vertex*
    get(
        const std::string& key
    ) const;

    /** Returns the object at the given position in the collection.
     * @throw ElementNotFoundException if the index is outside the bounds on the set
     */
    const Vertex*
    at(
        size_t pos
    ) const;

    /** Returns a random object, uniform probability */
    const Vertex*
    get_at_random(
    ) const;


    /** Returns the position of the input value in the collection, or -1 */
    int
    index_of(
        const Vertex* v
    ) const;


    bool
    erase(
        const Vertex * v
    );

    void
    attach(
        core::Observer<const Vertex>* obs
    );

    /**
     * Returns a short string summary of this store, for example including
     * the number of objects it contains.
     */
    virtual
    std::string
    summary(
    ) const;

};


}
}

#endif
