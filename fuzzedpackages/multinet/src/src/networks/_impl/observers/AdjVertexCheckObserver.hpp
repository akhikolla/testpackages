/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_DATASTRUCTURE_OBSERVERS_ADJVERTEXCHECKOBSERVER_H_
#define UU_NET_DATASTRUCTURE_OBSERVERS_ADJVERTEXCHECKOBSERVER_H_

#include "core/exceptions/NullPtrException.hpp"
#include "core/datastructures/observers/Observer.hpp"
#include "objects/Edge.hpp"

namespace uu {
namespace net {

/**
 * This observer checks if the vertices adjacent to new edges exist.
 */
template<typename V>
class AdjVertexCheckObserver :
    public core::Observer<const Edge>
{

  public:

    /**
     * Creates an observer with a pointer to the vertex store where vertices are expected to be.
     *
     */
    AdjVertexCheckObserver(
        V* vertices
    );

    /**
     * Informs the observer that a vertex has been added to its vertex store.
     */
    void
    notify_add(
        const Edge* e
    ) override;


    /**
     * Informs the observer that a vertex has been erased from its vertex store.
     */
    void
    notify_erase(
        const Edge* e
    ) override;


  private:
    /** Internal vertex store. */
    V* vertices_;

};



template<typename V>
AdjVertexCheckObserver<V>::
AdjVertexCheckObserver(
    V* vertices
) :
    vertices_(vertices)
{
    if (!vertices_)
    {
        throw core::NullPtrException("vertex store to be registered in the observer");
    }
}


template<typename V>
void
AdjVertexCheckObserver<V>::
notify_add(
    const Edge* e
)
{
    if (!e)
    {
        throw core::NullPtrException("edge passed to the observer");
    }

    if (!vertices_->contains(e->v1))
    {
        throw core::ElementNotFoundException("vertex " + e->v1->name);
    }

    if (!vertices_->contains(e->v2))
    {
        throw core::ElementNotFoundException("vertex " + e->v2->name);
    }
}


template<typename V>
void
AdjVertexCheckObserver<V>::
notify_erase(
    const Edge* e
)
{
    if (!e)
    {
        throw core::NullPtrException("edge passed to the observer");
    }
}


}
}

#endif
