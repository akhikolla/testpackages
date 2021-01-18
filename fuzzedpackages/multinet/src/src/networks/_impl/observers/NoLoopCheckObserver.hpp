/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_DATASTRUCTURE_OBSERVERS_NOLOOPCHECKOBSERVER_H_
#define UU_NET_DATASTRUCTURE_OBSERVERS_NOLOOPCHECKOBSERVER_H_

#include "core/exceptions/NullPtrException.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/datastructures/observers/Observer.hpp"
#include "objects/Edge.hpp"

namespace uu {
namespace net {

/**
 * This observer checks that the vertices at the end of new edges are not the same (that is, there are no loops).
 */
class NoLoopCheckObserver :
    public core::Observer<const Edge>
{

  public:

    /**
     * Creates an observer with a pointer to the vertex store where vertices are expected to be.
     *
     */
    NoLoopCheckObserver();

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

};


}
}

#endif
