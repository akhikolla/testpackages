/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_CORE_DATASTRUCTURES_OBSERVERS_SUBJECT_H_
#define UU_CORE_DATASTRUCTURES_OBSERVERS_SUBJECT_H_

#include <vector>
#include "core/datastructures/observers/Observer.hpp"
#include "core/exceptions/assert_not_null.hpp"

namespace uu {
namespace core {

/**
 * Observer/Subject pattern: a Subject can store several Observers, so that they can be
 * notified when something happens (typically an erase or add event).
 */
template <typename OT>
class
    Subject
{

  public:

    /**
     * Registers an observer that must be notified when vertices are erased, so that
     * they can propagate the action.
     *
     * For example, the attribute store will need to erase the attributes associated
     * to the erased vertex.
     */
    void
    attach(
        core::Observer<OT>* observer
    );

  protected:

    /* A vector of objects to be notified when a vertex is erased, so that
     * the observer can propagate the action to other stores.
     * For example, removing edges adjacent to the erased vertex.
     */
    std::vector<core::Observer<OT>*> observers;

};


template <typename OT>
void
Subject<OT>::
attach(
    core::Observer<OT>* obs
)
{
    assert_not_null(obs, "Subject::attach", "obs");
    observers.push_back(obs);
}

}
}

#endif
