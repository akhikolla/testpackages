/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_CORE_DATASTRUCTURES_OBSERVERS_OBSERVERSTORE_H_
#define UU_CORE_DATASTRUCTURES_OBSERVERS_OBSERVERSTORE_H_


#include <vector>
#include <memory>
#include "core/datastructures/observers/GenericObserver.hpp"

namespace uu {
namespace core {

/**
 * Observer/Subject pattern: a Subject can store several Observers, so that they can be
 * notified when something happens (typically an erase or add event).
 */
class ObserverStore
{
  public:

    void
    register_observer(
        std::unique_ptr<core::GenericObserver> obs
    );

  protected:

    std::vector<std::unique_ptr<GenericObserver>> observers_;

};

}
}

#endif
