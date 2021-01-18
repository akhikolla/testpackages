/**
 * History:
 * - 2018.03.09 file created
 */

#ifndef UU_CORE_DATASTRUCTURES_OBSERVERS_GENERICOBSERVER_H_
#define UU_CORE_DATASTRUCTURES_OBSERVERS_GENERICOBSERVER_H_

namespace uu {
namespace core {

/**
 * This is a generic observer, as used in the observer design pattern.
 * All Observers inherit from this class.
 * In this way, complex objects can store several observers of different types in a uniform way.
 */
class GenericObserver
{
  public:

    virtual
    ~GenericObserver() = default;

  protected:

    GenericObserver() {};

};



}
}

#endif
