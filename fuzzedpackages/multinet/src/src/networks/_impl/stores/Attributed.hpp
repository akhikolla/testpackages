/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_DATASTRUCTURES_STORES_ATTRIBUTED_H_
#define UU_NET_DATASTRUCTURES_STORES_ATTRIBUTED_H_

#include <memory>

namespace uu {
namespace net {

/**
 * Publicly inheriting from this class enables the usage of an attribute store.
 * In this way containers can internally store attributes and attribute values
 * for their elements.
 *
 * The attribute store becomes accessible through the attr method. The type of
 * store is defined by the template typename A.
 */
template <typename A>
class
    Attributed
{
  public:

    Attributed(
        std::unique_ptr<A> attr
    );

    /**
     * Returns a reference to the attribute store managing attributes for this attributed object.
     */
    A*
    attr(
    );

    /**
     * Returns a reference to the attribute store managing attributes for this attributed object.
     */
    const A*
    attr(
    ) const;

  protected:

    /** Internal attribute store, to handle attributes. */
    std::unique_ptr<A> attributes_;

};



template <typename A>
Attributed<A>::
Attributed(
    std::unique_ptr<A> attr
)
{
    attributes_ = std::move(attr);
}


template <typename A>
A*
Attributed<A>::attr(
)
{
    return attributes_.get();
}


template <typename A>
const A*
Attributed<A>::attr(
) const
{
    return attributes_.get();
}


}
}

#endif
