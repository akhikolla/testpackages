/**
 */

#ifndef UU_OLAP_SEL_SELECTION_H_
#define UU_OLAP_SEL_SELECTION_H_

#include <stddef.h>

namespace uu {
namespace core {
namespace sel {

class Indexes
{

  public:

    Indexes();

    /**  */
    virtual
    void
    eval(
        size_t size
    );

    /**  */
    virtual
    bool
    has_next(
    ) const;


    /**  */
    virtual
    size_t
    next(
    );

};

}
}
}

#endif
