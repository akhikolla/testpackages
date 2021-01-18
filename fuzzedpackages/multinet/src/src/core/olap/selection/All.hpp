/**
 */

#ifndef UU_OLAP_SEL_RANGE_H_
#define UU_OLAP_SEL_RANGE_H_

#include "core/olap/selection/Indexes.hpp"
#include <vector>

namespace uu {
namespace core {
namespace sel {

class All
    : public Indexes
{

  public:

    All();

    /**  */
    virtual
    void
    eval(
        size_t size
    ) override;

    /**  */
    virtual
    bool
    has_next(
    ) const override;


    /**  */
    virtual
    size_t
    next(
    ) override;

  private:

    size_t max_;
    size_t current_;
    bool has_next_;

};

}
}
}

#endif
