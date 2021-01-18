/**
 */

#ifndef UU_OLAP_SEL_RANGE_H_
#define UU_OLAP_SEL_RANGE_H_

#include "core/olap/selection/Indexes.hpp"
#include <vector>

namespace uu {
namespace core {
namespace sel {

class Range
    : public Indexes
{

  public:

    Range(size_t from, size_t to);

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

    size_t from_;
    size_t to_;
    size_t current_;
    bool has_next_;

};

}
}
}

#endif
