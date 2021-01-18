#include "core/olap/selection/IndexIterator.hpp"
#include <iostream>

namespace uu {
namespace core {
namespace sel {


IndexIterator::
IndexIterator(
    const std::vector<size_t>& size
)
{
    std::vector<std::vector<size_t>> indexes;

    for (size_t i=0; i<size.size(); i++)
    {
        indexes.push_back(std::vector<size_t>());

        for (size_t j=0; j<size.at(i); j++)
        {
            indexes.at(i).push_back(j);
        }
    }

    indexes_ = indexes;
}


IndexIterator::
IndexIterator(
    const std::vector<std::vector<size_t>>& indexes
) : indexes_(indexes)
{
}



IndexIterator::iterator
IndexIterator::
begin(
) const
{
    return iterator(indexes_);
}

IndexIterator::iterator
IndexIterator::
end(
) const
{
    return iterator(indexes_, {});
}

IndexIterator::iterator::
iterator(
    const std::vector<std::vector<size_t>>& indexes,
    const std::vector<size_t>& current
) : indexes_(indexes), current_(current)
{
}


IndexIterator::iterator::
iterator(
    const std::vector<std::vector<size_t>>& indexes
) : indexes_(indexes)
{
    // check there is no empty set of selectors
    for (auto dim: indexes)
    {
        if (dim.size() == 0)
        {
            current_ = {};
            return;
        }
    }

    // set the first index
    for (size_t i = 0; i < indexes_.size(); i++)
    {
        current_.push_back(0);
    }
}


std::vector<size_t>
IndexIterator::iterator::
operator*(
)
{
    std::vector<size_t> res;

    for (size_t idx = 0; idx < indexes_.size(); idx++)
    {
        auto val = indexes_.at(idx).at(current_.at(idx));
        res.push_back(val);
    }

    return res;
}


IndexIterator::iterator
IndexIterator::iterator::
operator++(
)
{
    // PREFIX

    for (size_t idx = 0; idx < indexes_.size(); idx++)
    {
        if (current_.at(idx) < indexes_.at(idx).size() - 1)
        {
            current_[idx]++;
            break;
        }

        else
        {
            if (idx == indexes_.size() - 1)
            {
                current_ = {};
                break;
            }

            current_[idx] = 0;

        }
    }

    return *this;
}


IndexIterator::iterator
IndexIterator::iterator::
operator++(
    int
)
{

    // POSTFIX
    auto tmp = IndexIterator::iterator(indexes_, current_);

    ++(*this);

    return tmp;
}


bool
IndexIterator::iterator::
operator==(
    const IndexIterator::iterator& rhs
)
{
    return (current_ == rhs.current_);
}


bool
IndexIterator::iterator::
operator!=(
    const IndexIterator::iterator& rhs
)
{
    return (current_ != rhs.current_);
}


}
}
}

