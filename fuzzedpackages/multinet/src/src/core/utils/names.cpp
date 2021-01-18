#include "core/utils/names.hpp"

#include <iostream>
#include <sstream>

namespace uu {
namespace core {

NameIterator::
NameIterator(
    const std::string& prefix,
    size_t num_names
) :
    prefix_(prefix),
    num_names_(num_names)
{
    num_digits_ = 1;

    num_names -= 1;
    num_names /= 10;

    while (num_names > 0)
    {
        num_digits_++;
        num_names /= 10;
    }

}


NameIterator::iterator::
iterator(
    const std::string& prefix,
    size_t num_names,
    size_t num_digits
):
    prefix_(prefix),
    num_names_(num_names),
    num_digits_(num_digits)
{
    current_ = 0;
}

NameIterator::iterator::
iterator(
    const std::string& prefix,
    size_t num_names,
    size_t num_digits,
    size_t current
):
    prefix_(prefix),
    num_names_(num_names),
    current_(current),
    num_digits_(num_digits)
{}


std::string
NameIterator::iterator::
operator*(
)
{

    std::string c = std::to_string(current_);
    size_t length = c.length();

    std::stringstream ss;

    ss << prefix_;

    for (size_t i = 0; i < num_digits_ - length; i++)
    {
        ss << "0";
    }

    ss << c;

    return ss.str();
}

NameIterator::iterator
NameIterator::iterator::
operator++(
)
{
    if (current_<= num_names_)
    {
        current_++;
    }

    return *this;
}

NameIterator::iterator
NameIterator::iterator::
operator++(
    int
)
{
    auto tmp = NameIterator::iterator(prefix_, num_names_, num_digits_, current_);

    ++(*this);

    return tmp;
}

bool
NameIterator::iterator::
operator==(
    const NameIterator::iterator& rhs
)
{
    return (prefix_ == rhs.prefix_
            && num_digits_ == rhs.num_digits_
            && current_ == rhs.current_);
}

bool
NameIterator::iterator::
operator!=(
    const NameIterator::iterator& rhs
)
{
    return (prefix_ != rhs.prefix_
            || num_digits_ != rhs.num_digits_
            || current_ != rhs.current_);
}

NameIterator::
iterator
NameIterator::
begin(
) const
{
    return NameIterator::iterator(prefix_, num_names_, num_digits_);
}

NameIterator::iterator
NameIterator::
end(
) const
{
    return NameIterator::iterator(prefix_, num_names_, num_digits_, num_names_);
}

}
}

