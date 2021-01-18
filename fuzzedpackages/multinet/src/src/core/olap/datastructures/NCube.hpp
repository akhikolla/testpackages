/**
 *
 * History:
 * - 2019.08.01 file created
 */

#ifndef UU_OLAP_NCUBE_H_
#define UU_OLAP_NCUBE_H_

#include "core/exceptions/ElementNotFoundException.hpp"
#include "core/exceptions/WrongParameterException.hpp"
#include "core/olap/datastructures/Cube.hpp"
#include <string>
#include <unordered_map>
#include <vector>

namespace uu {
namespace core {

template <class ENTRY_TYPE>
class NCube
    : public Cube<ENTRY_TYPE>
{
  private:

    using super = Cube<ENTRY_TYPE>;

  protected:

    using super::size_;
    using super::data_;

  public:

    /**
     * Creates a cube specifying its dimensions and the members for each dimension.
     */
    NCube(
        const std::vector<std::string>& dim,
        const std::vector<std::vector<std::string>>& members
    );

    /**
     * Creates a cube specifying its dimensions and the members for each dimension.
     */
    template <class Iterator>
    NCube(
        const std::vector<std::string>& dim,
        const std::vector<std::vector<std::string>>& members,
        Iterator begin,
        Iterator end
    );

    /**
     * Destroys the cube.
     */
    ~NCube()
    {

    }

    /**
     * Returns the dimensions of this cube.
     */
    const std::vector<std::string>&
    dim(
    ) const;

    /**
     * Returns the members of a dimension.
     */
    const std::vector<std::string>&
    members(
        const std::string& dim
    ) const;

    using super::operator[];

    /**
     * Returns the object at the given index in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    ENTRY_TYPE&
    operator[](
        const std::vector<std::string>& index
    );

    /**
     * Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    const ENTRY_TYPE&
    operator[](
        const std::vector<std::string>& index
    ) const;

    using super::at;

    /**
     * Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    ENTRY_TYPE&
    at(
        const std::vector<std::string>& index
    );


    /** Returns the object at the given position in the cube.
     * @throw ElementNotFoundException if the index is outside the bounds on the cube
     */
    const ENTRY_TYPE&
    at(
        const std::vector<std::string>& index
    ) const;

    /* index of a dimension */
    size_t
    index_of(
        const std::string& dim
    ) const;

    /* index of a member given a dimension */
    size_t
    index_of(
        const std::string& dim,
        const std::string& member
    ) const;

    /* computes a numerical index from a label-based index */
    std::vector<size_t>
    index_of(
        const std::vector<std::string>& index
    ) const;

    std::vector<size_t>
    size_of(
        const std::vector<std::vector<std::string>>& members
    ) const;

  private:

    /* Dimensions */
    std::vector<std::string> dim_;
    /* Index of each dimension (by name) */
    std::unordered_map<std::string,size_t> dim_idx_;

    /* Members, for each dimension */
    std::vector<std::vector<std::string>> members_;
    /* Index of each member (by name), for each dimension */
    std::vector<std::unordered_map<std::string, size_t>> members_idx_;


};

/* TEMPLATE CODE */

template <class ENTRY_TYPE>
NCube<ENTRY_TYPE>::
NCube(
    const std::vector<std::string>& dim,
    const std::vector<std::vector<std::string>>& members
) : super(size_of(members)), dim_(dim), members_(members)
{
    if (dim.size() != members.size())
    {
        throw core::WrongParameterException("parameters dim and members must have the same size");
    }

    size_t order = dim_.size();

    for (size_t i = 0; i < order; i++)
    {
        dim_idx_[dim_[i]] = i;
        members_idx_.push_back(std::unordered_map<std::string, size_t>());

        for (size_t j = 0; j < members[i].size(); j++)
        {
            auto member = members[i][j];
            members_idx_[i][member] = j;
        }
    }

    std::vector<size_t> size(order);

    for (size_t i = 0; i < order; i++)
    {
        size[i] = members_[i].size();
    }

    //super::resize(size);
}


template <class ENTRY_TYPE>
template <class Iterator>
NCube<ENTRY_TYPE>::
NCube(
    const std::vector<std::string>& dim,
    const std::vector<std::vector<std::string>>& members,
    Iterator begin,
    Iterator end
) : super(size_of(members), begin, end), dim_(dim), members_(members)
{
    // @todo remove code duplication if possible
    if (dim.size() != members.size())
    {
        throw core::WrongParameterException("parameters dim and members must have the same size");
    }

    size_t order = dim_.size();

    for (size_t i = 0; i < order; i++)
    {
        dim_idx_[dim_[i]] = i;
        members_idx_.push_back(std::unordered_map<std::string, size_t>());

        for (size_t j = 0; j < members[i].size(); j++)
        {
            auto member = members[i][j];
            members_idx_[i][member] = j;
        }
    }

    std::vector<size_t> size(order);

    for (size_t i = 0; i < order; i++)
    {
        size[i] = members_[i].size();
    }
}


template <class ENTRY_TYPE>
const std::vector<std::string>&
NCube<ENTRY_TYPE>::
dim(
) const
{
    return dim_;
}


template <class ENTRY_TYPE>
const std::vector<std::string>&
NCube<ENTRY_TYPE>::
members(
    const std::string& dim
) const
{
    auto f = dim_idx_.find(dim);

    if (f != dim_idx_.end())
    {
        // no need to check bounds
        return members_[f->second];
    }

    else
    {
        throw core::ElementNotFoundException("dimension");
    }
}


template <class ENTRY_TYPE>
ENTRY_TYPE&
NCube<ENTRY_TYPE>::
operator[](
    const std::vector<std::string>& index
)
{
    auto idx = index_of(index);
    return super::operator[](idx);
}

template <class ENTRY_TYPE>
const ENTRY_TYPE&
NCube<ENTRY_TYPE>::
operator[](
    const std::vector<std::string>& index
) const
{
    auto idx = index_of(index);
    return super::operator[](idx);
}


template <class ENTRY_TYPE>
ENTRY_TYPE&
NCube<ENTRY_TYPE>::
at(
    const std::vector<std::string>& index
)
{
    auto idx = index_of(index);
    return super::at(idx);
}

template <class ENTRY_TYPE>
const ENTRY_TYPE&
NCube<ENTRY_TYPE>::
at(
    const std::vector<std::string>& index
) const
{
    auto idx = index_of(index);
    return super::at(idx);
}


template <class ENTRY_TYPE>
std::vector<size_t>
NCube<ENTRY_TYPE>::
index_of(
    const std::vector<std::string>& members
) const
{

    if (dim_.size() != members.size())
    {
        throw core::WrongParameterException("parameter must have the same size as the cube order");
    }

    std::vector<size_t> res;

    for (size_t i = 0; i < members.size(); i++)
    {
        auto member = members[i];
        auto f = members_idx_[i].find(member);

        if (f != members_idx_[i].end())
        {
            res.push_back(f->second);
        }

        else
        {
            throw core::ElementNotFoundException("member " + member);
        }
    }

    return res;
}

/* index of a dimension */
template <class ENTRY_TYPE>
size_t
NCube<ENTRY_TYPE>::
index_of(
    const std::string& dim
) const
{
    return dim_idx_.at(dim);
}


/* index of a member given a dimension */
template <class ENTRY_TYPE>
size_t
NCube<ENTRY_TYPE>::
index_of(
    const std::string& dim,
    const std::string& member
) const
{
    size_t dim_idx = index_of(dim);
    return members_idx_.at(dim_idx).at(member);
}


template <class ENTRY_TYPE>
std::vector<size_t>
NCube<ENTRY_TYPE>::
size_of(
    const std::vector<std::vector<std::string>>& members
) const
{
    std::vector<size_t> size;

    for (auto m: members)
    {
        size.push_back(m.size());
    }

    return size;
}

}
}

#endif

