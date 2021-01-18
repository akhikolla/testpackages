/**
 * History:
 * - 2019.07.21 File created
 */

#include "core/exceptions/ElementNotFoundException.hpp"
#include "networks/_impl/olap/VCube.hpp"
#include "networks/_impl/olap/slice.hpp"

namespace uu {
namespace net {

VCube::
VCube(
    const std::string& name,
    const std::vector<std::string>& dimensions,
    const std::vector<std::vector<std::string>>& members
)
{
    name_ = name;

    std::vector<std::shared_ptr<VertexStore>> elements;
    size_t cube_size = 1;

    for (auto dim: members)
    {
        cube_size *= dim.size();
    }

    for (size_t i = 0; i < cube_size; i++)
    {
        elements.push_back(std::make_shared<VertexStore>());
    }

    cube_ = std::make_unique<core::CCube<VertexStore>>(dimensions, members, elements.begin(), elements.end());
}


template <typename Iterator>
VCube::
VCube(
    const std::string& name,
    const std::vector<std::string>& dimensions,
    const std::vector<std::vector<std::string>>& members,
    Iterator begin,
    Iterator end
)
{
    // @todo check Iterator type

    name_ = name;

    //cube_ = std::make_unique<core::CCube<VertexStore>>(dim, members, begin, end);

    std::vector<std::shared_ptr<VertexStore>> elements;
    size_t cube_size = 1;

    for (auto dim: members)
    {
        cube_size *= dim.size();
    }

    for (size_t i = 0; i < cube_size; i++)
    {
        elements.push_back(std::make_shared<VertexStore>());
    }

    cube_ = std::make_unique<core::CCube<VertexStore>>(dimensions, members, begin, end);
}

core::SortedRandomBag<const Vertex *>::iterator
VCube::
begin(
) const
{
    return cube_->elements()->begin();
}

core::SortedRandomBag<const Vertex *>::iterator
VCube::
end(
) const
{
    return cube_->elements()->end();
}


size_t
VCube::
size(
) const
{
    return cube_->elements()->size();
}

bool
VCube::
contains(
    const Vertex* v
) const
{
    return cube_->elements()->contains(v);
}


const Vertex*
VCube::
get(
    const std::string& key
) const
{
    return cube_->elements()->get(key);
}

const Vertex*
VCube::
at(
    size_t pos
) const
{
    return cube_->elements()->at(pos);
}

const Vertex*
VCube::
get_at_random(
) const
{
    return cube_->elements()->get_at_random();
}


int
VCube::
index_of(
    const Vertex* v
) const
{
    return cube_->elements()->index_of(v);
}


core::AttributeStore<Vertex>*
VCube::
attr(
)
{
    return cube_->attr();
}


const core::AttributeStore<Vertex>*
VCube::
attr(
) const
{
    return cube_->attr();
}


size_t
VCube::
order(
) const
{
    return cube_->order();
}


const std::vector<std::string>&
VCube::
dim(
) const
{
    return cube_->dim();
}


const std::vector<std::string>&
VCube::
members(
    const std::string& dim
) const
{
    return cube_->members(dim);
}


VertexStore*
VCube::
operator[](
    const std::vector<size_t>& index
)
{
    return cube_->operator[](index);
}


VertexStore*
VCube::
operator[](
    const std::vector<std::string>& index
)
{
    return cube_->operator[](index);
}


const VertexStore*
VCube::
operator[](
    const std::vector<size_t>& index
) const
{
    return cube_->operator[](index);
}


const VertexStore*
VCube::
operator[](
    const std::vector<std::string>& index
) const
{
    return cube_->operator[](index);
}


VertexStore*
VCube::
at(
    const std::vector<size_t>& index
)
{
    return cube_->at(index);
}


VertexStore*
VCube::
at(
    const std::vector<std::string>& index
)
{
    return cube_->at(index);
}

const VertexStore*
VCube::
at(
    const std::vector<size_t>& index
) const
{
    return cube_->at(index);
}


const VertexStore*
VCube::
at(
    const std::vector<std::string>& index
) const
{
    return cube_->at(index);
}

std::string
VCube::
to_string(
) const
{
    return name_;
}


void
VCube::
attach(
    core::Observer<const Vertex>* obs
)
{
    cube_->elements()->attach(obs);
}

/*
std::unique_ptr<VCube>
VCube::
create(
const std::string& name,
const std::vector<std::string>& dim,
const std::vector<std::vector<std::string>>& members
)
{
size_t num_entries = 1;

for (auto m: members)
{
    num_entries *= m.size();
}

std::vector<const std::shared_ptr<VertexStore>> stores;

for (size_t i = 0; i < num_entries; i++)
{
    stores.push_back(std::make_shared<VertexStore>());
}

return std::make_unique<VCube>(name, dim, members, stores.begin(), stores.end());
}
*/


void
VCube::
resize(
    const std::string& dimension,
    const std::string& member
)
{
    std::vector<std::shared_ptr<VertexStore>> elements;
    size_t num_new_elements = 1;

    for (auto d: dim())
    {
        if (d == dimension)
        {
            continue;
        }

        num_new_elements *= members(d).size();
    }

    for (size_t i = 0; i < num_new_elements; i++)
    {
        elements.push_back(std::make_shared<VertexStore>());
    }

    cube_->resize(dimension, member, elements.begin(), elements.end());
}


std::unique_ptr<VCube>
VCube::
vslice(
    const std::string& name,
    const std::vector<std::vector<size_t>>& indexes
)
{

    // get dimensions and members of the new cube

    auto dim_names = dim();

    std::vector<std::vector<std::string>> member_names(dim_names.size());

    for (size_t i = 0; i < dim_names.size(); i++)
    {
        auto m = members(dim_names.at(i));

        for (auto idx: indexes[i])
        {
            member_names[i].push_back(m.at(idx));

        }
    }


    auto iter = core::sel::IndexIterator(indexes);

    std::vector<std::shared_ptr<VertexStore>> elements;

    for (auto idx: iter)
    {
        auto vs = at(idx);
        elements.push_back(vs->shared_from_this());
    }

    auto res = std::make_unique<VCube>(name, dim_names, member_names, elements.begin(), elements.end());

    return res;

}

}
}

