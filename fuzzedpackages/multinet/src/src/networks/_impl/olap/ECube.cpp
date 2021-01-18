/**
 * History:
 * - 2019.07.21 File created
 */

#include "core/exceptions/ElementNotFoundException.hpp"
#include "networks/_impl/olap/ECube.hpp"

namespace uu {
namespace net {

ECube::
ECube(
    const std::string& name,
    VCube* vc1,
    VCube* vc2,
    EdgeDir dir,
    const std::vector<std::string>& dimensions,
    const std::vector<std::vector<std::string>>& members

)
{
    name_ = name;

    edge_directionality = dir;

    std::vector<std::shared_ptr<MDSimpleEdgeStore<VCube>>> elements;
    size_t cube_size = 1;

    for (auto dim: members)
    {
        cube_size *= dim.size();
    }

    for (size_t i = 0; i < cube_size; i++)
    {
        elements.push_back(std::make_shared<MDSimpleEdgeStore<VCube>>(vc1, vc2, dir));
    }

    cube_ = std::make_unique<core::CCube<MDSimpleEdgeStore<VCube>>>(dimensions, members, elements.begin(), elements.end());
}

core::SortedRandomBag<const MLEdge<Vertex, VCube>*>::iterator
ECube::
begin(
) const
{
    return cube_->elements()->begin();
}

core::SortedRandomBag<const MLEdge<Vertex, VCube> *>::iterator
ECube::
end(
) const
{
    return cube_->elements()->end();
}


size_t
ECube::
size(
) const
{
    return cube_->elements()->size();
}

bool
ECube::
contains(
    const MLEdge<Vertex, VCube>* v
) const
{
    return cube_->elements()->contains(v);
}


const MLEdge<Vertex, VCube>*
ECube::
get(
    const Vertex* v1,
    const VCube* l1,
    const Vertex* v2,
    const VCube* l2
) const
{
    return cube_->elements()->get(std::tuple<const Vertex*, const VCube*, const Vertex*, const VCube*>(v1,l1,v2,l2));
}

const MLEdge<Vertex, VCube>*
ECube::
at(
    size_t pos
) const
{
    return cube_->elements()->at(pos);
}

const MLEdge<Vertex, VCube>*
ECube::
get_at_random(
) const
{
    return cube_->elements()->get_at_random();
}


int
ECube::
index_of(
    const MLEdge<Vertex, VCube>* v
) const
{
    return cube_->elements()->index_of(v);
}


core::AttributeStore<MLEdge<Vertex, VCube>>*
        ECube::
        attr(
        )
{
    return cube_->attr();
}


const core::AttributeStore<MLEdge<Vertex, VCube>>*
        ECube::
        attr(
        ) const
{
    return cube_->attr();
}


size_t
ECube::
order(
) const
{
    return cube_->order();
}


const std::vector<std::string>&
ECube::
dim(
) const
{
    return cube_->dim();
}


const std::vector<std::string>&
ECube::
members(
    const std::string& dim
) const
{
    return cube_->members(dim);
}


MDSimpleEdgeStore<VCube>*
ECube::
operator[](
    const std::vector<size_t>& index
)
{
    return cube_->operator[](index);
}


MDSimpleEdgeStore<VCube>*
ECube::
operator[](
    const std::vector<std::string>& index
)
{
    return cube_->operator[](index);
}


const MDSimpleEdgeStore<VCube>*
ECube::
operator[](
    const std::vector<size_t>& index
) const
{
    return cube_->operator[](index);
}


const MDSimpleEdgeStore<VCube>*
ECube::
operator[](
    const std::vector<std::string>& index
) const
{
    return cube_->operator[](index);
}


MDSimpleEdgeStore<VCube>*
ECube::
at(
    const std::vector<size_t>& index
)
{
    return cube_->at(index);
}


MDSimpleEdgeStore<VCube>*
ECube::
at(
    const std::vector<std::string>& index
)
{
    return cube_->at(index);
}

const MDSimpleEdgeStore<VCube>*
ECube::
at(
    const std::vector<size_t>& index
) const
{
    return cube_->at(index);
}


const MDSimpleEdgeStore<VCube>*
ECube::
at(
    const std::vector<std::string>& index
) const
{
    return cube_->at(index);
}

std::string
ECube::
to_string(
) const
{
    return name_;
}


bool
ECube::
is_directed(
)
{
    return edge_directionality==EdgeDir::DIRECTED?true:false;
}


void
ECube::
attach(
    core::Observer<const MLEdge<Vertex, VCube>>* obs
)
{
    cube_->elements()->attach(obs);
}

/*
std::unique_ptr<ECube>
ECube::
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

std::vector<const std::shared_ptr<MDSimpleEdgeStore<VCube>>> stores;

for (size_t i = 0; i < num_entries; i++)
{
    stores.push_back(std::make_shared<MDSimpleEdgeStore<VCube>>());
}

return std::make_unique<ECube>(name, dim, members, stores.begin(), stores.end());
}
*/

}
}

