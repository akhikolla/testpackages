#include "networks/WeightedNetwork.hpp"

namespace uu {
namespace net {

WeightedNetwork::
WeightedNetwork(
    const std::string& name,
    EdgeDir dir,
    bool allows_loops
) : super(name, dir, allows_loops)
{

    auto w_attr = core::Attribute::create(kWEIGHT_ATTR_NAME, core::AttributeType::DOUBLE);

    edges()->attr()->add(std::move(w_attr));

}


void
WeightedNetwork::
set_weight(
    const Edge* e,
    double w
)
{
    edges()->attr()->set_double(e, kWEIGHT_ATTR_NAME, w);
}


core::Value<double>
WeightedNetwork::
get_weight(
    const Edge* e
) const
{
    return edges()->attr()->get_double(e, kWEIGHT_ATTR_NAME);
}


bool
WeightedNetwork::
is_weighted(
) const
{
    return true;
}


}
}

