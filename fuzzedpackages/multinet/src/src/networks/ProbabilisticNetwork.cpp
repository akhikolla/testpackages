#include "networks/ProbabilisticNetwork.hpp"

#include "core/exceptions/WrongParameterException.hpp"

namespace uu {
namespace net {

ProbabilisticNetwork::
ProbabilisticNetwork(
    const std::string& name,
    EdgeDir dir,
    bool allows_loops
) : super(name, dir, allows_loops)
{
    auto p_attr = core::Attribute::create(kPROB_ATTR_NAME, core::AttributeType::DOUBLE);

    edges()->attr()->add(std::move(p_attr));

}


void
ProbabilisticNetwork::
set_prob(
    const Edge* e,
    double p
)
{
    if (p<0 || p>1)
    {
        throw core::WrongParameterException("Edge probabilities must be between 0 and 1");
    }

    edges()->attr()->set_double(e, kPROB_ATTR_NAME, p);
}


core::Value<double>
ProbabilisticNetwork::
get_prob(
    const Edge* e
) const
{
    return edges()->attr()->get_double(e, kPROB_ATTR_NAME);
}


bool
ProbabilisticNetwork::
is_probabilistic(
) const
{
    return true;
}


}
}

