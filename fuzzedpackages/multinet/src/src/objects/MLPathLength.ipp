namespace uu {
namespace net {

template <typename M>
MLPathLength<M>::
MLPathLength(
    const M* mnet
) : mnet(mnet), total_length(0), ts(0) {}


template <typename M>
void
MLPathLength<M>::
step(
    const typename M::layer_type* layer1,
    const typename M::layer_type* layer2
)
{
    num_edges.inc(layer1,layer2);
    total_length++;
}

template <typename M>
long
MLPathLength<M>::
length(
) const
{
    return total_length;
}

template <typename M>
long
MLPathLength<M>::
length(
    const typename M::layer_type* layer
) const
{
    return num_edges.count(layer,layer);
}

template <typename M>
long
MLPathLength<M>::
length(
    const typename M::layer_type* from,
    const typename M::layer_type* to
) const
{
    return num_edges.count(from,to);
}

template <typename M>
ComparisonResult
MLPathLength<M>::
compare(
    const MLPathLength& other,
    ComparisonType comp
) const
{
    switch (comp)
    {
    case ComparisonType::FULL:
        return compare_full(other);

    case ComparisonType::SWITCH_COSTS:
        return compare_switch(other);

    case ComparisonType::MULTIPLEX:
        return compare_multiplex(other);

    case ComparisonType::SIMPLE:
        return compare_simple(other);
    }

    throw core::WrongParameterException("Wrong comparison type");
}

template <typename M>
ComparisonResult
MLPathLength<M>::
compare_full(
    const MLPathLength& other
) const
{
    bool canBeDominated = true;
    bool canDominate = true;

    if (mnet != other.mnet)
    {
        throw core::OperationNotSupportedException("Cannot compare distances on different networks");
    }

    for (auto layer1: *mnet->layers())
    {
        for (auto layer2: *mnet->layers())
        {
            long l1 = length(layer1,layer2);
            long l2 = other.length(layer1,layer2);

            if (l1 > l2)
            {
                canDominate = false;
            }

            else if (l1 < l2)
            {
                canBeDominated = false;
            }

            if (!canBeDominated && !canDominate)
            {
                return ComparisonResult::INCOMPARABLE;
            }
        }
    }

    if (canDominate && !canBeDominated)
    {
        return ComparisonResult::GREATER_THAN;
    }

    if (canBeDominated && !canDominate)
    {
        return ComparisonResult::LESS_THAN;
    }

    //if (canDominate && canBeDominated)
    return ComparisonResult::EQUAL;
}

template <typename M>
ComparisonResult
MLPathLength<M>::
compare_switch(
    const MLPathLength& other
) const
{
    bool canBeDominated = true;
    bool canDominate = true;

    if (mnet != other.mnet)
    {
        throw core::OperationNotSupportedException("Cannot compare distances on different networks");
    }

    long num_intralayer_steps1 = 0;
    long num_intralayer_steps2 = 0;

    for (auto layer1: *mnet->layers())
    {
        long l1 = length(layer1,layer1);
        num_intralayer_steps1 += l1;
        long l2 = other.length(layer1,layer1);
        num_intralayer_steps2 += l2;

        if (l1 > l2)
        {
            canDominate = false;
        }

        else if (l1 < l2)
        {
            canBeDominated = false;
        }

        if (!canBeDominated && !canDominate)
        {
            return ComparisonResult::INCOMPARABLE;
        }
    }

    long num_interlayer_steps1 = length()-num_intralayer_steps1;
    long num_interlayer_steps2 = other.length()-num_intralayer_steps2;

    if (num_interlayer_steps1 > num_interlayer_steps2)
    {
        canDominate = false;
    }

    else if (num_interlayer_steps1 < num_interlayer_steps2)
    {
        canBeDominated = false;
    }

    if (!canBeDominated && !canDominate)
    {
        return ComparisonResult::INCOMPARABLE;
    }

    if (canDominate && !canBeDominated)
    {
        return ComparisonResult::GREATER_THAN;
    }

    if (canBeDominated && !canDominate)
    {
        return ComparisonResult::LESS_THAN;
    }

    //if (canDominate && canBeDominated)
    return ComparisonResult::EQUAL;
}

template <typename M>
ComparisonResult
MLPathLength<M>::
compare_multiplex(
    const MLPathLength& other
) const
{
    bool canBeDominated = true;
    bool canDominate = true;

    if (mnet != other.mnet)
    {
        throw core::OperationNotSupportedException("Cannot compare distances on different networks");
    }

    long num_intralayer_steps1 = 0;
    long num_intralayer_steps2 = 0;

    for (auto layer1: *mnet->layers())
    {
        long l1 = length(layer1,layer1);
        num_intralayer_steps1 += l1;
        long l2 = other.length(layer1,layer1);
        num_intralayer_steps2 += l2;

        if (l1 > l2)
        {
            canDominate = false;
        }

        else if (l1 < l2)
        {
            canBeDominated = false;
        }

        if (!canBeDominated && !canDominate)
        {
            return ComparisonResult::INCOMPARABLE;
        }
    }

    if (canDominate && !canBeDominated)
    {
        return ComparisonResult::GREATER_THAN;
    }

    if (canBeDominated && !canDominate)
    {
        return ComparisonResult::LESS_THAN;
    }

    //if (canDominate && canBeDominated)
    return ComparisonResult::EQUAL;
}

template <typename M>
ComparisonResult
MLPathLength<M>::
compare_simple(
    const MLPathLength& other
) const
{
    bool canBeDominated = true;
    bool canDominate = true;

    if (mnet != other.mnet)
    {
        throw core::OperationNotSupportedException("Cannot compare distances on different networks");
    }

    long l1 = length();
    long l2 = other.length();

    if (l1 > l2)
    {
        canDominate = false;
    }

    else if (l1 < l2)
    {
        canBeDominated = false;
    }

    if (!canBeDominated && !canDominate)
    {
        return ComparisonResult::INCOMPARABLE;
    }

    if (canDominate && !canBeDominated)
    {
        return ComparisonResult::GREATER_THAN;
    }

    if (canBeDominated && !canDominate)
    {
        return ComparisonResult::LESS_THAN;
    }

    //if (canDominate && canBeDominated)
    return ComparisonResult::EQUAL;
}

template <typename M>
bool
MLPathLength<M>::
operator<(
    const MLPathLength& other
) const
{
    //return total_length < other.total_length;
    return ts < other.ts;
}

template <typename M>
bool
MLPathLength<M>::
operator>(
    const MLPathLength& other
) const
{
    //return total_length > other.total_length;
    return ts > other.ts;
}

template <typename M>
bool
MLPathLength<M>::
operator==(
    const MLPathLength& other
) const
{
    //return total_length == other.total_length;
    return ts == other.ts;
}

template <typename M>
bool
MLPathLength<M>::
operator!=(
    const MLPathLength& other
) const
{
    //return total_length != other.total_length;
    return ts != other.ts;
}

template <typename M>
std::string
MLPathLength<M>::
to_string(
) const
{
    std::string res;

    for (auto layer: *mnet->layers())
    {
        long l = length(layer,layer);
        res += " " + layer->name + ":" + std::to_string(l);
    }

    for (auto layer1: *mnet->layers())
    {
        for (auto layer2: *mnet->layers())
        {
            if (layer1==layer2)
            {
                continue;
            }

            long l = length(layer1,layer2);

            if (l!=0)
            {
                res += " " + layer1->name + "->" + layer2->name + ":" + std::to_string(l);
            }
        }
    }

    return res;
}

}
}

