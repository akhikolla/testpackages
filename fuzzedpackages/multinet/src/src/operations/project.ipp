namespace uu {
namespace net {

template <class M, class Target>
void
project_unweighted(
    const M* net,
    const typename M::layer_type* from,
    const typename M::layer_type* to,
    Target* target
)
{
    for (auto v: *to->vertices())
    {
        target->vertices()->add(v);
    }

    for (auto v: *from->vertices())
    {
        for (auto u1: *net->interlayer_edges()->neighbors(from, to, v, EdgeMode::INOUT))
        {
            for (auto u2: *net->interlayer_edges()->neighbors(from, to, v, EdgeMode::INOUT))
            {
                if (u1 <= u2)
                {
                    continue;
                }

                target->edges()->add(u1, u2);
            }
        }
    }

}



template <class M>
void
project_temporal(
    const M* net,
    const typename M::layer_type* from,
    const typename M::layer_type* to,
    TemporalNetwork* target,
    size_t delta_time
)
{
    for (auto v: *to->vertices())
    {
        target->vertices()->add(v);
    }

    for (auto v: *from->vertices())
    {
        for (auto e1: *net->interlayer_edges()->incident(from, to, v, EdgeMode::INOUT))
        {
            auto t1 = net->interlayer_edges()->attr()->get_time(e1, "t").value;
            // @todo check time

            for (auto e2: *net->interlayer_edges()->incident(from, to, v, EdgeMode::INOUT))
            {
                if (e1 <= e2)
                {
                    continue;
                }

                auto t2 = net->interlayer_edges()->attr()->get_time(e2, "t").value;
                // @todo check time

                auto time_diff = std::chrono::seconds(std::max(t1, t2) - std::min(t1, t2)).count();

                if (delta_time <= time_diff)
                {
                    continue;
                }

                auto u1 = (e1->v1 != v) ? e1->v1 : e1->v2;
                auto u2 = (e2->v1 != v) ? e2->v1 : e2->v2;

                if (u1 == u2)
                {
                    continue;
                }

                auto edge = target->edges()->add(u1, u2);
                target->set_time(edge, std::min(t1, t2));
            }
        }
    }

}

}
}
