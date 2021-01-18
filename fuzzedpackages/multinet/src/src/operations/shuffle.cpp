#include "operations/shuffle.hpp"

#include "networks/Network.hpp"
#include "objects/EdgeDir.hpp"
#include "core/attributes/conversion.hpp"
#include "core/attributes/Value.hpp"
#include "core/attributes/Time.hpp"
#include "core/exceptions/OperationNotSupportedException.hpp"
#include <vector>
#include <fstream>
#include <algorithm>

namespace uu {
namespace net {


void
shuffle(
    uu::net::OrderedMultiplexNetwork* net,
    size_t num
)
{
    for (auto&& l: *net->layers())
    {
        if (l->edges()->size() < 2)
        {
            continue;
        }

        auto link_swap_succes = false;

        for (size_t i = 0; i < num; i += 1)
        {
            size_t trials = 0;

            do
            {
                auto edge1 = l->edges()->get_at_random();
                auto edge2 = l->edges()->get_at_random();

                bool same_edge = (edge1 == edge2);

                //checks that the same edge has not been picked twice
                while (same_edge)
                {
                    edge1 = l->edges()->get_at_random();

                    same_edge = (edge1 == edge2);
                }

                // decide how the links should be swapped

                auto coin_flip = uu::core::getRandomInt(2);

                link_swap_succes = false;

                auto edge1_v1 = edge1->v1;
                auto edge1_v2 = edge1->v2;

                auto edge2_v1 = edge2->v1;
                auto edge2_v2 = edge2->v2;


                // checks for loops
                if (coin_flip == 0 && (edge1_v1 != edge2_v2 && edge1_v2 != edge2_v1))
                {

                    if (l->edges()->get(edge1_v1, edge2_v2) || l->edges()->get(edge2_v1, edge1_v2))
                    {
                        continue;
                    }

                    l->edges()->erase(edge1);
                    l->edges()->erase(edge2);

                    l->edges()->add(edge1_v1, edge2_v2);
                    l->edges()->add(edge2_v1, edge1_v2);

                    link_swap_succes = true;
                }

                //checks for loops
                else if (edge1_v1 != edge2_v1 && edge1_v2 != edge2_v2)
                {

                    if (l->edges()->get(edge1_v1, edge2_v1) || l->edges()->get(edge1_v2, edge2_v2))
                    {
                        continue;
                    }

                    l->edges()->erase(edge1);
                    l->edges()->erase(edge2);

                    l->edges()->add(edge1_v1, edge2_v1);
                    l->edges()->add(edge1_v2, edge2_v2);

                    link_swap_succes = true;
                }

            }
            while (!link_swap_succes && trials++<10);
        }
    }
}

}
}

