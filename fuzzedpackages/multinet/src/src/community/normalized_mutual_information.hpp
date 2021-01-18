/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_NET_COMMUNITY_NORMALIZEDMUTUALINFORMATION_H_
#define UU_NET_COMMUNITY_NORMALIZEDMUTUALINFORMATION_H_

#undef CS

namespace uu {
namespace net {


/**
 * Only for partitioning commiunity structures.
 * @param com1 a community structure
 * @param com1 a community structure
 * @param n number of vertices
 * @return the normalized mutual information between the two community structures
 * (1 if they are equal, down to 0)
 */
template <typename CS>
double
normalized_mutual_information(
    const CS& com1,
    const CS& com2,
    int n
)
{
    double entropy_c1 = 0;

    for (auto community: *com1)
    {
        int size1 = community->size();

        if (size1 == 0)
        {
            continue;
        }

        entropy_c1 -= (double)size1/n * std::log2((double)size1/n);
    }

    double entropy_c2 = 0;

    for (auto community: *com2)
    {
        int size2 = community->size();

        if (size2 == 0)
        {
            continue;
        }

        entropy_c2 -= (double)size2/n * std::log2((double)size2/n);
    }

    double info = 0;

    for (auto community1: *com1)
    {
        for (auto community2: *com2)
        {
            // intersection - @todo call some existing function...
            size_t common_nodes = 0;

            for (auto v: *community1)
            {
                if (community2->contains(v))
                {
                    common_nodes++;
                }
            }

            int size1 = community1->size();
            int size2 = community2->size();

            if (size1==0 || size2==0 || common_nodes==0)
            {
                continue;
            }

            info += (double)common_nodes/n * std::log2((double)n*common_nodes/(size1*size2));
        }
    }

    return info/((entropy_c1+entropy_c2)/2);
}

}
}

#endif
