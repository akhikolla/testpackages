#include "community/optimize_slices.hpp"

#include "networks/OrderedMultiplexNetwork.hpp"
#include "networks/Network.hpp"
#include "community/glouvain2.hpp"
#include "community/modularity.hpp"
#include "operations/shuffle.hpp"
#include "operations/slice.hpp"

namespace uu {
namespace net {

std::vector<double>
optimize_slices(
    const TemporalNetwork* original_net,
    size_t max_slices = 50)
{

    std::vector<double> res;

    //std::cerr << "testing network file: " << name << std::endl;

    for (size_t ts = 1; ts < max_slices; ++ts)
    {

        //std::cerr << "num slices: " << ts << std::endl;

        auto sliced_net = uu::net::slice_equal_time(original_net, ts);

        auto com = uu::net::glouvain2(sliced_net.get(), 1.0);

        double modul = ordered_modularity(sliced_net.get(), com.get(), 1.0);

        // get a high modularity result - try generalized_louvain 5 times and take the maximum
        for (size_t i=0; i<5; i++)
        {
            auto com_alt = uu::net::glouvain2(sliced_net.get(), 1.0);
            double modul_alt = ordered_modularity(sliced_net.get(), com.get(), 1.0);

            if (modul_alt > modul)
            {
                com = std::move(com_alt);
            }
        }

        auto sliced_rand_net = uu::net::slice_equal_time(original_net, ts);

        shuffle(sliced_rand_net.get(), original_net->edges()->size());

        auto rand_com = uu::net::glouvain2(sliced_rand_net.get(), 1.0);

        double rand_modul = ordered_modularity(sliced_rand_net.get(), com.get(), 1.0);

        // get a high modularity result - try generalized_louvain 5 times and take the maximum
        for (size_t i=0; i<5; i++)
        {
            auto rand_com_alt = uu::net::glouvain2(sliced_rand_net.get(), 1.0);
            double rand_modul_alt = ordered_modularity(sliced_rand_net.get(), rand_com.get(), 1.0);

            if (rand_modul_alt > rand_modul)
            {
                rand_com = std::move(rand_com_alt);
            }
        }

        auto modularity_before_shuffling = ordered_modularity(sliced_net.get(), com.get(), 1.0);
        auto modularity_after_shuffling = ordered_modularity(sliced_rand_net.get(), rand_com.get(), 1.0);

        auto normalized_modularity = modularity_before_shuffling - modularity_after_shuffling;

        res.push_back(normalized_modularity);
        /*
         std::cout
         << name << " "
         << ts << " "
         << com->size() << " "
         << rand_com->size() << " "
         << modularity_before_shuffling << " "
         << modularity_after_shuffling << " "
         << normalized_modularity
         << std::endl;
         */
    }

    return res;
}

}
}
