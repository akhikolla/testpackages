namespace uu {
namespace net {

template <typename G>
std::unique_ptr<CommunityStructure<Community<const Vertex*>>>
label_propagation(
    const G* net
)
{

    //NodeListSharedPtr nodes = mnet->get_nodes(layer); === net->vertices()

    std::unordered_map<const Vertex*, size_t> membership; // community membership
    std::vector<const Vertex*> order; // order vector to decide in which order to process the nodes
    // Initialize labels
    int l=0;

    for (auto node: *net->vertices())
    {
        membership[node] = l;
        order.push_back(node);
        l++;
    }

    while (true)
    {

        // Compute order of node processing
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

        //std::shuffle(order.begin(), order.end(), std::default_random_engine(seed));
        // NOT AVAILABLE IN ALL COMPILERS
        core::shuffle(order.begin(), order.end(), std::default_random_engine(seed));


        // re-assign labels
        for (auto node: order)
        {
            core::Counter<size_t> num_neighbors;

            for (auto neigh: *net->edges()->neighbors(node,EdgeMode::INOUT))
            {
                // TODO make it also for directed graphs?
                num_neighbors.inc(membership.at(neigh));
            }

            membership[node] = num_neighbors.max();
        }

        /* Check stopping condition */
        for (auto node: order)
        {

            core::Counter<int> num_neighbors;

            for (auto neigh: *net->edges()->neighbors(node,EdgeMode::INOUT))
            {
                num_neighbors.inc(membership.at(neigh));
            }

            if (num_neighbors.count(membership.at(node)) != num_neighbors.count(num_neighbors.max()))
            {
                continue;
            }
        }

        break;
    }

    // Build community structure
    return to_community_structure(membership);

}

}
}
