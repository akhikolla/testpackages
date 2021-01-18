#ifndef UU_MEASURES_DISTANCE_H_
#define UU_MEASURES_DISTANCE_H_

#include <unordered_map>
#include "core/exceptions/assert_not_null.hpp"
#include "core/utils/Counter.hpp"
#include "objects/Vertex.hpp"
#include "objects/EdgeMode.hpp"
#include "objects/MLPathLength.hpp"
#include "networks/_impl/containers/GenericObjectList.hpp"
#include "measures/neighborhood.hpp"

namespace uu {
namespace net {


template <typename M>
std::unordered_map<const Vertex*, std::set<MLPathLength<M>> >
        pareto_distance(
            const M* net,
            const Vertex* from
        );



template <typename M>
std::unordered_map<const Vertex*, std::set<MLPathLength<M>> >
        pareto_distance(
            const M* mnet,
            const Vertex* from
        )
{
    class TimestampComparator
    {
      public:
        int
        operator()(
            const std::pair<MLPathLength<M>,long>& lhs,
            const std::pair<MLPathLength<M>,long>& rhs
        )
        {
            return lhs.second < rhs.second;

        }
    };
    std::unordered_map<const Vertex*,std::set<std::pair<MLPathLength<M>,long>,TimestampComparator> > distances;
    // timestamps, used for efficiency reasons to avoid processing edges when no changes have occurred since the last iteration
    long ts = 0;
    std::unordered_map<const typename M::layer_type*, core::PairCounter<const Vertex*, const Vertex*>> last_updated;

    // initialize distance array - for every target vertex there is still no found path leading to it...
    for (auto actor: *mnet->actors())
    {
        distances[actor] = std::set<std::pair<MLPathLength<M>,long>,TimestampComparator>();
    } // ...except for the source node, reachable from itself via an empty path

    MLPathLength<M> empty(mnet);
    distances[from].insert(std::pair<MLPathLength<M>,long>(empty,ts));

    bool changes; // keep updating the paths until when no changes occur during one full scan of the edges

    do
    {
        changes = false;

        for (auto layer: *mnet->layers())
        {
            for (auto node_from: *layer->vertices())
            {
                for (auto node_to: *layer->edges()->neighbors(node_from,EdgeMode::OUT))
                {
                    ts++;
                    // last updated
                    long lastUpdate;

                    if (last_updated[layer].count(node_from,node_to)==0)
                    {
                        lastUpdate = -1;
                    }

                    else
                    {
                        lastUpdate = last_updated[layer].count(node_from,node_to);
                    }

                    last_updated[layer].set(node_from,node_to,ts);

                    //cout << ts << " " << node_from->actor->name << " on " << node_from->layer->name <<  " -> " << node_to->actor->name << " on " << node_to->layer->name << endl;

                    // if no tmp shortest paths exist to this, do nothing and continue
                    //if (distances[actor1].empty()) {
                    //    //cout << "No paths to " << mnet.getGlobalName(actor1) << endl;
                    //    continue;
                    //}
                    // expand each temporary distance to e.v1 and see if it generates a new shortest path to e.v2
                    for (auto dist: distances[node_from])
                    {
                        ts++;

                        if (dist.second < lastUpdate)
                        {
                            //cout << "Already processed: " << lastUpdate << ": "<< dist.second << endl;
                            // distance already processed
                            continue;// TODO for efficiency: paths are sorted e.v1 most recently updated, so we do not need to examine the others
                        }

                        // otherwise, extend the distance to reach e.v2
                        // TOADD: check it's not a cycle, for efficiency reasons (?)
                        // Extend
                        MLPathLength<M> extended_distance = dist.first;
                        extended_distance.ts = ts;
                        extended_distance.step(layer, layer);
                        //cout << "producing candidate: " << extended_distance << endl;

                        // compare the new distance with the other temporary distances to e.v2
                        bool should_be_inserted = true;
                        std::set<std::pair<MLPathLength<M>,long>,TimestampComparator> dominated; // here we store the distances that will be removed if dominated by the new one

                        for (auto previous: distances[node_to])
                        {
                            // check dominance, i.e., if this is a shorter distance
                            //cout << "comparison " << extended_distance << " vs. " << previous << ": ";

                            ComparisonResult dominance = extended_distance.compare(previous.first,ComparisonType::FULL);

                            switch (dominance)
                            {
                            case ComparisonResult::LESS_THAN: // stop here
                                should_be_inserted = false;
                                //cout << "dominated" << endl;
                                break;

                            case ComparisonResult::EQUAL:
                                // this means that the number of steps in each layer is the same.
                                // Only one of them will be kept
                                should_be_inserted = false;
                                //cout << "equal" << endl;
                                // go on with the others
                                break;

                            case ComparisonResult::INCOMPARABLE: // incomparable
                                //cout << "inc." << endl;
                                // go on with the others
                                break;

                            case ComparisonResult::GREATER_THAN: // dominates -> insert it in the list of paths to be removed
                                dominated.insert(previous);
                                //cout << " dominated - remove " << endl;
                                //debug("   - REMOV " + currentPath);
                                break;
                            }
                        }

                        if (should_be_inserted)
                        {
                            //cout << " INSERT NEW for " << actor2->name << " - " << extended_distance << endl;
                            distances[node_to].insert(std::pair<MLPathLength<M>,long>(extended_distance,ts));
                            //cout << "insert " << mnet.getGlobalName(actor2) << " - " << extended_distance << endl;
                            //cout << "add " << paths[toGlobalId].size() << "\n";
                            //cout << "New path " << fromGlobalId << " => "
                            //        << toGlobalId << extended_path << "\n";
                            changes = true;
                        }

                        // remove dominated paths
                        // ?? why not just remove?
                        std::set<std::pair<MLPathLength<M>,long>,TimestampComparator> diff;
                        std::set_difference(distances[node_to].begin(),
                                            distances[node_to].end(), dominated.begin(),
                                            dominated.end(), std::inserter(diff, diff.end()));
                        distances[node_to] = diff;
                    }
                }
            }
        }
    }
    while (changes);

    std::unordered_map<const Vertex*, std::set<MLPathLength<M>> > result;

    for (auto p: distances)
    {
        for (auto dist: p.second)
        {
            result[p.first].insert(dist.first);
            //cout << "new dist to " <<  mnet->get_actor(p.first)->name << ": " << dist.first.to_string() << " (" << result[mnet->get_actor(p.first)].size() << ")" << endl;
        }
    }

    return result;
    //cout << "here?\n";
}
}
}


#endif
