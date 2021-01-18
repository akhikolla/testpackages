/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_MNET_DATASTRUCTURES_GRAPHS_MULTILAYERNETWORKTYPE_H_
#define UU_MNET_DATASTRUCTURES_GRAPHS_MULTILAYERNETWORKTYPE_H_

namespace uu {
namespace net {

struct TMultilayerNetworkType
{
    bool is_vertex_aligned = false;
    bool is_ordered = false;
    bool is_directed = false;
    //bool is_weighted = false;
    //bool is_probabilistic = false;
    //bool is_temporal = false;
    //bool is_attributed = false;
    //bool allows_multi_edges = false;
    //bool allows_loops = false;

};

}
}

#endif
