/**
 * History:
 * - 2018.03.09 file created, following a restructuring of the previous library.
 */

#ifndef UU_MNET_IO_MULTILAYERIOFILESECTION_H_
#define UU_MNET_IO_MULTILAYERIOFILESECTION_H_

namespace uu {
namespace net {

enum class MultilayerIOFileSection
{
    DEFAULT,
    VERSION,
    TYPE,
    ACTOR_ATTRIBUTES,
    VERTEX_ATTRIBUTES,
    EDGE_ATTRIBUTES,
    LAYERS,
    ACTORS,
    VERTICES,
    INTRALAYER_VERTICES,
    EDGES
};


}
}

#endif
