#ifndef UU_GENERATION_EVOLUTIONMODEL_H_
#define UU_GENERATION_EVOLUTIONMODEL_H_


namespace uu {
namespace net {

typedef int evolution_strategy;

const int EVOLUTION_DEGREE=0;


/**********************************************************************/
/** Evolution models **************************************************/
/**********************************************************************/

template <typename M>
class
    EvolutionModel
{
  public:

    virtual
    ~EvolutionModel() = 0;

    virtual void
    init_step(
        M* mnet,
        typename M::layer_type* layer,
        GenericObjectList<Vertex>& available_actors
    ) = 0;

    virtual void
    internal_evolution_step(
        M* mnet,
        typename M::layer_type* layer,
        GenericObjectList<Vertex>& available_actors
    ) = 0;

    virtual void
    external_evolution_step(
        M* mnet,
        typename M::layer_type* target_layer,
        GenericObjectList<Vertex>& available_actors,
        const typename M::layer_type* ext_layer
    ) = 0;
};

template <typename M>
EvolutionModel<M>::
~EvolutionModel() {}






}
}

#endif
