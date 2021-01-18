#include <StatController.h>
#include <Stat.h>
#include <Stats.h>
#include <Offset.h>
#include <Offsets.h>
#include <Constraint.h>
#include <Constraints.h>
namespace lolog{



typedef boost::shared_ptr< AbstractStat<Directed> > DirStatPtr;
typedef boost::shared_ptr< std::map< std::string, DirStatPtr > > StatMapPtr;
typedef boost::shared_ptr< AbstractOffset<Directed> > DirOffsetPtr;
typedef boost::shared_ptr< std::map< std::string, DirOffsetPtr > > OffsetMapPtr;
template<> StatMapPtr StatController<Directed>::statMapPtr = StatMapPtr(new std::map< std::string, DirStatPtr >);
template<> OffsetMapPtr StatController<Directed>::offsetMapPtr =
        OffsetMapPtr(new std::map< std::string, DirOffsetPtr >);

typedef boost::shared_ptr< AbstractStat<Undirected> > UndirStatPtr;
typedef boost::shared_ptr< std::map< std::string, UndirStatPtr > > UndirStatMapPtr;
typedef boost::shared_ptr< AbstractOffset<Undirected> > UndirOffsetPtr;
typedef boost::shared_ptr< std::map< std::string, UndirOffsetPtr > > UndirOffsetMapPtr;
template<> UndirStatMapPtr StatController<Undirected>::statMapPtr =
        UndirStatMapPtr(new std::map< std::string, UndirStatPtr >);
template<> UndirOffsetMapPtr StatController<Undirected>::offsetMapPtr =
        UndirOffsetMapPtr(new std::map< std::string, UndirOffsetPtr >);


}


//[[Rcpp::export(name=".initStats")]]
void initStats(){
    using namespace lolog;
    /*
     * Directed network statistics
     */
    registerStatistic( DirStatPtr( new DirectedEdges() ) );
    registerStatistic( DirStatPtr( new DirectedTriangles() ) );
    registerStatistic( DirStatPtr( new DirectedMutual() ) );
    registerStatistic( DirStatPtr( new DirectedNodeMatch() ) );
    registerStatistic( DirStatPtr( new DirectedDegree() ) );
    registerStatistic( DirStatPtr( new DirectedStar() ) );
    registerStatistic( DirStatPtr( new DirectedNodeCov() ) );
    registerStatistic( DirStatPtr( new DirectedEdgeCovSparse() ) );
    registerStatistic( DirStatPtr( new DirectedGwesp() ) );
    registerStatistic( DirStatPtr( new DirectedGeoDist() ) );
    registerStatistic( DirStatPtr( new DirectedGwDegree() ) );
    registerStatistic( DirStatPtr( new DirectedGwdsp() ) );
    registerStatistic( DirStatPtr( new DirectedEsp() ) );
    registerStatistic( DirStatPtr( new DirectedPreferentialAttachment() ) );
    registerStatistic( DirStatPtr( new DirectedNodeFactor() ) );
    registerStatistic( DirStatPtr( new DirectedAbsDiff() ) );
    registerStatistic( DirStatPtr( new DirectedEdgeCov() ) );
    registerStatistic( DirStatPtr( new DirectedTwoPath() ) );
    ////////			Offsets				/////////

    //Make registration available outside lolog compilation unit
    R_RegisterCCallable("lolog",
            "registerDirectedStatistic",(DL_FUNC) &registerDirectedStatistic);
    R_RegisterCCallable("lolog",
            "registerDirectedOffset",(DL_FUNC) &registerDirectedOffset);

    /*
     * Undirected network statistics
     */
    registerStatistic( UndirStatPtr( new UndirectedEdges() ) );
    registerStatistic( UndirStatPtr( new UndirectedTriangles() ) );
    registerStatistic( UndirStatPtr( new UndirectedClustering() ) );
    registerStatistic( UndirStatPtr( new UndirectedTransitivity() ) );
    registerStatistic( UndirStatPtr( new UndirectedNodeLogMaxCov() ) );
    registerStatistic( UndirStatPtr( new UndirectedNodeMix() ) );
    registerStatistic( UndirStatPtr( new UndirectedDegree() ) );
    registerStatistic( UndirStatPtr( new UndirectedNodeMatch() ) );
    registerStatistic( UndirStatPtr( new UndirectedStar() ) );
    registerStatistic( UndirStatPtr( new UndirectedNodeCov() ) );
    registerStatistic( UndirStatPtr( new UndirectedEdgeCovSparse() ) );
    registerStatistic( UndirStatPtr( new UndirectedGwesp() ) );
    registerStatistic( UndirStatPtr( new UndirectedGeoDist() ) );
    registerStatistic( UndirStatPtr( new UndirectedGwDegree() ) );
    registerStatistic( UndirStatPtr( new UndirectedGwdsp() ) );
    registerStatistic( UndirStatPtr( new UndirectedEsp() ) );
    registerStatistic( UndirStatPtr( new UndirectedDegreeCrossProd() ) );
    registerStatistic( UndirStatPtr( new UndirectedPreferentialAttachment() ) );
    registerStatistic( UndirStatPtr( new UndirectedSharedNbrs() ) );
    registerStatistic( UndirStatPtr( new UndirectedNodeFactor() ) );
    registerStatistic( UndirStatPtr( new UndirectedAbsDiff() ) );
    registerStatistic( UndirStatPtr( new UndirectedEdgeCov() ) );
    registerStatistic( UndirStatPtr( new UndirectedTwoPath() ) );

    ////////			Offsets				/////////
    registerOffset( UndirOffsetPtr( new UndirectedBoundedDegreeConstraint() ) );
    //Make registration available outside lolog compilation unit
    R_RegisterCCallable("lolog",
            "registerUndirectedStatistic",(DL_FUNC) &registerUndirectedStatistic);
    R_RegisterCCallable("lolog",
            "registerUndirectedOffset",(DL_FUNC) &registerUndirectedOffset);

}


void registerDirectedStatistic(Rcpp::XPtr< lolog::AbstractStat<lolog::Directed> > ps){
    lolog::StatController<lolog::Directed>::addStat(
            boost::shared_ptr< lolog::AbstractStat<lolog::Directed> >(ps->vCloneUnsafe()));
}

void registerUndirectedStatistic(Rcpp::XPtr< lolog::AbstractStat<lolog::Undirected> > ps){
    lolog::StatController<lolog::Undirected>::addStat(
            boost::shared_ptr< lolog::AbstractStat<lolog::Undirected> >(ps->vCloneUnsafe()));
}

void registerDirectedOffset(Rcpp::XPtr< lolog::AbstractOffset<lolog::Directed> > ps){
    lolog::StatController<lolog::Directed>::addOffset(
            boost::shared_ptr< lolog::AbstractOffset<lolog::Directed> >(ps->vCloneUnsafe()));
}

void registerUndirectedOffset(Rcpp::XPtr< lolog::AbstractOffset<lolog::Undirected> > ps){
    lolog::StatController<lolog::Undirected>::addOffset(
            boost::shared_ptr< lolog::AbstractOffset<lolog::Undirected> >(ps->vCloneUnsafe()));
}






