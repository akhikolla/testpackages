#include <test_Constraints.h>

#include <Rcpp.h>

#include <BinaryNet.h>
#include <Stat.h>
#include <Stats.h>
#include <Offsets.h>
#include <Constraint.h>
#include <Constraints.h>
#include <Model.h>
#include <VarAttrib.h>
#include <tests.h>
namespace lolog {

namespace tests{


template<class Engine>
void testBoundedDegree(){
    GetRNGstate();
    IntegerMatrix tmp(0,2);
    BinaryNet<Engine> net(tmp,30);
    boost::shared_ptr< AbstractStat<Engine> > ed(new Stat<Engine, Edges<Engine> >());
    Rcpp::List ll;
    ll.push_back(2);
    ll.push_back(10);
    boost::shared_ptr< AbstractOffset<Engine> > off(
            new Constraint<Engine,BoundedDegree<Engine> >(ll));
    Model<Engine> model(net);
    model.addStatPtr(ed);
    model.addOffsetPtr(off);
    model.calculate();
    model.setThetas(std::vector<double>(1,0));

    EXPECT_TRUE(model.offset().at(0) < -100000)

    PutRNGstate();
}


void testConstraints(){
    RUN_TEST(testBoundedDegree<Undirected>());
}

}

} /* namespace lolog */
