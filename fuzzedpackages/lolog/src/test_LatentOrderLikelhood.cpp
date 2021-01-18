#include <Rcpp.h>

#include <BinaryNet.h>
#include <Stat.h>
#include <Stats.h>
#include <Constraint.h>
#include <Model.h>
#include <VarAttrib.h>
#include <LatentOrderLikelihood.h>
#include <Ranker.h>
#include <test_LatentOrderLikelihood.h>
#include <tests.h>
namespace lolog {
namespace tests {

template<class Engine>
void lt() {
    using namespace std;
    using namespace lolog;
    vector<int> vals(30, 1);
    vals[3] = 3;
    vals[2] = 2;
    vals[4] = 2;
    vector<string> labels(3, "a");
    labels[1] = "b";
    labels[2] = "c";
    DiscreteAttrib attr;
    attr.setName("fact");
    attr.setLabels(labels);

    vector<int> vals1(30, 1);
    vals1[3] = 2;
    vals1[2] = 2;
    vals1[4] = 2;
    vector<string> labels1(2, "a");
    labels1[1] = "b";
    DiscreteAttrib attr1;
    attr1.setName("out");
    attr1.setLabels(labels1);

    IntegerMatrix tmp(0, 2);
    BinaryNet<Engine> net(tmp, 30);
    GetRNGstate();
    Language call1("set.seed", wrap(1.0));
    call1.eval();
    for (int i = 0; i < 30; i++) {
        pair<int, int> dyad = net.randomDyad();
        net.addEdge(dyad.first, dyad.second);
    }
    net.addEdge(1, 2);
    net.addEdge(3, 2);
    net.addEdge(1, 3);
    net.addDiscreteVariable(vals, attr);
    net.addDiscreteVariable(vals1, attr1);

    vector<double> vals2;
    for (int i = 0; i < 30; i++)
        vals2.push_back(Rf_runif(-90, 90));
    ContinAttrib attr2;
    attr2.setName("contin");
    attr2.setLowerBound(-90);
    attr2.setUpperBound(90);
    net.addContinVariable(vals2, attr2);

    vals2.clear();
    for (int i = 0; i < 30; i++)
        vals2.push_back(Rf_runif(-180, 180));
    ContinAttrib attr3;
    attr3.setName("contin1");
    attr3.setLowerBound(-180);
    attr3.setUpperBound(180);
    net.addContinVariable(vals2, attr3);

    //calculate network statistics
    vector<int> q;
    q.push_back(3);
    vector<int> degrees;
    degrees.push_back(2);
    degrees.push_back(3);
    degrees.push_back(14);
    Rcpp::List deg;
    deg.push_back(degrees);
    boost::shared_ptr<Stat<Engine, Edges<Engine> > > ed(
            new Stat<Engine, Edges<Engine> >());
    boost::shared_ptr<AbstractStat<Engine> > stat;
    Rcpp::List fact;
    fact.push_back("fact");
    //string fact = "fact";
    Rcpp::List log;
    log.push_back("out");
    log.push_back("fact");
    Rcpp::List hom;
    hom.push_back("fact");
    hom.push_back(0 + UNDIRECTED);
    hom.push_back(true);
    hom.push_back(false);
    stat = boost::shared_ptr<Stat<Engine, Triangles<Engine> > >(
            new Stat<Engine, Triangles<Engine> >());

    //create model
    Model<Engine> model(net);
    model.addStatPtr(ed);
    model.addStatPtr(stat);
    //model.addOffsetPtr(off);

    model.calculate();

    double llik = model.logLik();
    //Language call2("print",wrap(model.terms()));
    //call2.eval();
    //cout << "\nllik: " << llik << "\n";
    //cout << "\nnedges:" << net.nEdges()<<"\n";

    std::vector<int> ord(30);
    for (int i = 0; i < 30; i++) {
        ord[i] = i / 5;
    }
    model.setVertexOrderVector(ord);

    LatentOrderLikelihood<Engine> lol = LatentOrderLikelihood<Engine>(model);

    List result = lol.variationalModelFrame(1, .005);

    EXPECT_TRUE(model.getVertexOrderVector().size() == 30);
    lol.generateNetwork();

    model.setVertexOrderVector(std::vector<int>());
    EXPECT_TRUE(model.getVertexOrderVector().size() == 0);

    lol = LatentOrderLikelihood<Engine>(model);
    lol.generateNetwork();

    //Language call("print",result);
    //call.eval();
    PutRNGstate();
}

void rnker() {
    //Rcpp::Environment base_env("package:base");
    //Rcpp::Function set_seed_r = base_env["set.seed"];
    //set_seed_r(10);
    GetRNGstate();
    vector<int> vals1(5, 1);
    vals1[3] = 2;
    vals1[2] = 3;
    vals1[4] = 3;
    vector<int> ranks(5, 1);
    //std::cout << "Rank: Average\n";
    rank(vals1, ranks, "average");
    //for (uint i = 0; i < ranks.size(); ++i)
    //  std::cout << vals1[i] << " " << ranks[i] << std::endl;
    //std::cout << "Rank: Random\n";
    rank(vals1, ranks, "random");
    //for (uint i = 0; i < ranks.size(); ++i)
    //  std::cout << vals1[i] << " " << ranks[i] << std::endl;

    //std::cout << "Rank: Order\n";
    order(vals1, ranks);
    //for (uint i = 0; i < ranks.size(); ++i)
    //  std::cout << vals1[i] << " " << ranks[i] << std::endl;

    PutRNGstate();
}

void testLatent() {
    RUN_TEST(lt<Undirected>());
    RUN_TEST(lt<Directed>());
    RUN_TEST(rnker());

}

}
}

