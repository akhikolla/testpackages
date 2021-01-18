#include <Rcpp.h>
#include <BinaryNet.h>
#include <Stat.h>
#include <Stats.h>
#include <Constraint.h>
#include <Model.h>
#include <VarAttrib.h>
#include <tests.h>

namespace lolog{

namespace tests{



template<class Engine>
void changeStatTest(std::string statName){
    using namespace std;
    using namespace lolog;
    vector<int> vals(30,1);
    vals[3] = 3;
    vals[2] = 2;
    vals[4] = 2;
    vector<string> labels(3,"a");
    labels[1] = "b";
    labels[2] = "c";
    DiscreteAttrib attr;
    attr.setName("fact");
    attr.setLabels(labels);

    vector<int> vals1(30,1);
    vals1[3] = 2;
    vals1[2] = 2;
    vals1[4] = 2;
    vector<string> labels1(2,"a");
    labels1[1] = "b";
    DiscreteAttrib attr1;
    attr1.setName("out");
    attr1.setLabels(labels1);

    IntegerMatrix tmp(0,2);
    BinaryNet<Engine> net(tmp,30);
    GetRNGstate();
    Language call1("set.seed",wrap(1.0));
    call1.eval();
    for(int i=0;i<30;i++){
        pair<int,int> dyad = net.randomDyad();
        net.addEdge(dyad.first,dyad.second);
    }
    net.addDiscreteVariable(vals,attr);
    net.addDiscreteVariable(vals1,attr1);

    vector<double>vals2;
    for(int i=0; i<30; i++)
        vals2.push_back(Rf_runif(-90,90));
    ContinAttrib attr2;
    attr2.setName("contin");
    attr2.setLowerBound(-90);
    attr2.setUpperBound(90);
    net.addContinVariable(vals2,attr2);

    vals2.clear();
    for(int i=0; i<30; i++)
        vals2.push_back(Rf_runif(-180,180));
    ContinAttrib attr3;
    attr3.setName("contin1");
    attr3.setLowerBound(-180);
    attr3.setUpperBound(180);
    net.addContinVariable(vals2,attr3);

    //calculate network statistics
    vector<int> q;
    q.push_back(3);
    vector<int> degrees;
    degrees.push_back(2);
    degrees.push_back(3);
    degrees.push_back(14);
    Rcpp::List deg;
    deg.push_back(degrees);
    boost::shared_ptr< Stat<Engine,Edges<Engine> > > ed(new Stat<Engine,Edges<Engine> >());
    boost::shared_ptr< AbstractStat<Engine> > stat;
    Rcpp::List fact;
    fact.push_back("fact");
    //string fact = "fact";
    Rcpp::List log;
    log.push_back("out");
    log.push_back("fact");
    Rcpp::List hom;
    hom.push_back("fact");
    hom.push_back(0+UNDIRECTED);
    hom.push_back(true);
    hom.push_back(false);
    if(statName == "NodeMatch")
        stat = boost::shared_ptr< Stat<Engine,NodeMatch<Engine> > >(
                new Stat<Engine, NodeMatch<Engine> >(fact));
    else if(statName == "Triangles")
        stat = boost::shared_ptr< Stat<Engine, Triangles<Engine> > >(
                new Stat<Engine, Triangles<Engine> >());
    else if(statName == "Clustering")
        stat = boost::shared_ptr< Stat<Engine, Clustering<Engine> > >(
                new Stat<Engine, Clustering<Engine> >());
    else if(statName == "Transitivity")
        stat = boost::shared_ptr< Stat<Engine, Transitivity<Engine> > >(
                new Stat<Engine, Transitivity<Engine> >());
    else if(statName == "Degree")
        stat = boost::shared_ptr< Stat<Engine, Degree<Engine> > >(
                new Stat<Engine, Degree<Engine> >(deg));
    else if(statName == "DegreeCrossProd")
        stat = boost::shared_ptr< Stat<Engine, DegreeCrossProd<Engine> > >(
                new Stat<Engine, DegreeCrossProd<Engine> >());
    else if(statName == "Star"){
        vector<int> tmp;
        tmp.push_back(2);
        tmp.push_back(3);
        Rcpp::List l;
        l.push_back(tmp);
        l.push_back("in");
        stat = boost::shared_ptr< Stat<Engine, Star<Engine> > >(
                new Stat<Engine, Star<Engine> >(l));
    }else if(statName == "NodeCov"){
        Rcpp::List ncpar;
        ncpar.push_back("contin");
        stat = boost::shared_ptr< Stat<Engine, NodeCov<Engine> > >(
                new Stat<Engine, NodeCov<Engine> >(ncpar));
    }else if(statName == "NodeCov (discrete)"){
        Rcpp::List ncpar;
        ncpar.push_back("fact");
        stat = boost::shared_ptr< Stat<Engine, NodeCov<Engine> > >(
                new Stat<Engine, NodeCov<Engine> >(ncpar));
    }else if(statName == "NodeFactor"){
        Rcpp::List ncpar;
        ncpar.push_back("fact");
        stat = boost::shared_ptr< Stat<Engine, NodeFactor<Engine> > >(
                new Stat<Engine, NodeFactor<Engine> >(ncpar));
    }else if(statName == "Gwesp"){
        Rcpp::List ncpar;
        ncpar.push_back(.5);
        stat = boost::shared_ptr< Stat<Engine, Gwesp<Engine> > >(
                new Stat<Engine, Gwesp<Engine> >(ncpar));
    }else if(statName == "GeoDist"){
        Rcpp::List ncpar;
        ncpar.push_back("contin1");
        ncpar.push_back("contin");
        stat = boost::shared_ptr< Stat<Engine, GeoDist<Engine> > >(
                new Stat<Engine, GeoDist<Engine> >(ncpar));
    }else if(statName == "Gwdegree"){
        Rcpp::List ncpar;
        ncpar.push_back(.5);
        stat = boost::shared_ptr< Stat<Engine, GwDegree<Engine> > >(
                new Stat<Engine, GwDegree<Engine> >(ncpar));
    }else if(statName == "Gwdsp"){
        Rcpp::List ncpar;
        ncpar.push_back(.5);
        stat = boost::shared_ptr< Stat<Engine, Gwdsp<Engine> > >(
                new Stat<Engine, Gwdsp<Engine> >(ncpar));
    }else if(statName == "Esp"){
        Rcpp::List ncpar;
        ncpar.push_back(1);
        stat = boost::shared_ptr< Stat<Engine, Esp<Engine> > >(
                new Stat<Engine, Esp<Engine> >(ncpar));
    }else if(statName == "TwoPath"){
      stat = boost::shared_ptr< Stat<Engine, TwoPath<Engine> > >(
        new Stat<Engine, TwoPath<Engine> >());
    }else
        Rf_error("changeStatTest: unknown stat");
    vector<int> tmp1;
    tmp1.push_back(4);
    tmp1.push_back(1);
    tmp1.push_back(0);

    //create model
    Model<Engine> model(net);
    model.addStatPtr(ed);
    model.addStatPtr(stat);
    //model.addOffsetPtr(off);
    model.calculate();

    stat->vTheta().at(0) = 0.0;
    //Language call3("print",wrap(model.statistics()));
    //call3.eval();

    vector<int> togVars(1,0);
    togVars.push_back(1);

    vector<int> order(30,1);
    for(int i=0;i<30;i++) order[i] = i;
    for(int i=0;i<30;i++){
        pair<int,int> dyad = net.randomDyad();
        model.dyadUpdate(dyad.first,dyad.second, order, dyad.first);
        net.toggle(dyad.first,dyad.second);
    }

    for(int i=0;i<300;i++){
        pair<int,int> dyad = net.randomDyad();
        model.dyadUpdate(dyad.first,dyad.second, order, dyad.first);
        if(Rf_runif(0.0,1.0) < .5)
            model.rollback();
        else
            net.toggle(dyad.first,dyad.second);
    }

    vector<double> mcmcStats = model.statistics();
    model.calculateStatistics();
    vector<double> realStats = model.statistics();

    for(int i=0;i<realStats.size();i++){
        //cout << i << " " << mcmcStats.at(i) << " " << realStats.at(i) << " ";
        EXPECT_NEAR((mcmcStats.at(i) + .0001)/(realStats.at(i) + .0001),1.0);
    }
    //Language call4("print",wrap(mh.generateSampleStatistics(100,100,100)));
    //  call4.eval();
    PutRNGstate();
}


void testStats(){

    RUN_TEST(changeStatTest<Directed>("NodeMatch"));
    RUN_TEST(changeStatTest<Directed>("Degree"));
    RUN_TEST(changeStatTest<Directed>("Star"));
    RUN_TEST(changeStatTest<Directed>("NodeCov"));
    RUN_TEST(changeStatTest<Directed>("NodeCov (discrete)"));
    RUN_TEST(changeStatTest<Directed>("Gwesp"));
    RUN_TEST(changeStatTest<Directed>("Gwdegree"));
    RUN_TEST(changeStatTest<Directed>("Triangles"));
    RUN_TEST(changeStatTest<Directed>("Esp"));
    RUN_TEST(changeStatTest<Directed>("NodeFactor"));
    RUN_TEST(changeStatTest<Directed>("TwoPath"));

    RUN_TEST(changeStatTest<Undirected>("Triangles"));
    RUN_TEST(changeStatTest<Undirected>("Clustering"));
    RUN_TEST(changeStatTest<Undirected>("Transitivity"));
    RUN_TEST(changeStatTest<Undirected>("NodeMatch"));
    RUN_TEST(changeStatTest<Undirected>("Degree"));
    RUN_TEST(changeStatTest<Undirected>("Star"));
    RUN_TEST(changeStatTest<Undirected>("NodeCov"));
    RUN_TEST(changeStatTest<Undirected>("NodeCov (discrete)"));
    RUN_TEST(changeStatTest<Undirected>("Gwesp"));
    RUN_TEST(changeStatTest<Undirected>("GeoDist"));
    RUN_TEST(changeStatTest<Undirected>("Gwdegree"));
    RUN_TEST(changeStatTest<Undirected>("Gwdsp"));
    RUN_TEST(changeStatTest<Undirected>("Esp"));
    RUN_TEST(changeStatTest<Undirected>("DegreeCrossProd"));
    RUN_TEST(changeStatTest<Undirected>("NodeFactor"));
    RUN_TEST(changeStatTest<Undirected>("TwoPath"));

}


}
}
