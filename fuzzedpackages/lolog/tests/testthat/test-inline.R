context("inline tests")


#library(lolog)
library(BH)
library(testthat)



test_that("Inline", {
  skip_on_cran()
  # This creates a function in C++ to create an empty network of size n
  # and expose it to R.
  src <- "
  lolog::BinaryNet<lolog::Directed> makeEmptyNetwork(const int n){
    Rcpp::IntegerMatrix tmp(0,2);
    lolog::BinaryNet<lolog::Directed> net(tmp, n);
    return net;
  }
  "
  Rcpp::registerPlugin("lolog",inlineLologPlugin)
  emptyNetwork <- cppFunction(src,plugin="lolog")
  net <- emptyNetwork(10L)
  expect_true(all(!net[1:10,1:10]))
  
  
  # Test inline statistic registration
  src <- '
  // [[Rcpp::depends(lolog, BH)]]
  #include <lolog.h>
  using namespace lolog;

  template<class Engine>
  class Edges2 : public BaseStat<Engine>{
  public:
  Edges2(){}
  Edges2(List params){}
  
  std::string name(){ return "edges2";}
  
  std::vector<std::string> statNames(){
    std::vector<std::string> statnames(1,"edges2");
    return statnames;
  }
  
  void calculate(const BinaryNet<Engine>& net){
    std::vector<double> v(1,2.0*net.nEdges());
    this->stats=v;
    this->lastStats = std::vector<double>(1,0.0);
    this->thetas = std::vector<double>(1, 0.0);
  }
  
  void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
    BaseOffset<Engine>::resetLastStats();
    BaseOffset<Engine>::update(net.hasEdge(from,to) ? -2.0 : 2.0, 0);
  }
  
  bool isOrderIndependent(){return true;}
  
  bool isDyadIndependent(){return true;}
  };
  
  typedef Stat<Directed, Edges2<Directed> > DirectedEdges2;
  typedef Stat<Undirected, Edges2<Undirected> > UndirectedEdges2;

  // [[Rcpp::export]]
  void registerEdges2(){
	Rcpp::XPtr< AbstractStat<Undirected> > ps1(new UndirectedEdges2());
	  REGISTER_UNDIRECTED_STATISTIC(ps1);

	  Rcpp::XPtr< AbstractStat<Directed> > ps2(new DirectedEdges2());
	  REGISTER_DIRECTED_STATISTIC(ps2);
  }
  '
  sourceCpp(code=src)
  data(ukFaculty)
  f1 <- lolog(ukFaculty~edges)
  registerEdges2()
  f2 <- lolog(ukFaculty~edges2)
  expect_true(abs(f1$theta - f2$theta*2) < .000000000001)
})
