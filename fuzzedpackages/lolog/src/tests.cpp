#include <tests.h>
#include <Rcpp/iostream/Rstreambuf.h>
#include <test_BinaryNet.h>
#include <test_Constraints.h>
#include <test_LatentOrderLikelihood.h>
#include <test_Stats.h>
#include <test_ParamParser.h>

namespace lolog{
namespace tests{

std::string testContext;
std::map< std::string, void(*)()> testFunctions;

void addTestFunction(std::string name, void(*test)() ){
    testFunctions.insert(std::make_pair(name, test));
}

void registerLologTests(){
    addTestFunction("testBinaryNet", testBinaryNet);
    addTestFunction("testStats", testStats);
    addTestFunction("testConstraints", testConstraints);
    addTestFunction("testLatent", testLatent);
    addTestFunction("testParamParser", testParamParser);
}


void runLologTests(){
#ifdef INSIDE
    Rcpp::Rcout << "\n\t";
#endif
    registerLologTests();
    //testBinaryNet();
    //testStats();
    //testConstraints();
    //testLatent();
    std::map< std::string, void(*)()>::iterator it;
    for ( it = testFunctions.begin(); it != testFunctions.end(); it++ ){
        testContext = it->first;
        it->second();
    }
#ifdef INSIDE
    Rcpp::Rcout << "All C++ Tests Complete\n";
#endif
}


}
}
