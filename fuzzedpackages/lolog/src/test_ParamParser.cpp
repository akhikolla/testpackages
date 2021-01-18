
#include <tests.h>
#include <test_ParamParser.h>
#include <ParamParser.h>

namespace lolog{
namespace tests{

void testParsing(){
    Rcpp:List params = Rcpp::List::create(1, Rcpp::Named("a")=2, Rcpp::Named("b") = "ss");
    CharacterVector v = params.names();
    std::string val;
    val = v.at(0);
    ParamParser p  = ParamParser("test", params);

    int first = p.parseNext("ll", 3);
    EXPECT_TRUE(first == 1);

    int second = p.parseNext("a", 1);
    EXPECT_TRUE(second == 2);


    std::string third = p.parseNext<std::string>("b");
    EXPECT_TRUE(third == "ss");

    std::string fourth = p.parseNext("other", "default");
    EXPECT_TRUE(fourth == "default");

    p.end();
}


void testParamParser(){
    RUN_TEST(testParsing());
}

}
}
