#ifndef TESTS_H_
#define TESTS_H_
#include <assert.h>
#include <Rcpp.h>
#include <Rcpp/iostream/Rstreambuf.h>
#include <string>
#include "util.h"

namespace lolog{
namespace tests{


extern std::string testContext;
extern std::map< std::string, void(*)()> testFunctions;

void addTestFunction(std::string name, void(*test)() );

void registerLologTests();


/*
 * Some (rather simple) unit test macros named for consistency with the testthat package.
 * Would it be better to use a full-on framework?
 */
#define EXPECT_TRUE(x) try{ \
		if(!(x)){ \
			Rcpp::Rcerr << "Test Failed. (" << #x <<") : line " \
				<< __LINE__ << " of file " << __FILE__ <<"\n"; \
			::Rf_error("failed");\
		}\
	}catch(int e){ \
		Rcpp::Rcerr << "Test produced error. (" << #x <<") : line " \
		<< __LINE__ << " of file " << __FILE__ <<"\n"; \
		::Rf_error("failed");\
	}

#define EXPECT_EQUAL(x,y) try{ \
		if(!((x)==(y))){ \
			Rcpp::Rcerr << "Test Failed. (" << #x <<") : line " \
				<< __LINE__ << " of file " << __FILE__ <<"\n"; \
			::Rf_error("failed");\
		}\
	}catch(int e){ \
		Rcpp::Rcerr << "Test produced error. (" << #x <<") : line " \
		<< __LINE__ << " of file " << __FILE__ <<"\n"; \
		::Rf_error("failed");\
	}

#define EXPECT_NEAR(x,y) try{ \
		if(!near((x),(y))){ \
			Rcpp::Rcerr << "Test Failed. (" << #x <<") : line " \
				<< __LINE__ << " of file " << __FILE__ <<"\n" << (x) << "!=" << (y)<<" "; \
			::Rf_error("failed");\
		}\
	}catch(int e){ \
		Rcpp::Rcerr << "Test produced error. (" << #x <<") : line " \
		<< __LINE__ << " of file " << __FILE__ <<"\n" << (x) << "!=" << (y); \
		::Rf_error("failed");\
	}

#ifdef INSIDE
#define RUN_TEST(x) try{ \
		Rcpp::Rcout << testContext <<" : "<< #x << " : "; \
		x;\
		Rcpp::Rcout << "Complete" << "\n\t"; \
	}catch(int e){ \
		Rcpp::Rcout << "Test produced error. (" << #x <<") : line " \
		<< __LINE__ << " of file " << __FILE__ <<"\n"; \
		::Rf_error("failed");\
	}
#endif

#ifndef INSIDE
#define RUN_TEST(x) try{                                      \
x;                                                            \
}catch(int e){                                                \
  Rcpp::Rcout << "Test produced error. (" << #x <<") : line " \
              << __LINE__ << " of file " << __FILE__ <<"\n";  \
  ::Rf_error("failed");                                       \
}
#endif

/*!
 * Runs all tests
 */
void runLologTests();


}
}
#endif /* TESTS_H_ */
