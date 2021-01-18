/*
 * ParamParser.h
 *
 *  Created on: Apr 4, 2018
 *      Author: ianfellows
 */

#ifndef PARAMPARSER_H_
#define PARAMPARSER_H_
#include <string>
#include <Rcpp.h>

namespace lolog {

/**
 * For use with constraints/offsets/stats. Helper in parsing parameters passed down from R in a list.
 *
 * Parsing rules:
 * 1. Unnamed passed parameters are matched positionally until a named parameter is reached.
 * 2. All names must match exactly.
 * 3. No unknown items are allowed in passed parameters.
 *
 * Usage:
 * //create a parser
 * p = ParamParser("statName", passedParams);
 *
 * //one by one pull out the stats parameter values
 * int a = p.parseNext<int>("firstParamName") // a required int parameter
 * std::string b = p.parseNext("secondParameter", "default") // a string parameter with a default value
 *
 * p.end() // throws an error if there are any passed parameters that have not been extracted.
 *
 */
class ParamParser {
    std::string name;
    Rcpp::List params;
    int nUnnamedParsed;
    int totalParsed;
protected:
    template<class T>
    T parseNext(std::string paramName, T defaultValue, bool allowDefault){
        T ret(defaultValue);
        int s = params.size();
        if(s <= nUnnamedParsed){
            if(!allowDefault)
                ::Rf_error(("Error in " + name + ": To few parameters.").c_str());
            return ret;
        }
        std::string pName;
        CharacterVector names;
        if(!::Rf_isNull(params.names())){
            names = params.names();
            pName = names.at(nUnnamedParsed);
        }else
            pName = "";
        if(pName == ""){
            try{
                ret = Rcpp::as<T>(params.at(nUnnamedParsed));
                totalParsed++;
            }catch(...){
                ::Rf_error(("Error in " + name + ": Invalid type for parameter " + paramName).c_str());
            }
            nUnnamedParsed++;
        }else{
            bool found = false;
            for(int i=nUnnamedParsed; i < s; i++){
                pName = names.at(i);
                found = pName == paramName;
                if(found){
                    try{
                        ret = Rcpp::as<T>(params.at(i));
                        totalParsed++;
                    }catch(...){
                        ::Rf_error(("Error in " + name + ": Invalid type for parameter " + paramName).c_str());
                    }
                    continue;
                }
            }
            if(!found && !allowDefault){
                ::Rf_error(("Error in " + name + ":  Required parameter " + paramName + " missing").c_str());
            }
        }
        return ret;
    }

public:
    ParamParser() : name("ParamParser"), params(), nUnnamedParsed(0), totalParsed(0){};

    ParamParser(std::string funcName) : name(funcName), params(), nUnnamedParsed(0), totalParsed(0){};

    ParamParser(std::string funcName, Rcpp::List passedParamValues) : name(funcName), params(passedParamValues), nUnnamedParsed(0), totalParsed(0){};

    virtual ~ParamParser(){};

    template<class T>
    T parseNext(std::string paramName, T defaultValue){
        return parseNext(paramName, defaultValue, true);
    }

    template<class T>
    T parseNext(std::string paramName){
        return parseNext(paramName, T(), false);
    }

    EdgeDirection parseNextDirection(std::string paramName, EdgeDirection defaultValue){
        std::string defaultDir = defaultValue == UNDIRECTED ? "undirected" : (defaultValue == IN ? "in" : "out");
        std::string par = parseNext(paramName, defaultDir, true);

        if(par == "in")
            return IN;
        if(par == "out")
            return OUT;
        if(par == "undirected")
            return UNDIRECTED;

        ::Rf_error(("Error in " + name + ":  Required parameter " + paramName + " must be 'in', 'out', or 'undirected'").c_str());
        return UNDIRECTED;
    }

    EdgeDirection parseNextDirection(std::string paramName){
        std::string par = parseNext(paramName, "", false);

        if(par == "in")
            return IN;
        if(par == "out")
            return OUT;
        if(par == "undirected")
            return UNDIRECTED;

        ::Rf_error(("Error in " + name + ":  Required parameter " + paramName + " must be 'in', 'out', or 'undirected'").c_str());
        return UNDIRECTED;
    }

    void end(){
        if(totalParsed != params.size()){
            Rf_error(("Either unknown or duplicate parameters passed to " + name).c_str());
        }
    }
};

}

#endif /* PARAMPARSER_H_ */

