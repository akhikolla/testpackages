/*
 * VarAttrib.h
 *
 *  Created on: Jun 9, 2011
 *      Author: ianfellows
 */

#ifndef VARATTRIBH_
#define VARATTRIBH_

#include <vector>
#include <string>
#include <Rcpp.h>

namespace lolog{

/*!
 * Vertex variable attributed
 */
class VarAttrib {
public:
    enum VarClass {DOUBLE,INTEGER, CATEGORICAL};
    VarAttrib(){
        type = VarAttrib::DOUBLE;
        name = "";
    }
    virtual ~VarAttrib(){}
    virtual bool isDouble(){
        return type==VarAttrib::DOUBLE;
    }
    virtual bool isCategorical(){
        return type==VarAttrib::INTEGER;
    }
    virtual bool isInteger(){
        return type==VarAttrib::CATEGORICAL;
    }
    virtual std::string getName(){
        return name;
    }
    virtual void setName(std::string newName){
        name = newName;
    }
protected:
    VarClass type;
    std::string name;
};

class DiscreteAttrib : public VarAttrib{
protected:
    std::vector<std::string> labs;
    bool hasLb,hasUb;
    int lb,ub;
public:
    DiscreteAttrib(){
        type = VarAttrib::INTEGER;
        hasUb=false;
        hasLb=false;
        ub=0;
        lb=0;
        name="";
    }
    virtual ~DiscreteAttrib(){}
    virtual void setLabels(std::vector<std::string> l){
        labs = l;
    }
    virtual const std::vector<std::string>& labels() const{
        return labs;
    }
    virtual bool hasLowerBound(){
        return hasLb;
    }
    virtual bool hasUpperBound(){
        return hasUb;
    }
    virtual int lowerBound(){
        return lb;
    }
    virtual int upperBound(){
        return ub;
    }
    virtual void setLowerBound(int lower){
        try{
            if(hasUb && ub<lower)
                throw std::range_error("lower bound can not be set to be larger than upper bound");
            hasLb=true;
            lb=lower;
        } catch( std::exception& _Ex__ ) {
            forward_exception_to_r( _Ex__ );
        }
    }
    virtual void setUpperBound(int upper){
        try{
            if(hasLb && upper<lb)
                std::range_error("upper bound can not be set to be larger than lower bound");
            hasUb=true;
            ub=upper;
        } catch( std::exception& _Ex__ ) {
            forward_exception_to_r( _Ex__ );
        }
    }
    virtual void removeBound(bool upper){
        if(upper){
            hasUb=false;
        }else{
            hasLb=false;
        }
    }
};

class ContinAttrib : public VarAttrib{
protected:
    bool hasLb,hasUb;
    double lb,ub;
public:
    ContinAttrib(){
        type = VarAttrib::DOUBLE;
        hasUb=false;
        hasLb=false;
        ub=0;
        lb=0;
        name="";
    }
    virtual ~ContinAttrib(){}
    virtual bool hasLowerBound(){
        return hasLb;
    }
    virtual bool hasUpperBound(){
        return hasUb;
    }
    virtual double lowerBound(){
        return lb;
    }
    virtual double upperBound(){
        return ub;
    }
    virtual void setLowerBound(double lower){
        try{
            if(hasUb && ub<lower)
                throw std::range_error("lower bound can not be set to be larger than upper bound");
            hasLb=true;
            lb=lower;
        } catch( std::exception& _Ex__ ) {
            forward_exception_to_r( _Ex__ );
        }
    }
    virtual void setUpperBound(double upper){
        try{
            if(hasLb && upper<lb)
                std::range_error("upper bound can not be set to be larger than lower bound");
            hasUb=true;
            ub=upper;
        } catch( std::exception& _Ex__ ) {
            forward_exception_to_r( _Ex__ );
        }
    }
    virtual void removeBound(bool upper){
        if(upper){
            hasUb=false;
        }else{
            hasLb=false;
        }
    }
};

}

#endif /* VARATTRIBH_ */
