//  From c++ Cookbook by D. Ryan Stephens (Author) , Christopher Diggins (Author) , Jonathan Turkanis (Author) 

#include "TypeConversion.h"
//constructor of type conversion

typeConversion::typeConversion(const int initialIntiger, string &initialString)
{

    Intiger = initialIntiger;
    String = initialString;
    stringstream ss;
    ss.str("");
    ss << Intiger;

    initialString = ss.str();

    ss.clear();

}
typeConversion::typeConversion(const double initialIntiger, string &initialString)
{

    Double = initialIntiger;
    String = initialString;
    stringstream ss;
    ss.str("");
    ss << Double;
    initialString = ss.str();


    ss.clear();

}
typeConversion::typeConversion(const unsigned int initialIntiger,string &initialString)
{
    Intiger = initialIntiger;
    String = initialString;
    stringstream ss;
    ss.str("");
    ss << Intiger;

    initialString = ss.str();

    ss.clear();


}

typeConversion::~typeConversion()
{

}
