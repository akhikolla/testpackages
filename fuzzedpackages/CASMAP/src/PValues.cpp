/*
 * PValues.cpp
 *
 *  Created on: Sep 8, 2016
 *      Author: mabaker
 */

#include <fstream>

#include "PValues.h"
#include "Exception.h"

using namespace std;

namespace SignificantPattern
{

    PValues::PValues ()
    : pValueVector()
    {
    }

    PValues::~PValues ()
    {
    }

    void PValues::addPValue(double pValue)
    {
        pValueVector.push_back(pValue);
    }

    void PValues::writeToFile(const std::string& filename)
    {
        ofstream file;
        file.exceptions ( ifstream::failbit | ifstream::badbit );
        try
        {
            file.open(filename.c_str());
        }
        catch (const std::ios_base::failure& e)
        {
            throw Exception("Failed opening " + filename + " for writing");
        }
        for (size_t i=0; i<pValueVector.size(); ++i)
            file << pValueVector[i] << endl;
        file.close();

    }

} /* namespace SignificantPattern */
