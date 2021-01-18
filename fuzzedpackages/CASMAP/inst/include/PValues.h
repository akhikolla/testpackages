/*
 * PValues.h
 *
 *  Created on: Sep 8, 2016
 *      Author: mabaker
 */

#ifndef SRC_PVALUES_H
#define SRC_PVALUES_H

#include <vector>
#include <string>

namespace SignificantPattern
{

    class PValues
    {
    private:
        std::vector<double> pValueVector;
    public:
        PValues ();
        virtual ~PValues ();
        std::vector<double>& getPValueVector() { return pValueVector; }

        void addPValue(double pValue);
        void writeToFile(const std::string& filename);
    };

} /* namespace SignificantPattern */

#endif /* SRC_PVALUES_H */
