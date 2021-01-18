/*
 * Summary.cpp
 *
 *  Created on: Sep 9, 2016
 *      Author: mabaker
 */

#include "Summary.h"
#include "Exception.h"

using namespace std;

namespace SignificantPattern
{

/* -------------------------------------------------------------------------- */
/* ----------------------------- Summary ------------------------------------ */
/* -------------------------------------------------------------------------- */
    void Summary::writeToFileStream(std::ofstream& summary_file) const
    {
        summary_file << "DATASET CHARACTERISTICS:" << endl;
        summary_file << "\tN = " << N << ", n = " << n << ", L = " << L << endl;
        summary_file << "RESULTS:" <<  endl;
        double percent_of_total = ((double) (100*numFeaturesProcessed))/((double) getNumFeatureSetsTotal());
        if(percent_of_total < 1e-12) // no abs() for when underflow happened
            summary_file << "Number of " << getFeaturesString() << " processed: " << numFeaturesProcessed << " (less than 1e-12% of total)." << endl;
        else
            summary_file << "Number of " << getFeaturesString() << " processed: " << numFeaturesProcessed << " (" << percent_of_total << "% of total)." << endl;
        if(L_max==0)
            summary_file << "Maximum "  << getFeatureString() << " length to be processed: unlimited" << endl;
        else
            summary_file << "Maximum " <<  getFeatureString() << " length to be processed: " << L_max << endl;
        summary_file << "Associated testability threshold: " << std::scientific << delta << endl;
        summary_file << "Number of testable "<<  getFeaturesString() << ": " << m << endl;
        summary_file << "Corrected significance threshold at level " << std::scientific << alpha << ": " << std::scientific << delta_opt << endl;
        summary_file << "Number of significantly associated " <<  getFeaturesString() << " found: " << numSignificantFeatures << endl;
        writeExtrasToFileStream(summary_file);
    }

    void Summary::writeToFile (const std::string& filename) const
    {
        //string filename = basefilename + "_summary.txt";
        ofstream summary_file;
        summary_file.exceptions ( ifstream::failbit | ifstream::badbit );
        try
        {
            summary_file.open(filename.c_str());
        }
        catch (const std::ios_base::failure& e)
        {
            throw Exception("Failed opening " + filename);
        }

        writeToFileStream(summary_file);

        summary_file.close();

    }



/* -------------------------------------------------------------------------- */
/* --------------------------- SummaryInt ---------------------------------- */
/* -------------------------------------------------------------------------- */
    void SummaryInt::writeExtrasToFileStream(std::ofstream& summary_file) const
    {
        super::writeExtrasToFileStream(summary_file);
        summary_file << "Maximum testable "  << getFeatureString() << " length: " << maxTestableIntervalLength << endl;
    }



/* -------------------------------------------------------------------------- */
/* --------------------------- SummaryFais ---------------------------------- */
/* -------------------------------------------------------------------------- */
    void SummaryFais::writeExtrasToFileStream(std::ofstream& summary_file) const
    {
        super::writeExtrasToFileStream(summary_file);
        summary_file << "Resultant testability region: L [" << sl1 << "," << sl2 << "] U [" << su1 << "," << su2 << "]" << endl;
    }



/* -------------------------------------------------------------------------- */
/* ---------------------------- SummaryWy ----------------------------------- */
/* -------------------------------------------------------------------------- */
    void SummaryWy::writeExtrasToFileStream(std::ofstream& summary_file) const
    {
        super::writeExtrasToFileStream(summary_file);
        summary_file << "FWER at testability threshold: " << std::scientific << FWER << endl;
        summary_file << "FWER at corrected significance threshold: " << std::scientific << FWER_opt << endl;
        //summary_file << "MINIMUM P-VALS (" << J << " PERMUTATIONS)" << endl;
        //for (int j=0; j<J-1; ++j)
        //	summary_file << min_pval[j] << ',';
        //summary_file << endl;
    }



/* -------------------------------------------------------------------------- */
/* --------------------------- SummaryFacs ---------------------------------- */
/* -------------------------------------------------------------------------- */
    void  SummaryFacs::writeExtrasToFileStream(std::ofstream& summary_file) const
    {
        super::writeExtrasToFileStream(summary_file);
        summary_file << "Number of closed " << getFeaturesString() << " processed: " << getNumItemsetsClosedProcessed() << endl;
    }

} /* namespace SignificantPattern */
