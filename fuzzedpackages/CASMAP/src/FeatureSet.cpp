/*
 * FeatureSet.cpp
 *
 *
 *  Created on: Sep 8, 2016
 *      Author: mabaker
 */
//TODO: get rid of defaultfloat

#include <iostream> // hexfloat, defaultfloat
#include <sstream>

#include "FeatureSet.h"
#include "Exception.h"
#include <vector>
#include <iterator>

namespace SignificantPattern
{

    FeatureSet::FeatureSet () : alphaVector(), pValueVector() {}
    FeatureSet::~FeatureSet () {}

    const std::string FeatureSet::COL_SEP = "\t";
    //const std::string FeatureSet::HEADER_PROPS = "p-value" + FeatureSet::COL_SEP + "a";
    //const std::string FeatureSet::HEADER_PROPS = "p-value" + FeatureSet::COL_SEP + "score" + FeatureSet::COL_SEP + "OR";
    const std::string FeatureSet::HEADER_PROPS = "p-value";

    void FeatureSet::addFeatureProps(longint alpha, double pValue)
    {
        alphaVector.push_back(alpha);
        pValueVector.push_back(pValue);
    }

    std::string const& FeatureSet::getHeaderProps() const {
        return HEADER_PROPS;
    }
    std::string const FeatureSet::getLineProps(size_t i) const {
        std::stringstream ss;
        ss.setf(std::ios_base::scientific, std::ios_base::floatfield);
        ss << pValueVector[i];

        return ss.str();
    }

    void FeatureSet::writeHeaderToFile(std::ofstream& file) const
    {
        file << getHeaderProps() << COL_SEP << getHeaderFeature() << std::endl;
    }
    void FeatureSet::writeLineToFile(std::ofstream& file, size_t i) const
    {
        file << getLineProps(i) << COL_SEP << getLineFeature(i) << std::endl;
    }

    void FeatureSet::writeToFile(const std::string& filename) const
    {
        std::ofstream file;
        file.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
        try
        {
            file.open(filename.c_str());
        }
        catch (const std::ios_base::failure& e)
        {
            throw Exception("Failed opening " + filename + " for writing");
        }
        writeHeaderToFile(file);
        for (size_t i=0; i<getLength(); ++i) writeLineToFile(file, i);
        file.close();

    }





    ItemsetSet::ItemsetSet () : FeatureSet(), itemsetsVector() {}
    ItemsetSet::~ItemsetSet () {}

    const std::string ItemsetSet::ITEMS_SEP = ";"; //LP: output
    //const std::string ItemsetSet::HEADER_FEATURE = "itemset";
    const std::string ItemsetSet::HEADER_FEATURE = "index_1;...;index_N";

    void ItemsetSet::addFeature(const std::vector<longint> itemset, longint alpha, double pValue)
    {
        addFeatureProps(alpha, pValue);
        itemsetsVector.push_back(itemset);
    }

    std::string const& ItemsetSet::getHeaderFeature() const {
        return HEADER_FEATURE;
    }
    std::string const ItemsetSet::getLineFeature(size_t i) const {
        std::stringstream ss;
        std::vector<longint> itemset = itemsetsVector[i];
        size_t n = itemset.size();
        for(size_t i = 0; i < n-1; i++)
            ss << itemset[i] << ITEMS_SEP;
        ss << itemset[n-1];
        return ss.str();
    }


    ItemsetSetWithOddsRatio::ItemsetSetWithOddsRatio() : ItemsetSet(), oddsRatioVector(), scoreVector() {}
    ItemsetSetWithOddsRatio::~ItemsetSetWithOddsRatio() {}

    const std::string ItemsetSetWithOddsRatio::HEADER_PROPS_WITH_ODDS_RATIO = ItemsetSetWithOddsRatio::HEADER_PROPS + ItemsetSetWithOddsRatio::COL_SEP + "score" + ItemsetSetWithOddsRatio::COL_SEP + "OR";

    void ItemsetSetWithOddsRatio::addFeature(const std::vector<longint> itemset,
                                              longint alpha, double pValue) {
        addFeature(itemset, alpha, -1, -1, pValue);
    }
    void ItemsetSetWithOddsRatio::addFeature(const std::vector<longint> itemset,
                                              longint alpha, double score, double oddsRatio,
                                              double pValue) {
        super::addFeature(itemset, alpha, pValue);
        scoreVector.push_back(score);
        oddsRatioVector.push_back(oddsRatio);
    }

    std::string const& ItemsetSetWithOddsRatio::getHeaderProps() const {
        return HEADER_PROPS_WITH_ODDS_RATIO;
    }
    std::string const ItemsetSetWithOddsRatio::getLineProps(size_t i) const {
        std::stringstream ss;
        ss << super::getLineProps(i) << COL_SEP;
        ss.unsetf(std::ios_base::floatfield);
        ss << scoreVector[i] << COL_SEP
           << oddsRatioVector[i];
        return ss.str();
    }



    IntervalSet::IntervalSet () : FeatureSet(), startVector(), endVector() {}
    IntervalSet::~IntervalSet () {}

    //const std::string IntervalSet::HEADER_FEATURE = "start" + IntervalSet::COL_SEP + "end";
    const std::string IntervalSet::HEADER_FEATURE = "index_1;...;index_N";

    void IntervalSet::addFeature(longint start, longint end, longint alpha, double pValue)
    {
        addFeatureProps(alpha, pValue);
        startVector.push_back(start);
        endVector.push_back(end);
    }

    void IntervalSet::getLAndTauVectors(std::vector<longint>& lVector, std::vector<longint>& tauVector) const
    {
        for (size_t i=0; i<getLength(); ++i)
        {
            longint tau = startVector[i];
            longint l = endVector[i] - tau + 1;
            lVector.push_back(l);
            tauVector.push_back(tau);
        }
    }

    std::string const& IntervalSet::getHeaderFeature() const {
        return HEADER_FEATURE;
    }

    std::string const IntervalSet::getLineFeature(size_t i) const {
    	//LP: output
		//Compute 'index_1;index_2;index_3;...;index_N'
		std::vector<int> myVec;
		std::ostringstream oss;

		for(int j=startVector[i]; j<= endVector[i]; j++){
			myVec.push_back(j);
		}
		if (!myVec.empty())
		  {
			// Convert all but the last element to avoid a trailing ";"
			std::copy(myVec.begin(), myVec.end()-1, std::ostream_iterator<int>(oss, ";"));

			// Now add the last element with no delimiter
			oss << myVec.back();
		  }

        std::stringstream ss;
        ss << oss.str();
        //ss << startVector[i] << COL_SEP << endVector[i];
        return ss.str();
    }



    IntervalSetWithFreq::IntervalSetWithFreq() : IntervalSet(), xVector() {}
    IntervalSetWithFreq::~IntervalSetWithFreq() {}

    const std::string IntervalSetWithFreq::HEADER_PROPS_WITH_FREQ = IntervalSetWithFreq::HEADER_PROPS + IntervalSetWithFreq::COL_SEP + "x";

    void IntervalSetWithFreq::addFeature(longint start, longint end,
                                        longint alpha, double pValue) {
        addFeature(start, end, alpha, -1, pValue);
    }
    void IntervalSetWithFreq::addFeature(longint start, longint end,
                                        longint alpha, longint x,
                                        double pValue) {
        super::addFeature(start, end, alpha, pValue);
        xVector.push_back(x);
    }

    std::string const& IntervalSetWithFreq::getHeaderProps() const {
        return HEADER_PROPS_WITH_FREQ;
    }
    std::string const IntervalSetWithFreq::getLineProps(size_t i) const {
        std::stringstream ss;
        ss << super::getLineProps(i) << COL_SEP << xVector[i];
        return ss.str();
    }

    IntervalSetWithOddsRatio::IntervalSetWithOddsRatio() : IntervalSet(), oddsRatioVector(), scoreVector() {}
    IntervalSetWithOddsRatio::~IntervalSetWithOddsRatio() {}

    const std::string IntervalSetWithOddsRatio::HEADER_PROPS_WITH_ODDS_RATIO = IntervalSetWithOddsRatio::HEADER_PROPS + IntervalSetWithOddsRatio::COL_SEP + "score" + IntervalSetWithOddsRatio::COL_SEP + "OR";

    void IntervalSetWithOddsRatio::addFeature(longint start, longint end,
                                         longint alpha, double pValue) {
        addFeature(start, end, alpha, -1, -1, pValue);
    }
    void IntervalSetWithOddsRatio::addFeature(longint start, longint end,
                                         longint alpha, double score, double oddsRatio,
                                         double pValue) {
        super::addFeature(start, end, alpha, pValue);
        scoreVector.push_back(score);
        oddsRatioVector.push_back(oddsRatio);
    }

    std::string const& IntervalSetWithOddsRatio::getHeaderProps() const {
        return HEADER_PROPS_WITH_ODDS_RATIO;
    }
    std::string const IntervalSetWithOddsRatio::getLineProps(size_t i) const {
        std::stringstream ss;
        ss << super::getLineProps(i) << COL_SEP;
        ss.unsetf(std::ios_base::floatfield);

        ss << scoreVector[i] << COL_SEP
           << oddsRatioVector[i];
        return ss.str();
    }

} /* namespace SignificantPattern */
