/*
 * Summary.h
 *
 *  Created on: Sep 9, 2016
 *      Author: mabaker
 */

#ifndef SRC_SUMMARY_H_
#define SRC_SUMMARY_H_

#include <cmath> // pow
#include <fstream> // ofstream
#include <string>

#include "types.h"

namespace SignificantPattern
{

    class Summary
    {
    private:
        longint N;
        longint n;
        longint L;
        longint m; // number of testable features
        longint numSignificantFeatures;
        longint numFeaturesProcessed;
        longint L_max; // maximum feature set size to be processed, or 0 for unlimited
        double delta; // Associated testability region
        double delta_opt; // Corrected testability region
        double alpha; // Corrected significance threshold

    protected:
        virtual const std::string getFeatureString() const = 0;
        virtual const std::string getFeaturesString() const = 0;

        void writeToFileStream(std::ofstream& file) const;
        inline virtual void writeExtrasToFileStream(std::ofstream& file) const {};

    public:
        inline Summary()
            : N(0), n(0), L(0), m(0), numSignificantFeatures(0),
              numFeaturesProcessed(0), L_max(0), delta(0), delta_opt(0),
              alpha(0) {};
        inline virtual ~Summary() {};

        inline longint getN() const { return N; }
        inline longint getn() const { return n; }
        inline longint getL() const { return L; }
        inline longint getm() const { return m; }
        inline longint getNumSignificantFeatures() const { return numSignificantFeatures; }
        inline longint getNumFeaturesProcessed() const { return numFeaturesProcessed; }
        virtual inline longint getNumFeatureSetsTotal() const { return (longint) L*((L+1)/2);}
        inline longint getL_max() const { return L_max; }
        inline double getDelta() const { return delta; }
        inline double getDelta_opt() const { return delta_opt; }
        inline double getAlpha() const { return alpha; }

        inline void setN(longint val) { N = val; };
        inline void setn(longint val) { n = val; };
        inline void setL(longint val) { L = val; };
        inline void setm(longint val) { m = val; };
        inline void setNumSignificantFeatures(longint val) { numSignificantFeatures = val; };
        inline void setNumFeaturesProcessed(longint val) { numFeaturesProcessed = val; };
        inline void setL_max(longint val) { L_max = val; };
        inline void setDelta(double val) { delta = val; };
        inline void setDelta_opt(double val) { delta_opt = val; };
        inline void setAlpha(double val) { alpha = val; };

        virtual void writeToFile(const std::string& filename) const;
    };

    class SummaryInt : public Summary
    {
    private:
        // super class pattern for code independence of changes in inheritance
        typedef Summary super;

        longint maxTestableIntervalLength;
    protected:
        virtual inline const std::string getFeatureString() const override { return "interval"; }
        virtual inline const std::string getFeaturesString() const override { return "intervals"; }

        virtual void writeExtrasToFileStream(std::ofstream& file) const override;
    public:
        inline SummaryInt()
            : maxTestableIntervalLength(0) {};
        inline virtual ~SummaryInt() {};

        inline longint getMaxTestableIntervalLength() const { return maxTestableIntervalLength; }

        inline void setMaxTestableIntervalLength(longint val) { maxTestableIntervalLength = val; };
    };

    class SummaryFais : public SummaryInt
    {
    private:
        // super class pattern for code independence of changes in inheritance
        typedef SummaryInt super;

        longint sl1, sl2, su1, su2; // Resultant testability region

    protected:
        virtual void writeExtrasToFileStream(std::ofstream& file) const override;

    public:
        inline SummaryFais()
            : sl1(0), sl2(0), su1(0), su2(0) {};
        inline virtual ~SummaryFais() {};

        inline longint getSl1() const { return sl1; }
        inline longint getSl2() const { return sl2; }
        inline longint getSu1() const { return su1; }
        inline longint getSu2() const { return su2; }

        inline void setTestabilityRegion(longint sl1_, longint sl2_, longint su1_, long su2_) {
            sl1 = sl1_; sl2 = sl2_; su1 = su1_; su2 = su2_;
        }
    };

    class SummaryWy : public SummaryFais
    {
    private:
        // super class pattern for code independence of changes in inheritance
        typedef SummaryFais super;

        double FWER;
        double FWER_opt;
        double* min_pval;
        longint J;

    protected:
        virtual void writeExtrasToFileStream(std::ofstream& file) const override;

    public:
        inline SummaryWy()
            : FWER(0), FWER_opt(0), min_pval(0), J(0) {};
        inline virtual ~SummaryWy() {};

        inline double getFWER() const { return FWER; }
        inline double getFWER_opt() const { return FWER_opt; }
        inline double* getMin_pval() const { return min_pval; }
        inline longint getJ() const { return J; }

        inline void setFWER(double val) { FWER = val; }
        inline void setFWER_opt(double val) { FWER_opt = val; }
        inline void setMin_pval(double* val) { min_pval = val; }
        inline void setJ(longint val) { J = val; }
    };

    class SummaryIset : public Summary
    {
    private:
        // super class pattern for code independence of changes in inheritance
        typedef Summary super;

    protected:
        inline const std::string getFeatureString() const override { return "itemset"; };
        inline const std::string getFeaturesString() const override { return "itemsets"; };
    };

    class SummaryFacs : public SummaryIset
    {
    private:
        // super class pattern for code independence of changes in inheritance
        typedef SummaryIset super;

        longint numItemsetsClosedProcessed;

    protected:
        void writeExtrasToFileStream(std::ofstream& file) const override;

    public:
        inline longint getNumItemsetsClosedProcessed() const {
            return numItemsetsClosedProcessed;
        }
        inline void setNumItemsetsClosedProcessed(longint val) {
            numItemsetsClosedProcessed = val;
        }

        inline longint getNumFeatureSetsTotal() const override {
            return (longint) pow(  ((long double) 2), ((long double) getL())  );
            //need to cast to long double, because pow(int, long int) not defined
        };
    };

} /* namespace SignificantPattern */

#endif /* SRC_SUMMARY_H_ */
