/*
 * FeatureSet.h
 *
 *  Created on: Sep 8, 2016
 *      Author: mabaker
 */

#ifndef SRC_FEATURESET_H_
#define SRC_FEATURESET_H_

#include <fstream>
#include <string>
#include <vector>

#include "types.h"

namespace SignificantPattern
{

    /**
     * Base class of a container of significant features of a binary sequence
     * and its properties (p-value, or cell count).
     */
    class FeatureSet
    {
    private:
        /**
         * positive labels (cell) count
         */
        std::vector<longint> alphaVector;
        /**
         * p-values
         */
        std::vector<double> pValueVector;

    protected:
        static const std::string HEADER_PROPS;

        void addFeatureProps(longint alpha, double pValue);

        virtual std::string const& getHeaderProps() const;
        virtual std::string const getLineProps(size_t i) const;
        virtual std::string const& getHeaderFeature() const = 0;
        virtual std::string const getLineFeature(size_t i) const = 0;
        virtual void writeHeaderToFile(std::ofstream& file) const;
        virtual void writeLineToFile(std::ofstream& file, size_t i) const;

    public:
       FeatureSet ();
        virtual ~FeatureSet ();

        /**
         * columns separator
         */
        static const std::string COL_SEP;

        inline std::vector<longint> const& getAlphaVector() const { return alphaVector; }
        inline std::vector<double> const& getPValueVector() const { return pValueVector; }

        inline size_t getLength() const { return pValueVector.size(); }
        virtual void writeToFile(const std::string& filename) const;
    };

    /**
     * Itemsets (sets of positions) of binary sequence as features.
     */
    class ItemsetSet :  public FeatureSet
    {
    private:
        /**
         * super class pattern for code independence of changes in inheritance
         */
        typedef FeatureSet super;

        /**
         * itemsets
         */
        std::vector< std::vector<longint> > itemsetsVector;

    protected:
        static const std::string HEADER_FEATURE;

        virtual std::string const& getHeaderFeature() const override;
        virtual std::string const getLineFeature(size_t i) const override;

    public:
        ItemsetSet ();
        virtual ~ItemsetSet ();

        /**
         * items separator
         */
        static const std::string ITEMS_SEP;

        inline std::vector<std::vector<longint>> const &getItemsetsVector() const { return itemsetsVector; }

        virtual void addFeature(const std::vector<longint> itemset, longint alpha, double pValue);
    };

    /**
     * Intervals (sets of start and end pairs) of binary sequence as features,
     * with their score and odds ratio.
     */
    class ItemsetSetWithOddsRatio : public ItemsetSet
    {
    private:
        /**
         * super class pattern for code independence of changes in inheritance
         */
        typedef ItemsetSet super;

        /**
         * odds ratio of intervals
         */
        std::vector<double> oddsRatioVector;

        /**
         * scores of intervals
         */
        std::vector<double> scoreVector;

    protected:
        static const std::string HEADER_PROPS_WITH_ODDS_RATIO;

        virtual std::string const& getHeaderProps() const override;
        virtual std::string const getLineProps(size_t i) const override;

    public:
        ItemsetSetWithOddsRatio();
        virtual ~ItemsetSetWithOddsRatio();

        inline std::vector<double> const &getOddsRatioVector() const { return oddsRatioVector; }
        inline std::vector<double> const &getScoreVector() const { return scoreVector; }
        virtual void addFeature(const std::vector<longint> itemset, longint alpha,
                                double pValue) override;
        void addFeature(const std::vector<longint> itemset, longint alpha, double score,
                        double oddsRatio, double pValue);
    };

    /**
     * Intervals (sets of start and end pairs) of binary sequence as features.
     */
    class IntervalSet :  public FeatureSet
    {
    private:
        /**
         * super class pattern for code independence of changes in inheritance
         */
        typedef FeatureSet super;

        /**
         * start positions of intervals
         */
        std::vector<longint> startVector;
        /**
         * end positions of an intervals
         */
        std::vector<longint> endVector;

    protected:
        static const std::string HEADER_FEATURE;

        virtual std::string const& getHeaderFeature() const override;
        virtual std::string const getLineFeature(size_t i) const override;

    public:
        IntervalSet ();
        virtual ~IntervalSet ();
        inline std::vector<longint> const& getStartVector() const { return startVector; }
        inline std::vector<longint> const& getEndVector() const { return endVector; }

        void getLAndTauVectors(std::vector<longint>& lVector, std::vector<longint>& tauVector) const;

        virtual void addFeature(longint start, longint end, longint alpha, double pValue);
    };

    /**
     * Intervals (sets of start and end pairs) of binary sequence as features,
     * with their frequencies counts.
     */
    class IntervalSetWithFreq : public IntervalSet
    {
    private:
        /**
         * super class pattern for code independence of changes in inheritance
         */
        typedef IntervalSet super;

        /**
         * frequencies of intervals
         */
        std::vector<longint> xVector;

    protected:
        static const std::string HEADER_PROPS_WITH_FREQ;

        virtual std::string const& getHeaderProps() const override;
        virtual std::string const getLineProps(size_t i) const override;

    public:
        IntervalSetWithFreq();
        virtual ~IntervalSetWithFreq();

        inline std::vector<longint> const &getXVector() const { return xVector; }
        virtual void addFeature(longint start, longint end, longint alpha,
                                double pValue) override;
        void addFeature(longint start, longint end, longint alpha, longint x,
                        double pValue);
    };

    /**
     * Intervals (sets of start and end pairs) of binary sequence as features,
     * with their score and odds ratio.
     */
    class IntervalSetWithOddsRatio : public IntervalSet
    {
    private:
        /**
         * super class pattern for code independence of changes in inheritance
         */
        typedef IntervalSet super;

        /**
         * odds ratio of intervals
         */
        std::vector<double> oddsRatioVector;

        /**
         * scores of intervals
         */
        std::vector<double> scoreVector;

    protected:
        static const std::string HEADER_PROPS_WITH_ODDS_RATIO;

        virtual std::string const& getHeaderProps() const override;
        virtual std::string const getLineProps(size_t i) const override;

    public:
        IntervalSetWithOddsRatio();
        virtual ~IntervalSetWithOddsRatio();

        inline std::vector<double> const &getOddsRatioVector() const { return oddsRatioVector; }
        inline std::vector<double> const &getScoreVector() const { return scoreVector; }
        virtual void addFeature(longint start, longint end, longint alpha,
                                double pValue) override;
        void addFeature(longint start, longint end, longint alpha, double score,
                        double oddsRatio, double pValue);
    };

} /* namespace SignificantPattern */

#endif /* SRC_FEATURESET_H_ */
