#ifndef TRACE_H
#define TRACE_H


#include "../mixtureDefinition.h"


#include <iostream>
#include <vector>
#include <cctype>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

class Trace {
	private:

		unsigned numCodonSpecificParamTypes;

		std::vector<std::vector<double>> stdDevSynthesisRateTrace; //mixture, samples
		std::vector<double> stdDevSynthesisRateAcceptanceRateTrace; //samples TODO: Correctly sized for the time being,
        //however, it will need to be changed at some point when there are some adjustments to hyper parameter acceptance/rejection
		std::vector<std::vector<std::vector<double>>>synthesisRateAcceptanceRateTrace; //order: expressionCategory, gene, sample
		std::vector<std::vector<double>> codonSpecificAcceptanceRateTrace;//order: codon, sample
                std::vector<std::vector<double>> nseSpecificAcceptanceRateTrace;//order: codon, sample
		std::vector<std::vector<std::vector<float>>> synthesisRateTrace;//order: expression category, gene, samples
		std::vector<std::vector<unsigned>> mixtureAssignmentTrace;//order: numGenes, samples
		std::vector<std::vector<double>> mixtureProbabilitiesTrace;//order: numMixtures, samples
		std::vector<std::vector<std::vector<std::vector<float>>>> codonSpecificParameterTrace; //order: paramType, category, numParam, samples
		//std::vector<std::vector<std::vector<double>>> codonSpecificParameterTraceTwo; //order: category, numParam, samples
		std::vector<mixtureDefinition> *categories;


		//ROC Trace:
		std::vector<std::vector <double>> synthesisOffsetTrace;
		std::vector<std::vector <double>> synthesisOffsetAcceptanceRateTrace;
		std::vector<std::vector <double>> observedSynthesisNoiseTrace;


		//FONSE Trace:
                std::vector<double> initiationCostTrace;
                std::vector<double> initiationCostAcceptanceRateTrace;

		//PANSE Trace:
		std::vector<std::vector <double>> partitionFunctionTrace;
		std::vector<double> partitionFunctionTraceAcceptanceRateTrace;

		//--------------------------------------//
		//------ Initialization Functions ------//
		//--------------------------------------//
		void initializeSharedTraces(unsigned samples, unsigned num_genes, unsigned numSelectionCategories, unsigned numMixtures,
			std::vector<mixtureDefinition> &_categories, unsigned maxGrouping,std::vector<double> init_phi, std::vector<unsigned> init_mix_assign, unsigned numObservedPhiSets,bool estimateSynthesisRate = true);

		void initStdDevSynthesisRateTrace(unsigned numSelectionCategories, unsigned samples);
		void initSynthesisRateAcceptanceRateTrace(unsigned num_genes, unsigned numExpressionCategories);
		void initSynthesisRateTrace(unsigned samples, unsigned num_genes, unsigned numExpressionCategories,std::vector<double> init_phi,bool estimateSynthesisRate = true);
		void initMixtureAssignmentTrace(unsigned samples, unsigned num_genes,std::vector<unsigned> init_mix_assign);
		void initMixtureProbabilitiesTrace(unsigned samples, unsigned numMixtures);
		void initCodonSpecificParameterTrace(unsigned samples, unsigned numMutationCategories, unsigned numParam, unsigned paramType);


		//ROC Specific:
		void initSynthesisOffsetTrace(unsigned samples, unsigned numPhiGroupings);
		void initObservedSynthesisNoiseTrace(unsigned samples, unsigned numPhiGroupings);


		//FONSE Specific:
                void initInitiationCostTrace(unsigned samples);

		//PANSE Specific:
		void initPartitionFunctionTrace(unsigned samples, unsigned numPartitionFunctionsGroupings);

	public:
		//Constructors & Destructors:
		Trace();
		virtual ~Trace();
		Trace(unsigned _numCodonSpecificParamTypes);


		//Initialization Functions:
		void initializePATrace(unsigned samples, unsigned num_genes, unsigned numAlphaCategories,
			unsigned numLambdaPrimeCategories, unsigned numParam, unsigned numMixtures,
			std::vector<mixtureDefinition> &_categories, unsigned maxGrouping,unsigned numObservedPhiSets, std::vector<double> init_phi,
                        std::vector<unsigned> init_mix_assign, bool estimateSynthesisRate = true);


		void initializeROCTrace(unsigned samples, unsigned num_genes, unsigned numMutationCategories,
			unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures, std::vector<mixtureDefinition> &_categories,
			unsigned maxGrouping, unsigned numObservedPhiSets,std::vector<double> init_phi,
                        std::vector<unsigned> init_mix_assign, bool estimateSynthesisRate = true);


		void initializeFONSETrace(unsigned samples, unsigned num_genes, unsigned numMutationCategories,
			unsigned numSelectionCategories, unsigned numParam, unsigned numMixtures,
			std::vector<mixtureDefinition> &_categories, unsigned maxGrouping,unsigned numObservedPhiSets,std::vector<double> init_phi,
                        std::vector<unsigned> init_mix_assign, bool estimateSynthesisRate = true);


		void initializePANSETrace(unsigned samples, unsigned num_genes, unsigned numAlphaCategories,
			unsigned numLambdaPrimeCategories, unsigned numParam, unsigned numMixtures,
			std::vector<mixtureDefinition> &_categories, unsigned maxGrouping,unsigned numObservedPhiSets,std::vector<double> init_phi,
                        std::vector<unsigned> init_mix_assign, bool estimateSynthesisRate = true);


		//------------------------------//
		//------ Getter Functions ------//
		//------------------------------//
        std::vector<double> getStdDevSynthesisRateTrace(unsigned selectionCategory);
        std::vector<double> getExpectedSynthesisRateTrace();
        std::vector<double> getStdDevSynthesisRateAcceptanceRateTrace();
        std::vector<std::vector<std::vector<float>>> getSynthesisRateTrace();
        std::vector<double> getSynthesisRateAcceptanceRateTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex);
        std::vector<std::vector<std::vector<double>>> getSynthesisRateAcceptanceRateTrace();
        std::vector<double> getCodonSpecificAcceptanceRateTraceForAA(std::string aa);
        std::vector<double> getCodonSpecificAcceptanceRateTraceForCodon(std::string codon);
        std::vector<float> getSynthesisRateTraceForGene(unsigned geneIndex); //will build the trace appropriately based on what cat you are in
        std::vector<float> getSynthesisRateTraceByMixtureElementForGene(unsigned mixtureElement, unsigned geneIndex);
        std::vector<unsigned> getMixtureAssignmentTraceForGene(unsigned geneIndex);
        std::vector<double> getMixtureProbabilitiesTraceForMixture(unsigned mixtureIndex);
        std::vector<std::vector<unsigned>> getMixtureAssignmentTrace();
        std::vector<std::vector<double>> getMixtureProbabilitiesTrace();
        std::vector<std::vector<double>> getCodonSpecificAcceptanceRateTrace();
        unsigned getSynthesisRateCategory(unsigned mixtureElement);
        unsigned getCodonSpecificCategory(unsigned mixtureElement, unsigned paramType);
		std::vector<std::vector<std::vector<std::vector<float>>>>* getCodonSpecificParameterTrace();



        std::vector<double> getSynthesisOffsetTrace(unsigned index);
        std::vector<double> getSynthesisOffsetAcceptanceRateTraceForIndex(unsigned index);
        std::vector<double> getObservedSynthesisNoiseTrace(unsigned index);

        //ROC Specific:
        std::vector<float> getCodonSpecificParameterTraceByMixtureElementForCodon(unsigned mixtureElement, std::string& codon,
                unsigned paramType, bool withoutReference = true);
	std::vector<float> getCodonSpecificParameterTraceByGeneElementForCodon(unsigned geneIndex, std::string& codon,
		unsigned paramType, bool withoutReference = true);
        
        std::vector<std::vector<std::vector<float>>> getCodonSpecificParameterTraceByParamType(unsigned paramType);
        std::vector<std::vector<double>> getSynthesisOffsetAcceptanceRateTrace();


        //FONSE Specific:
        std::vector<double> getInitiationCostTrace();
        std::vector<double> getInitiationCostAcceptanceRateTrace();

        //PANSE Specific:
        std::vector<double> getPartitionFunctionTrace(unsigned mixtureIndex);
        std::vector<double> getPartitionFunctionAcceptanceRateTrace();

        //------------------------------//
		//------ Update Functions ------//
		//------------------------------//
        void updateStdDevSynthesisRateTrace(unsigned sample, double stdDevSynthesisRate, unsigned synthesisRateCategory);
        void updateStdDevSynthesisRateAcceptanceRateTrace(double acceptanceLevel);
        void updateSynthesisRateAcceptanceRateTrace(unsigned category, unsigned geneIndex, double acceptanceLevel);
        void updateCodonSpecificAcceptanceRateTrace(unsigned codonIndex, double acceptanceLevel);
        void updateSynthesisRateTrace(unsigned sample, unsigned geneIndex, std::vector<std::vector <double>> &currentExpressionLevel);
        void updateMixtureAssignmentTrace(unsigned sample, unsigned geneIndex, unsigned value);
        void updateMixtureProbabilitiesTrace(unsigned samples, std::vector<double> &categoryProbabilities);


        //ROC Specific:
        void updateCodonSpecificParameterTraceForAA(unsigned sample, std::string aa, std::vector<std::vector<double>> &curParam, unsigned paramType);
        void updateSynthesisOffsetTrace(unsigned index, unsigned sample, double value);
        void updateSynthesisOffsetAcceptanceRateTrace(unsigned index, double value);
        void updateObservedSynthesisNoiseTrace(unsigned index, unsigned sample, double value);


        //FONSE Specific:
        void updateInitiationCostTrace(unsigned sample,double value);
        void updateInitiationCostAcceptanceRateTrace(double value);


        //PANSE Specific:
        void updateCodonSpecificParameterTraceForCodon(unsigned sample, std::string codon, std::vector<std::vector<double>> &curParam, unsigned paramType);
        void updateNseRateSpecificAcceptanceRateTrace(unsigned codonIndex, double acceptanceLevel);
        void updatePartitionFunctionTrace(unsigned index, unsigned sample, double value);
        void updatePartitionFunctionAcceptanceRateTrace(double value);

        void resizeNumberCodonSpecificParameterTrace(unsigned _numCodonSpecificParamTypes);

        //R Section:
#ifndef STANDALONE
        //Getter Functions:
        std::vector<double> getSynthesisRateAcceptanceRateTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex);//R WRAPPER
        std::vector<float> getSynthesisRateTraceForGeneR(unsigned geneIndex);//R WRAPPER
        std::vector<float> getSynthesisRateTraceByMixtureElementForGeneR(unsigned mixtureElement, unsigned geneIndex);//R WRAPPER
        std::vector<unsigned> getMixtureAssignmentTraceForGeneR(unsigned geneIndex);//R WRAPPER
        std::vector<double> getMixtureProbabilitiesTraceForMixtureR(unsigned mixtureIndex);//R WRAPPER
        std::vector<std::vector<double>> getStdDevSynthesisRateTraces();
        unsigned getNumberOfMixtures();
        std::vector<std::vector<double>> getNseRateSpecificAcceptanceRateTrace();




        //Setter Functions:
        void setStdDevSynthesisRateTraces(std::vector<std::vector<double>> _stdDevSynthesisRateTrace);
        void setStdDevSynthesisRateAcceptanceRateTrace(std::vector<double> _stdDevSynthesisRateAcceptanceRateTrace);
        void setSynthesisRateTrace(std::vector<std::vector<std::vector<float>>> _synthesisRateTrace);
        void setSynthesisRateAcceptanceRateTrace(std::vector<std::vector<std::vector<double>>>_synthesisRateAcceptanceRateTrace);
        void setMixtureAssignmentTrace(std::vector<std::vector<unsigned>> _mixtureAssignmentTrace);
        void setMixtureProbabilitiesTrace(std::vector<std::vector<double>> _mixtureProbabilitiesTrace);
        void setCodonSpecificAcceptanceRateTrace(std::vector<std::vector<double>> _cspAcceptanceRateTrace);
        void setCategories(std::vector<mixtureDefinition> &_categories);
        void setNseRateSpecificAcceptanceRateTrace(std::vector<std::vector<double>> _nseAcceptanceRateTrace);


        //ROC Specific:updateSynthesisOffsetAcceptanceRateTrace
	std::vector<float> getCodonSpecificParameterTraceByMixtureElementForCodonR(unsigned mixtureElement, std::string& codon, unsigned paramType,
		        bool withoutReference);
        std::vector<std::vector<double>> getSynthesisOffsetTraceR();
        std::vector<std::vector<double>> getObservedSynthesisNoiseTraceR();


        void setSynthesisOffsetTrace(std::vector<std::vector <double> > _NoiseOffsetTrace);
        void setSynthesisOffsetAcceptanceRateTrace(std::vector<std::vector <double> > _NoiseOffsetAcceptanceRateTrace);
        void setObservedSynthesisNoiseTrace(std::vector<std::vector <double> > _ObservedSynthesisNoiseTrace);
        void setCodonSpecificParameterTrace(std::vector<std::vector<std::vector<float>>> _parameterTrace, unsigned paramType);

        //FONSE specific:
        void setInitiationCostTrace(std::vector <double> _InitiationCostTrace);
        //PANSE Specific:

        std::vector<double> getPartitionFunctionTraceR(unsigned mixtureIndex);
        void setPartitionFunctionTraces(std::vector<std::vector <double> > _PartitionFunctionTrace);
        std::vector<std::vector<double>> getPartitionFunctionTraces();

        bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);

#endif
};
#endif // TRACE_H
