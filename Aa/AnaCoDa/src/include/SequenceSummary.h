#ifndef SequenceSummary_H
#define SequenceSummary_H


#include "Utility.h"


#include <string>
#include <map>
#include <algorithm>
#include <cctype>
#include <vector>
#include <array>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

class SequenceSummary
{
	private:

		std::array<unsigned, 64> ncodons; //64 for the number of codons.
		std::array<unsigned, 22> naa; //22 for the number of amino acids.
		std::vector <std::vector <unsigned>> codonPositions; // used in FONSEModel.
        // outer index is the codonID, size of 64 for number of codons
        // inner index is the position of each occurrence of the codonID specified.

        std::vector <std::vector <unsigned>> RFPCount;
		// outer index is the RFPCount for the category specified via index
		// inner index is number of position

		std::vector <std::array <unsigned, 64>> sumRFPCount;
		// outer index is the RFPCount for the category specified via index
		// inner index, 64, is the number of codons

		std::vector <unsigned> positionCodonID;
		// index is the number of position, where a value (codonID) is set


	public:

		//Static Member Variables:
		static const std::string Ser2;
		static const std::vector<std::string> AminoAcidArray;
		static const std::string codonArray[];
		static const std::string codonArrayParameter[];
		static const std::map<std::string, unsigned> aaToIndex;
		static const std::map<std::string, unsigned> codonToIndexWithReference;
		static const std::map<std::string, unsigned> codonToIndexWithoutReference;


		//Constructors & Destructors:
		explicit SequenceSummary();
		SequenceSummary(const std::string& sequence);
		SequenceSummary(const SequenceSummary& other);
		SequenceSummary& operator=(const SequenceSummary& other);
		bool operator==(const SequenceSummary& other) const;
		virtual ~SequenceSummary(); //TODO:Why is this virtual????


		//Data Manipulation Functions (All tested):
		unsigned getAACountForAA(std::string aa);
		unsigned getAACountForAA(unsigned aaIndex);
		unsigned getCodonCountForCodon(std::string& codon);
		unsigned getCodonCountForCodon(unsigned codonIndex);
		std::vector <unsigned> *getCodonPositions(std::string codon);
		std::vector <unsigned> *getCodonPositions(unsigned index);


		//RFP Functions (for PA and PANSE models) (All tested):
        //RFP Count has an inner vector indexed by position Sum RFP Count has an inner index of Codon Type
		void initRFPCount(unsigned numCategories);
		std::vector <unsigned> getRFPCount(unsigned RFPCountColumn = 0u);
		unsigned getSingleRFPCount(unsigned position, unsigned RFPCountColumn = 0u);
		void setRFPCount(std::vector <unsigned> arg, unsigned RFPCountColumn = 0u);
		unsigned getSumTotalRFPCount(unsigned RFPCountColumn = 0u);

        void initSumRFPCount(unsigned numCategories);
		std::array <unsigned, 64> getSumRFPCount(unsigned RFPCountColumn = 0u);
		void setSumRFPCount(std::array <unsigned, 64> arg, unsigned RFPCountColumn = 0u);

        //These functions deal with a single codon at a time
        unsigned getCodonSpecificSumRFPCount(std::string codon, unsigned RFPCountColumn = 0u);
        unsigned getCodonSpecificSumRFPCount(unsigned codonIndex, unsigned RFPCountColumn = 0u);
        void setCodonSpecificSumRFPCount(unsigned codonIndex, unsigned value, unsigned RFPCountColumn = 0u);

        //Poisitonal information about Codons
        std::vector <unsigned> getPositionCodonID(); //Used in PANSE for getting codon positions over gene
        void setPositionCodonID(std::vector <unsigned> arg);

		//Other Functions (All tested):
		void clear();
		bool processSequence(const std::string& sequence);
        bool processPA(std::vector <std::vector <int>> table);
        bool processPANSE(std::vector <std::vector <int>> table);

		//Static Functions:
		static unsigned AAToAAIndex(std::string aa); //Moving to CT
		static void AAIndexToCodonRange(unsigned aaIndex, unsigned& start, unsigned& end, bool forParamVector = false); //Moving to CT
		static void AAToCodonRange(std::string aa, unsigned& start, unsigned& end, bool forParamVector = false); //Moving to CT
		static std::vector<std::string> AAToCodon(std::string aa, bool forParamVector = false); //Moving to CT, but used in R currently
		static std::string codonToAA(std::string& codon); //Moving to CT
		static unsigned codonToIndex(std::string& codon, bool forParamVector = false); //Moving to CT
		static unsigned codonToAAIndex(std::string& codon); //Moving to CT
		static std::string indexToAA(unsigned aaIndex); //Moving to CT
		static std::string indexToCodon(unsigned index, bool forParamVector = false); //Moving to CT
		static unsigned GetNumCodonsForAA(std::string& aa, bool forParamVector = false); //Moving to CT
		static char complimentNucleotide(char ch); //Tested
		static std::vector<std::string> aminoAcids(); //Moving to CT, but used in R currently
		static std::vector<std::string> codons(); //Moving to CT, but used in R currently


	protected:
};

#endif // SequenceSummary_H
