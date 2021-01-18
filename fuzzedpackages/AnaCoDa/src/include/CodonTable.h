#ifndef CodonTable_H
#define CodonTable_H

#include <string>
#include <map>
#include <algorithm>
#include <cctype>
#include <vector>
#include <array>
#include "Utility.h"

class CodonTable
{
    private:

        unsigned tableId;
        bool splitAA;
        std::vector<std::vector<unsigned>> codon_mapping; // dim: AA, codon



    public:
        static const std::string Ser2;
        static const std::string Ser1; // necessary for codon table 12
        static const std::string Thr4_1; // necessary for codon table 3
        static const std::string Thr4_2; // necessary for codon table 3
        static const std::string Leu1; // necessary for codon table 16, 22

		static const std::string AminoAcidArray[]; // dim: table, AA
		static const unsigned numCodonsPerAAForTable[25][26]; // dim: table, AA
		static const std::string codonTableDefinition[25];
		static const std::string codonArray[];
		static const std::string codonArrayParameter[];
		static const std::map<std::string, unsigned> aaToIndex;
		static const std::map<std::string, unsigned> codonToIndexWithReference;
		static const std::map<std::string, unsigned> codonToIndexWithoutReference;


        //Constructors & destructors:
        explicit CodonTable();
        CodonTable(unsigned _tableId, bool _splitAA);
        virtual ~CodonTable();
        CodonTable(const CodonTable& other);
        CodonTable& operator=(const CodonTable& other);

        void setupCodonTable();
        unsigned AAToAAIndex(std::string aa);
        unsigned getNumCodons(std::string aa);
        unsigned getNumCodons(unsigned aa);

};

#endif
