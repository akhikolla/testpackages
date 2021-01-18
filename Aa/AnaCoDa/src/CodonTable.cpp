#include "include/CodonTable.h"
#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif


//---- OPERATOR AND CONSTRUCTOR ------
CodonTable::CodonTable()
{
    tableId = 1; //standard codon table by NCBI
    splitAA = true;
}


CodonTable::CodonTable(unsigned _tableId, bool _splitAA) : tableId(_tableId), splitAA(_splitAA)
{
	if (tableId == 7 || tableId == 8 || tableId == 15 || tableId == 17 || tableId == 18 || tableId == 19 || tableId == 20 || tableId > 25 || tableId < 1)
	{
		tableId = 1; //standard codon table by NCBI
		my_printError("Warning: Invalid codon table: % using default codon table (NCBI codon table 1)\n", tableId);
	}
}


CodonTable::~CodonTable()
{
	//dtor
}


CodonTable::CodonTable(const CodonTable& other)
{
	tableId = other.tableId;
	splitAA = other.splitAA;
}


CodonTable& CodonTable::operator=(const CodonTable& rhs)
{
	if (this == &rhs) return *this; // handle self assignment
	tableId = rhs.tableId;
	splitAA = rhs.splitAA;
	return *this;
}

// STATIC MEMBERS

const std::string CodonTable::Ser2 = "Z";
const std::string CodonTable::Ser1 = "J";
const std::string CodonTable::Thr4_1 = "O";
const std::string CodonTable::Thr4_2 = "B";
const std::string CodonTable::Leu1 = "U";

const std::string CodonTable::AminoAcidArray[] = {
	"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", CodonTable::Leu1, "M", "N", "P", "Q", "R",
	CodonTable::Ser1, CodonTable::Ser2, "S", "T", CodonTable::Thr4_1, CodonTable::Thr4_2, "V", "W", "Y", "X"};

/* const std::map<std::string, unsigned> CodonTable::aaToIndex = {{"A", 0}, {"C", 1}, {"D", 2}, {"E", 3}, {"F", 4},
	{"G", 5}, {"H", 6}, {"I", 7}, {"K", 8}, {"L", 9}, {"M", 10}, {"N", 11}, {"P", 12}, {"Q", 13}, {"R", 14}, {"S", 15},
	{"T", 16}, {"V", 17}, {"W", 18}, {"Y", 19}, {CodonTable::Ser2, 20}, {"X", 21}};
*/

const std::map<std::string, unsigned> CodonTable::aaToIndex = {{"A", 0}, {"C", 1}, {"D", 2}, {"E", 3}, {"F", 4},
 {"G", 5}, {"H", 6}, {"I", 7}, {"K", 8}, {"L",9}, {CodonTable::Leu1, 10}, {"M", 11}, {"N", 12}, {"P", 13}, {"Q", 14}, {"R", 15},
															   {CodonTable::Ser1, 16}, {CodonTable::Ser2, 17}, {"S", 18},
 {"T", 19}, {CodonTable::Thr4_1, 20}, {CodonTable::Thr4_2, 21}, {"V", 22}, {"W", 23}, {"Y", 24}, {"X", 25}};

const std::map<std::string, unsigned> CodonTable::codonToIndexWithReference = {{"GCA", 0}, {"GCC", 1}, {"GCG", 2},
	{"GCT", 3}, {"TGC", 4}, {"TGT", 5}, {"GAC", 6}, {"GAT", 7}, {"GAA", 8}, {"GAG", 9}, {"TTC", 10}, {"TTT", 11},
	{"GGA", 12}, {"GGC", 13}, {"GGG", 14}, {"GGT", 15}, {"CAC", 16}, {"CAT", 17}, {"ATA", 18}, {"ATC", 19}, {"ATT", 20},
	{"AAA", 21}, {"AAG", 22}, {"CTA", 23}, {"CTC", 24}, {"CTG", 25}, {"CTT", 26}, {"TTA", 27}, {"TTG", 28}, {"ATG", 29},
	{"AAC", 30}, {"AAT", 31}, {"CCA", 32}, {"CCC", 33}, {"CCG", 34}, {"CCT", 35}, {"CAA", 36}, {"CAG", 37}, {"AGA", 38},
	{"AGG", 39}, {"CGA", 40}, {"CGC", 41}, {"CGG", 42}, {"CGT", 43}, {"TCA", 44}, {"TCC", 45}, {"TCG", 46}, {"TCT", 47},
	{"ACA", 48}, {"ACC", 49}, {"ACG", 50}, {"ACT", 51}, {"GTA", 52}, {"GTC", 53}, {"GTG", 54}, {"GTT", 55}, {"TGG", 56},
	{"TAC", 57}, {"TAT", 58}, {"AGC", 59}, {"AGT", 60}, {"TAA", 61}, {"TAG", 62}, {"TGA", 63}};

const std::map<std::string, unsigned> CodonTable::codonToIndexWithoutReference = {{"GCA", 0}, {"GCC", 1},
	{"GCG", 2}, {"TGC", 3}, {"GAC", 4}, {"GAA", 5}, {"TTC", 6}, {"GGA", 7}, {"GGC", 8}, {"GGG", 9}, {"CAC", 10},
	{"ATA", 11}, {"ATC", 12}, {"AAA", 13}, {"CTA", 14}, {"CTC", 15}, {"CTG", 16}, {"CTT", 17}, {"TTA", 18}, {"AAC", 19},
	{"CCA", 20}, {"CCC", 21}, {"CCG", 22}, {"CAA", 23}, {"AGA", 24}, {"AGG", 25}, {"CGA", 26}, {"CGC", 27}, {"CGG", 28},
	{"TCA", 29}, {"TCC", 30}, {"TCG", 31}, {"ACA", 32}, {"ACC", 33}, {"ACG", 34}, {"GTA", 35}, {"GTC", 36}, {"GTG", 37},
	{"TAC", 38}, {"AGC", 39}};

const std::string CodonTable::codonArray[] =
	{"GCA", "GCC", "GCG", "GCT", "TGC", "TGT", "GAC", "GAT", "GAA", "GAG",
	"TTC", "TTT", "GGA", "GGC", "GGG", "GGT", "CAC", "CAT", "ATA", "ATC",
	"ATT", "AAA", "AAG", "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "ATG",
	"AAC", "AAT", "CCA", "CCC", "CCG", "CCT", "CAA", "CAG", "AGA", "AGG",
	"CGA", "CGC", "CGG", "CGT", "TCA", "TCC", "TCG", "TCT", "ACA", "ACC",
	"ACG", "ACT", "GTA", "GTC", "GTG", "GTT", "TGG", "TAC", "TAT", "AGC",
	"AGT", "TAA", "TAG", "TGA"};

const std::string CodonTable::codonArrayParameter[] =
	{"GCA", "GCC", "GCG", "TGC", "GAC",
	"GAA", "TTC", "GGA", "GGC", "GGG",
	"CAC", "ATA", "ATC", "AAA", "CTA",
	"CTC", "CTG", "CTT", "TTA", "AAC",
	"CCA", "CCC", "CCG", "CAA", "AGA",
	"AGG", "CGA", "CGC", "CGG", "TCA",
	"TCC", "TCG", "ACA", "ACC", "ACG",
	"GTA", "GTC", "GTG", "TAC", "AGC"};

// TODO NOTE: THERE IS NO CODON TABLE 7, 8, 15, 17, 18, 19, 20 ACCORDING TO NCBI !
// http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
const std::string CodonTable::codonTableDefinition[] = {"1. The Standard Code", "2. The Vertebrate Mitochondrial Code",
	"3. The Yeast Mitochondrial Code", "4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
	"5. The Invertebrate Mitochondrial Code", "6. The Ciliate, Dasycladacean and Hexamita Nuclear Code", "7. Invalid Codon Table", "8. Invalid Codon Table",
	"9. The Echinoderm and Flatworm Mitochondrial Code", "10. The Euplotid Nuclear Code",
	"11. The Bacterial, Archaeal and Plant Plastid Code", "12. The Alternative Yeast Nuclear Code", "13. The Ascidian Mitochondrial Code",
	"14. The Alternative Flatworm Mitochondrial Code", "15. Invalid Codon Table", "16. Chlorophycean Mitochondrial Code",
	"17. Invalid Codon Table", "18. Invalid Codon Table", "19. Invalid Codon Table", "20. Invalid Codon Table",
	"21. Trematode Mitochondrial Code", "22. Scenedesmus obliquus Mitochondrial Code", "23. Thraustochytrium Mitochondrial Code",
	"24. Pterobranchia Mitochondrial Code",	"25. Candidate Division SR1 and Gracilibacteria Code"};


// {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "Z", "X"};
const unsigned CodonTable::numCodonsPerAAForTable[25][26] = {
		{4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,3}, // 1. The Standard Code
		{4,2,2,2,2,4,2,2,2,6,0,2,2,4,2,4,0,2,4,4,0,0,4,2,2,4}, // 2. The Vertebrate Mitochondrial Code
		{4,2,2,2,2,4,2,2,2,2,0,2,2,4,2,4,0,2,4,0,4,4,4,2,2,2}, // 3. The Yeast Mitochondrial Code
		{4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,2,2,2}, // 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
		{4,2,2,2,2,4,2,2,2,6,0,2,2,4,2,4,0,4,4,4,0,0,4,2,2,2}, // 5. The Invertebrate Mitochondrial Code
		{4,2,2,2,2,4,2,3,2,6,0,1,2,4,4,6,0,2,4,4,0,0,4,1,2,1}, // 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 7. Invalid Codon Table
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 8. Invalid Codon Table
		{4,2,2,2,2,4,2,3,1,6,0,1,3,4,2,4,0,4,4,4,0,0,4,2,2,2}, // 9. The Echinoderm and Flatworm Mitochondrial Code
		{4,3,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2}, // 10. The Euplotid Nuclear Code
		{4,2,2,2,2,4,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,3}, // 11. The Bacterial, Archaeal and Plant Plastid Codee
		{4,2,2,2,2,4,2,3,2,5,0,1,2,4,2,6,1,2,4,4,0,0,4,1,2,3}, // 12. The Alternative Yeast Nuclear Code
		{4,2,2,2,2,6,2,2,2,6,0,2,2,4,2,4,0,2,4,4,0,0,4,2,2,2}, // 13. The Ascidian Mitochondrial Code
		{4,2,2,2,2,4,2,3,1,6,0,1,3,4,2,4,0,4,4,4,0,0,4,2,3,1}, // 14. The Alternative Flatworm Mitochondrial Code
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 15. Invalid Codon Table
		{4,2,2,2,2,4,2,3,2,6,1,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2}, // 16. Chlorophycean Mitochondrial Code
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 17. Invalid Codon Table
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 18. Invalid Codon Table
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 19. Invalid Codon Table
		{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 20. Invalid Codon Table
		{4,2,2,2,2,4,2,2,1,6,0,2,3,4,2,4,0,4,4,4,0,0,4,2,2,2}, // 21. Trematode Mitochondrial Code
		{4,2,2,2,2,4,2,3,2,6,1,1,2,4,2,6,0,2,3,4,0,0,4,1,2,3}, // 22. Scenedesmus obliquus Mitochondrial Code
		{4,2,2,2,2,4,2,3,2,5,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,4}, // 23. Thraustochytrium Mitochondrial Code
		{4,2,2,2,2,4,2,3,3,6,0,1,2,4,2,4,0,3,4,4,0,0,4,2,2,2}, // 24. Pterobranchia Mitochondrial Code
		{4,2,2,2,2,5,2,3,2,6,0,1,2,4,2,6,0,2,4,4,0,0,4,1,2,2}, // 25. Candidate Division SR1 and Gracilibacteria Code
};

// --- CODON TABLE SPECIFIC MAPPER FUNCTIONS ------

unsigned CodonTable::AAToAAIndex(std::string aa)
{
        return CodonTable::aaToIndex.find(aa) -> second;
}

unsigned CodonTable::getNumCodons(std::string aa)
{
	return getNumCodons(AAToAAIndex(aa));
}

unsigned CodonTable::getNumCodons(unsigned aa)
{
	return numCodonsPerAAForTable[tableId][aa];
}

void CodonTable::setupCodonTable()
{
	unsigned numAA = 21;
	if (tableId >= 1 && tableId <= 6  && !splitAA) numAA = 20; //If not splitting AAs, tables 1-6 all have 20 AA, excluding stop codes
    else if (tableId >= 9 && tableId <= 14  && !splitAA) numAA = 20;
    else if (tableId == 16 && !splitAA) numAA = 20;
    else if (tableId >= 21 && tableId <= 25 && !splitAA) numAA = 20;


	codon_mapping.resize(numAA);

    unsigned aaIndex = 0;
    unsigned filled = 0;
    while (filled != numAA)
    {
        unsigned numCodons = getNumCodons(aaIndex);
        if (numCodons != 0)
        {
            if (aaIndex == 0) //A
            {
                for (unsigned i = 0; i < 4; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 1) //C
            {
                for (unsigned i = 4; i < 6; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
                if (tableId == 10)
                {
                    codon_mapping[filled].push_back(63);
                }
            }
            else if (aaIndex == 2) //D
            {
                for (unsigned i = 6; i < 8; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 3) //E
            {
                for (unsigned i = 8; i < 10; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 4) //F
            {
                for (unsigned i = 10; i < 12; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 5) //G
            {
                for (unsigned i = 12; i < 16; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
                if (tableId == 13)
                {
                    codon_mapping[filled].push_back(38);
                    codon_mapping[filled].push_back(39);
                }
                else if (tableId == 25)
                {
                    codon_mapping[filled].push_back(63);
                }
            }
            else if (aaIndex == 6) //H
            {
                for (unsigned i = 16; i < 18; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 7) //I
            {
                if (tableId == 1 || tableId == 4 || tableId == 6 || tableId == 9 || tableId == 10 || tableId == 11
                        || tableId == 12 || tableId == 14 || tableId == 16 || tableId >= 22)
                {
                    codon_mapping[filled].push_back(18);
                }
                for (unsigned i = 19; i < 21; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 8) //K
            {
                if ((tableId >= 1 && tableId <= 6) || (tableId >= 10 && tableId <= 13) || (tableId == 16) ||
                        (tableId >= 22))
                {
                    for (unsigned i = 21; i < 23; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                }
                else if (tableId == 9 || tableId == 14 || tableId == 21)
                {
                    codon_mapping[filled].push_back(22);
                }

                if (tableId == 24)
                {
                    codon_mapping[filled].push_back(39);
                }
            }
            else if (aaIndex == 9) //L
            {
                if (tableId != 3 && tableId != 23)
                {
                    for (unsigned i = 23; i < 29; i++)
                    {
                        if (i == 25 && tableId != 12)
                        {
                            codon_mapping[filled].push_back(i);
                        }
                        else
                        {
                            codon_mapping[filled].push_back(i);
                        }
                    }
                    if ((tableId == 16 || tableId == 22) && !splitAA)
                    {
                        codon_mapping[filled].push_back(62);
                    }
                }
                else if (tableId == 3)
                {
                    for (unsigned i = 27; i < 29; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                }
                else if (tableId == 23)
                {
                    for (unsigned i = 23; i < 27; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                    codon_mapping[filled].push_back(28);
                }
            }
            else if (aaIndex == 10) //Leu1
            {
                if (splitAA)
                {
                    codon_mapping[filled].push_back(62);
                }
                else
                {
                    filled--;
                }
            }
            else if (aaIndex == 11) //M
            {
                codon_mapping[filled].push_back(29);
                if (tableId == 2 || tableId == 3 || tableId == 5 || tableId == 13 || tableId == 21)
                {
                    codon_mapping[filled].push_back(18);
                }
            }
            else if (aaIndex == 12) //N
            {
                for (unsigned i = 30; i < 32; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
                if (tableId == 9 || tableId == 14 || tableId == 21)
                {
                    codon_mapping[filled].push_back(21);
                }
            }
            else if (aaIndex == 13) //P
            {
                for (unsigned i = 32; i < 36; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 14) //Q
            {
                for (unsigned i = 36; i < 38; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
                if (tableId == 6)
                {
                    codon_mapping[filled].push_back(61);
                    codon_mapping[filled].push_back(62);
                }
            }
            else if (aaIndex == 15) //R
            {
                if (tableId == 1 || tableId == 3 || tableId == 4 || tableId == 6 || (tableId >= 10 && tableId <= 12) || tableId == 16 ||
                        tableId == 22 || tableId == 23 || tableId == 25)
                {
                    codon_mapping[filled].push_back(38);
                    codon_mapping[filled].push_back(39);
                }

                if (tableId == 3)
                {
                    codon_mapping[filled].push_back(42);
                    codon_mapping[filled].push_back(43);
                }
                else
                {
                    for (unsigned i = 40; i < 44; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                }
            }
            else if (aaIndex == 16) //Ser1
            {
                if (splitAA)
                {
                    codon_mapping[filled].push_back(25); //table 12
                }
                else
                {
                    filled--;
                }
            }
            else if (aaIndex == 17) //Ser2
            {
                if (splitAA)
                {
                    codon_mapping[filled].push_back(59);
                    codon_mapping[filled].push_back(60);

                    if (tableId == 5 || tableId == 9 || tableId == 14 || tableId == 21)
                    {
                        codon_mapping[filled].push_back(38);
                        codon_mapping[filled].push_back(39);
                    }
                    else if (tableId == 24)
                    {
                        codon_mapping[filled].push_back(38);
                    }
                }
                else
                {
                    filled--;
                }
            }
            else if (aaIndex == 18) //Ser (Ser4)
            {
                if (tableId != 22)
                {
                    for (unsigned i = 44; i < 48; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                }
                else
                {
                    for (unsigned i = 45; i < 48; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                }

                if (!splitAA)
                {
                    if (tableId == 12)
                    {
                        codon_mapping[filled].push_back(25);
                    }

                    codon_mapping[filled].push_back(59);
                    codon_mapping[filled].push_back(60);

                    if (tableId == 5 || tableId == 9 || tableId == 14 || tableId == 21)
                    {
                        codon_mapping[filled].push_back(38);
                        codon_mapping[filled].push_back(39);
                    }
                    else if (tableId == 24)
                    {
                        codon_mapping[filled].push_back(38);
                    }
                }
            }
            else if (aaIndex == 19) //T
            {
                if (tableId != 3)
                {
                    for (unsigned i = 48; i < 52; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                }
                else if (tableId == 3 && !splitAA)
                {
                    for (unsigned i = 48; i < 52; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }

                    for (unsigned i = 23; i < 27; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                }
            }
            else if (aaIndex == 20) //Thr1
            {
                if (splitAA)
                {
                    for (unsigned i = 48; i < 52; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                }
                else
                {
                    filled--;
                }
            }
            else if (aaIndex == 21) //Thr2
            {
                if (splitAA)
                {
                    for (unsigned i = 23; i < 27; i++)
                    {
                        codon_mapping[filled].push_back(i);
                    }
                }
                else
                {
                    filled--;
                }
            }
            else if (aaIndex == 22) //V
            {
                for (unsigned i = 52; i < 56; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
            }
            else if (aaIndex == 23) //W
            {
                codon_mapping[filled].push_back(56);
                if ((tableId >= 2 && tableId <= 5) || (tableId == 9) || (tableId == 13) || (tableId == 14) || (tableId == 21))
                {
                    codon_mapping[filled].push_back(63);
                }
            }
            else if (aaIndex == 24) //Y
            {
                for (unsigned i = 56; i < 58; i++)
                {
                    codon_mapping[filled].push_back(i);
                }
                if (tableId == 14) codon_mapping[filled].push_back(61);
            }
            //ignoring stop amino acid for now
            filled++;
        }
        aaIndex++;
    }
}



