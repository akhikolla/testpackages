#include "include/Genome.h"

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif





//--------------------------------------------------//
// ---------- Constructors & Destructors ---------- //
//--------------------------------------------------//

/* Gene constructor (RCPP EXPOSED)
 * Arguments: None
 * Blank constructor for Gene class. Sets the sequence, id, and
 * description fields to empty strings.
*/
Genome::Genome()
{
	//ctor
}


Genome& Genome::operator=(const Genome& rhs)
{
	if (this == &rhs) return *this; // handle self assignment
	genes = rhs.genes;
	simulatedGenes = rhs.simulatedGenes;
	numGenesWithPhi = rhs.numGenesWithPhi;
	RFPCountColumnNames = rhs.RFPCountColumnNames;
	prev_genome_size = rhs.prev_genome_size;
	//assignment operator
	return *this;
}


Genome::~Genome()
{
	//dtor
}


bool Genome::operator==(const Genome& other) const
{
	bool match = true;

 	//Do a ! operation because only the gene comparison is implemented, not the != operator.
	if (!(this->genes == other.genes)) { match = false; }
	if (!(this->simulatedGenes == other.simulatedGenes)) { match = false; }
	if (this->numGenesWithPhi != other.numGenesWithPhi) { match = false; }
	if (this->RFPCountColumnNames != other.RFPCountColumnNames) { match = false; }

	return match;
}


//----------------------------------------//
//---------- File I/O Functions ----------//
//----------------------------------------//

/* readFasta (RCPP EXPOSED)
 * Arguments: string filename, boolean to determine if we are appending to an existing Fasta sequence
 * (if not set to true, will default to clearing file; defaults to false)
 * Takes input in Fasta format from file and saves to genome.
*/
void Genome::readFasta(std::string filename, bool append)
{
	prev_genome_size = genes.size();
	try
	{
		if (!append)
			clear();
		std::ifstream Fin;
		Fin.open(filename.c_str());
		if (Fin.fail())
			my_printError("ERROR: Error in Genome::readFasta: Can not open Fasta file %\n", filename);
		else
		{
			//my_print("File opened\n");
			bool fastaFormat = false;
			std::string buf;
			int newLine;

			Gene tmpGene;
			std::string tempSeq = "";
			while (true)
			{
			    #ifndef STANDALONE
			    Rcpp::checkUserInterrupt();
			    #endif

				// read a new line in every cycle
				std::getline(Fin, buf);
				/* The new line is marked with an integer, which corresponds
					 to one of the following cases:

					 1. New sequence line, which starts with "> ".
					 2. Sequence line, which starts with any character except for ">"
					 3. End of file, detected by function eof().
				 */

				if (buf[0] == '>' ) newLine=1;
				else if (Fin.eof()) newLine=3;
				else newLine=2;

				if ( newLine == 1 )
				{ // this is a start of a new chain.
					if (!fastaFormat)
					{
						// if it is the first chain, just build a new chain object
						tmpGene.clear();
						fastaFormat = true;
					}
					else
					{
						// otherwise, need to store the old chain first, then build a new chain
						tmpGene.setSequence(tempSeq);
						addGene(tmpGene);

						tmpGene.clear();
						tempSeq = "";
					}
					tmpGene.setDescription( buf.substr(1,buf.size()-1) );
					std::size_t pos = buf.find(' ') - 1;
					tmpGene.setId( buf.substr(1,pos) );
				}

				if ( newLine == 2 && fastaFormat )
				{ // sequence line
					tempSeq.append(buf);
				}

				if ( newLine == 3 )
				{ // end of file
					//my_print("EOF reached\n");
					if ( !fastaFormat )
						throw std::string("Genome::readFasta throws: ") + std::string(filename)
							  + std::string(" is not in Fasta format.");
					else
					{
						// otherwise, need to store the old chain first, then to
						// build a new chain
						tmpGene.setSequence(tempSeq);
						addGene(tmpGene);
						break;
					}
				}
			} // end while
		} // end else
		Fin.close();
	} // end try
	catch(char* pMsg)
	{
		my_printError("\nException: %\n", pMsg);
	}
}


/* writeFasta (RCPP EXPOSED)
 * Arguments: filename to write to,
 * boolean specifying if the genome is simulated or not (default non-simulated).
 * Writes a genome in Fasta format to given file
*/
void Genome::writeFasta (std::string filename, bool simulated)
{
	try {
		std::ofstream Fout;
		Fout.open(filename.c_str());
		if (Fout.fail())
			my_printError("Error in Genome::writeFasta: Can not open output Fasta file %\n", filename);
		else
		{
			unsigned sized = simulated ? (unsigned)simulatedGenes.size() : (unsigned)genes.size();

			for (unsigned i = 0u; i < sized; i++)
			{
				Gene *currentGene = simulated ? &simulatedGenes[i] : &genes[i];

				Fout << ">" << currentGene->getDescription() << "\n";
				for (unsigned j = 0u; j < currentGene->length(); j++)
				{
					Fout << currentGene->getNucleotideAt(j);
					if ((j + 1) % 60 == 0) Fout << std::endl;
				}
				Fout << std::endl;
			}
		} // end else
		Fout.close();
	} // end try
	catch(char* pMsg)
	{
		my_printError("\nException: %\n", pMsg);
	}
}


// Recall: Because this is simulated genome data from PAModel's simulateGenome function, it will read in data that
// disregards position. Codon counts are summed up when setSequence calls processSequence.
// Note that for debug purposes this reads genome into the non-simulated genome TODO fix that
void Genome::readSimulatedGenomeFromPAModel(std::string filename)
{
	std::ifstream Fin;
	Fin.open(filename.c_str());
	if (Fin.fail())
		my_printError("Error in Genome::readSimulatedGenomeFromPAModel: Can not open file %\n", filename);

	std::string tmp;
	// Trash the first line, but we know it is in a simulatedGenomeFromPAModel format: GeneID,Codon,Codon_Count,RFPCount
	if (!std::getline(Fin, tmp))
		my_printError("Error in Genome::readSimulatedGenomeFromPAModel: RFP file % has no header.\n", filename);

	std::string prevID = "";
	Gene tmpGene;
	bool first = true;
	std::string seq = "";
	totalRFPCount = 0;

	while (std::getline(Fin, tmp))
	{
            #ifndef STANDALONE
            Rcpp::checkUserInterrupt();
            #endif

    // Remove whitespace from the string
    tmp.erase(std::remove_if(tmp.begin(), tmp.end(),
              std::bind(std::isspace<char>, std::placeholders::_1, std::locale::classic())), tmp.end());

		// Find the GeneID
		std::size_t pos = tmp.find(',');
		std::string ID = tmp.substr(0, pos);

		if (first)
		{
			prevID = ID;
			first = false;
		}
		if (ID != prevID)
		{
			//my_print("I am clearing the Gene sequence!\n");
			tmpGene.setId(prevID);
			tmpGene.setDescription("No description for PA(NSE) Model");
			tmpGene.setSequence(seq);
			addGene(tmpGene, true); //add to genome
			tmpGene.clear();
			seq = "";
		}

        // Now find codon string: Follows prior comma, guaranteed to be of size 3
        std::string codon = tmp.substr(pos + 1, 3);
        codon[0] = (char)std::toupper(codon[0]);
        codon[1] = (char)std::toupper(codon[1]);
        codon[2] = (char)std::toupper(codon[2]);

        std::size_t pos2 = pos + 4;
        pos = tmp.find(",", pos2 + 1);
        std::string value = tmp.substr(pos2 + 1, pos - (pos2 + 1));
        unsigned codonCount = (unsigned)std::atoi(value.c_str());

        // Concatenate the sequence based on this number of codons
        for (unsigned i = 0u; i < codonCount; i++)
            seq.append(codon);

        //Read rfp count
        value = tmp.substr(pos + 1);
		unsigned rfpcount = (unsigned)std::atoi(value.c_str());

		unsigned codonIndex = SequenceSummary::codonToIndex(codon);
		tmpGene.geneData.setCodonSpecificSumRFPCount(codonIndex, rfpcount, 0);
		prevID = ID;
	}

	tmpGene.setId(prevID);
	tmpGene.setDescription("No description for PA(NSE) Model");
	tmpGene.setSequence(seq);
	addGene(tmpGene, true); //add to genome
	totalRFPCount += tmpGene.geneData.getSumTotalRFPCount(0);

	Fin.close();
}


/* readRFPData (RCPP EXPOSED)
 * Arguments: string filename, boolean to determine if we are appending to the existing genome
 * (if not set to true, will default to clearing genome data; defaults to false)
 * Read in a PA-formatted file: GeneID,Position (1-indexed),Codon,RFPCount(s) (may be multiple)
 * The positions are not necessarily in the right order.
 * Ignores ambiguously-positioned codons (marked with negative position).
 * RFPCounts are not given as a positive number are set to -1 and are stored, but not used in calculation.
 * There may be more than one RFPCount, and thus the header is important.
*/
void Genome::readRFPData(std::string filename, bool append, bool positional)
{
	prev_genome_size = genes.size();
	try {
		if (!append) clear();
		totalRFPCount = 0;
		std::ifstream Fin;
		Fin.open(filename.c_str());

		if (Fin.fail())
			my_printError("Error in Genome::readRFPData: Can not open RFPData file %\n", filename);
		else
		{
			// Analyze the header line
			std::string tmp;

            if (!std::getline(Fin, tmp))
                my_printError("Error in Genome::readRFPData: RFPData file % has no header.\n", filename);

			// Ignore first 3 commas: ID, position, codon
			std::size_t pos = tmp.find(',');
			pos = tmp.find(',', pos + 1);
			pos = tmp.find(',', pos + 1);

			// While there are more commas, there are more categories of RFP counts
			unsigned numCategories = 0;
			std::size_t pos2;
			while (pos != std::string::npos)
			{
				numCategories++;
				pos2 = tmp.find(',', pos + 1);
				addRFPCountColumnName(tmp.substr(pos + 1, pos2 - (pos + 1)));
				pos = pos2;
			}
			unsigned tableWidth = 2 + numCategories; // Table size: position + codonID + each category

			// -------------------------------------------------------------------------//
			// --------------- Now for each RFPCount category, record it. ------------- //
			// ----- Treat data as a table: push back vectors of size tableWidth. ----- //
			// -------------------------------------------------------------------------//
			std::string prevID = "";
			Gene tmpGene;
			int position, possValue;
			prevID = "";
			bool first = true;

			std::vector <std::vector <int>> table; // Dimensions: nRows of tableWidth-sized vectors

			//Now for each line associated with a gene ID, set the string appropriately
			while (std::getline(Fin, tmp))
			{

	            #ifndef STANDALONE
	            Rcpp::checkUserInterrupt();
	            #endif

		        // Remove whitespace from the string
		        tmp.erase(std::remove_if(tmp.begin(), tmp.end(),
	                  std::bind(std::isspace<char>, std::placeholders::_1, std::locale::classic())), tmp.end());

				pos = tmp.find(',');
				std::string ID = tmp.substr(0, pos);

				// Make pos2 find next comma to find position integer
				pos2 = tmp.find(',', pos + 1);

				// Set to 0-indexed value for convenience of vector calculation
				// Error-checking note: Of course, atoi function returns 0 if not an integer
				// -> leads to 0 - 1 == -1 -> codon ignored (as intended).
				position = std::atoi(tmp.substr(pos + 1, pos2 - (pos + 1)).c_str()) - 1;

				// Position integer: Ensure that if position is negative, ignore the codon.
				if (position > -1) // for convenience of calculation; ambiguous positions are marked by -1.
				{
					std::vector <int> tableRow(tableWidth);
					unsigned tableIndex = 0;
					tableRow[tableIndex] = position;

					if (first)
					{
						prevID = ID;
						first = false;
					}
					//This is when we start at a new geneID
					if (ID != prevID)
					{
						tmpGene.setId(prevID);
						tmpGene.setDescription("No description for PA(NSE) Model");
						if(positional) 
						{
							tmpGene.setPANSESequence(table);
						}
						else 
						{
							tmpGene.setPASequence(table);
						}
						addGene(tmpGene, false); //add to genome
						totalRFPCount += tmpGene.geneData.getSumTotalRFPCount(0);
						tmpGene.clear();
						table.clear();
					}
					// Now find codon string: Follows prior comma, guaranteed to be of size 3
					std::string codon = tmp.substr(pos2 + 1, 3);
					codon[0] = (char)std::toupper(codon[0]);
					codon[1] = (char)std::toupper(codon[1]);
					codon[2] = (char)std::toupper(codon[2]);

					tableIndex ++;
					tableRow[tableIndex] = SequenceSummary::codonToIndex(codon);
					// Note: May be an invalid Codon read in, but this is resolved when the sequence is set and processed.

					// Skip to end RFPCount(s), if any
					pos = tmp.find(',', pos2 + 1);
					tableIndex ++;

					// While there are more commas, there are more categories of RFP counts
					while (pos != std::string::npos)
					{
						pos2 = tmp.find(',', pos + 1);

						// RFPCount is input as either its integer value (> -1) or as -1 (stored, not calculated).
						// Also accounts for if the value is "NA" or a string (converted to -1).
						if (sscanf(tmp.substr(pos + 1, pos2 - (pos + 1)).c_str(), "%d", &possValue) == 1)
						{
							if (possValue > -1) tableRow[tableIndex] = possValue;
							else tableRow[tableIndex] = -1;
						}
                        else
						{
							tableRow[tableIndex] = -1;
						}

						pos = pos2;
						tableIndex++;
					}

					prevID = ID;
					table.push_back(tableRow);
				}
			}

            // Ensure that at least one entry was read in
            if (prevID != "")
            {
                tmpGene.setId(prevID);
                tmpGene.setDescription("No description for PA(NSE) Model");
                tmpGene.setPASequence(table);
                totalRFPCount += tmpGene.geneData.getSumTotalRFPCount(0);
                addGene(tmpGene, false); //add to genome
            }
		} // end else
		Fin.close();
	} // end try
	catch(char* pMsg)
	{
		my_printError("\nException: %\n", pMsg);
	}
}


/* writeRFPData (RCPP EXPOSED)
 * Arguments: string filename, boolean to specify if we are printing simulated genes or not (default non-simulated)
 * Write a PA-formatted file: GeneID,Position (1-indexed),Codon,RFPCount(s) (may be multiple)
 * The positions will be printed in ascending order.
*/
void Genome::writeRFPData(std::string filename, bool simulated)
{
	std::ofstream Fout;
	Fout.open(filename.c_str());
	if (Fout.fail())
		my_printError("Error in Genome::writeRFPData: Can not open output RFPData file %\n", filename);
	else
	{
		unsigned numGenes = simulated ? (unsigned)simulatedGenes.size() : (unsigned)genes.size();

		if (!simulated)
		{
			Fout << "GeneID,Position,Codon";

			// For each category name, print for header
			std::vector<std::string> RFPCountColumnNames = getRFPCountColumnNames();
			unsigned numCategories = (unsigned)RFPCountColumnNames.size();
			for (unsigned category = 0u; category < numCategories; category++)
				Fout << "," << RFPCountColumnNames[category];

			Fout << "\n";

			for (unsigned geneIndex = 0u; geneIndex < numGenes; geneIndex++)
			{
				Gene *currentGene = &genes[geneIndex];
				std::vector<unsigned> positionCodonID = currentGene->geneData.getPositionCodonID();
				unsigned numPositions = (unsigned)positionCodonID.size();

				for (unsigned position = 0u; position < numPositions; position++)
				{
					unsigned codonID = positionCodonID[position];
					std::string codon = SequenceSummary::codonArray[codonID];

					// Print position + 1 because it's externally one-indexed
					Fout << currentGene->getId() << "," << position + 1 << "," << codon;

					for (unsigned category = 0; category < numCategories; category++)
						Fout << "," << currentGene->geneData.getSingleRFPCount(position, category);

					Fout << "\n";
				}
			}
		}
		else
		{
			Fout << "GeneID,Position,Codon,RFPCount\n";
			for (unsigned geneIndex = 0u; geneIndex < numGenes; geneIndex++)
			{
                //unsigned position = 1u;
				Gene *currentGene = &simulatedGenes[geneIndex];
                SequenceSummary *sequenceSummary = currentGene->getSequenceSummary();
                std::vector <unsigned> positions = sequenceSummary->getPositionCodonID();
                std::vector <unsigned> rfpCounts = sequenceSummary->getRFPCount(0);
				for (unsigned positionIndex = 0u; positionIndex < positions.size(); positionIndex++)
				{
					unsigned codonID = positions[positionIndex];
    				std::string codon = SequenceSummary::codonArray[codonID];

    				Fout << currentGene->getId() << "," << positionIndex + 1 << "," << codon << "," << rfpCounts[positionIndex] << "\n";

				}
			}
		}
	}
	Fout.close();
}


/* readObservedPhiValues
 * Arguments: string filename, bool byId
 *
 * User may choose to store Phi values either by ID (if byId is true) or by index (if false).
 * In by index, each gene read in should correspond numerically to the genome.
 * Values that less than or equal to 0 or not a number will be converted to -1 and not
 * counted as a phi value.
 * NOTE: IF AN ERROR FILE IS READ IN, numGenesWithPhi IS STILL INITIALIZED WITH 0'S.
*/
void Genome::readObservedPhiValues(std::string filename, bool byId)
{
	std::ifstream input;
	std::string tmp;
	unsigned numPhi = 0;
	bool exitFunction = false;

	input.open(filename);
	if (input.fail())
		my_printError("ERROR in Genome::readObservedPhiValues: Can not open file %\n", filename);
	else
	{
		//Trash the header line
		if (!std::getline(input, tmp))
			my_printError("Error in Genome::readObservedPhiValues File: File % has no header.\n", filename);

		if (genes.size() == 0)
			my_printError("ERROR: Genome is empty, function will not execute!\n");
		else
		{
			if (byId)
			{
				//Mapping is done so the genes can be found
				std::map<std::string, Gene *> genomeMapping;
				for (unsigned i = 0; i < genes.size(); i++)
				{
					genomeMapping.insert(make_pair(genes[i].getId(), &genes[i]));
				}

                bool first = true;
                while (std::getline(input, tmp))
				{

            #ifndef STANDALONE
            Rcpp::checkUserInterrupt();
            #endif
                    std::size_t pos = tmp.find(',');
                    std::string geneID = tmp.substr(0, pos);
                    std::map<std::string, Gene *>::iterator it;
                    it = genomeMapping.find(geneID);

                    if (it == genomeMapping.end())
                        my_printError("WARNING: Gene % not found!\n", geneID);

                    else //gene is found
                    {
                        std::string val = "";
                        bool notDone = true;
						double dval;

						// Loop over each comma-separated value
                        while (notDone)
						{
							std::size_t pos2 = tmp.find(',', pos + 1);

							// End of string, end loop
							if (pos2 == std::string::npos)
							{
								val = tmp.substr(pos + 1, tmp.size() + 1);
								notDone = false;
							}
							else
							{
								val = tmp.substr(pos + 1, pos2 - (pos + 1));
								pos = pos2;
							}
							dval = std::atof(val.c_str());

							//If the value is negative, nan, or 0, we set it to -1 and throw a warning message.
							if (dval <=  0 || std::isnan(dval))
							{
								dval = -1;
								my_printError("WARNING! Invalid, negative, or 0 phi value given; values should not be on the log scale. Missing value flag stored.\n");
							}
							it->second->observedSynthesisRateValues.push_back(dval);
						}

						// If this is the first value, initialize the size of numGenesWithPhi
						if (first)
						{
							first = false;
							numPhi = (unsigned) it->second->observedSynthesisRateValues.size();
							numGenesWithPhi.resize(numPhi, 0);
						}
						else if (it->second->observedSynthesisRateValues.size() != numPhi)
						{
                            my_printError("ERROR: Gene % has a different number of phi values given other genes: \n", geneID);
                            my_printError("Gene % has % ", geneID, it -> second -> observedSynthesisRateValues.size());
                            my_printError("while others have %. Exiting function.\n", numPhi);
							exitFunction = true;
							for (unsigned a = 0; a < getGenomeSize(); a++)
							{
								genes[a].observedSynthesisRateValues.clear();
							}
							break;
						}
                    }
                    if (exitFunction) break;
                }
				// If number of phi values match those of other given genes, execution continues normally.
				// Proceed to increment numGenesWithPhi.
                if (!exitFunction)
				{
                    for (unsigned i = 0; i < getGenomeSize(); i++)
					{
						Gene *gene = &(getGene(i));
                        if (gene->observedSynthesisRateValues.size() != numPhi)
						{
							my_printError("WARNING: Gene # % (%) does not have any phi values. ", i, gene->getId());
                            my_printError("Please check your file to make sure every gene has a phi value. Filling empty genes ");
                            my_printError("with Missing Value Flag for calculations.\n");
                           	gene->observedSynthesisRateValues.resize(numPhi, -1);

                        }

						// Finally increment numGenesWithPhi based on stored observedSynthesisRateValues
						for (unsigned j = 0; j < numPhi; j++)
						{
							if (gene->observedSynthesisRateValues[j] != -1)
								numGenesWithPhi[j]++;
						}
                    }
                }
            }
				// By index
            else
			{
				//unsigned geneIndex = 0;
				unsigned geneIndex=prev_genome_size;
				bool first = true;

				while (std::getline(input, tmp))
				{
					if (geneIndex >= genes.size())
					{
						my_printError("ERROR: GeneIndex exceeds the number of genes in the genome. Exiting function.");
						break;
					}

					std::size_t pos = tmp.find(',');
					std::string val = "";
					bool notDone = true;
					double dval;

					// Loop over each comma-separated value
					while (notDone)
					{
						std::size_t pos2 = tmp.find(',', pos + 1);

						// End of string, end loop
						if (pos2 == std::string::npos)
						{
							val = tmp.substr(pos + 1, tmp.size() + 1);
							notDone = false;
						}
						else
						{
							val = tmp.substr(pos + 1, pos2 - (pos + 1));
							pos = pos2;
						}
						dval = std::atof(val.c_str());

						//If the value is negative, nan, or 0, we set it to -1 and throw a warning message.
						if (dval <= 0 || std::isnan(dval))
						{
							dval = -1;
							my_printError("WARNING! Invalid, negative, or 0 phi value given; values should not be on the log scale. Missing value flag stored.\n");
						}
						genes[geneIndex].observedSynthesisRateValues.push_back(dval);
					}

					// If this is the first value, initialize the size of numGenesWithPhi
					if (first)
					{
						first = false;
						numPhi = (unsigned) genes[geneIndex].observedSynthesisRateValues.size();
						numGenesWithPhi.resize(numPhi, 0);
					}
					else if (numPhi != genes[geneIndex].observedSynthesisRateValues.size())
					{
                        my_printError("ERROR: Gene % has a different number of phi values given other genes: \n", geneIndex);
                        my_printError("Gene % has % ", geneIndex, genes[geneIndex].observedSynthesisRateValues.size());
                        my_printError("while others have %. Exiting function.\n", numPhi);
						exitFunction = true;
						for (unsigned a = 0; a < getGenomeSize(); a++)
						{
							genes[a].observedSynthesisRateValues.clear();
						}
						break;
					}
					geneIndex++;
					if (exitFunction) break;
				}
				// If number of phi values match those of other given genes, execution continues normally.
				// Proceed to increment numGenesWithPhi.
				if (!exitFunction)
				{
					for (unsigned i = 0; i < getGenomeSize(); i++)
					{
						Gene *gene = &(getGene(i));
						if (gene->observedSynthesisRateValues.size() != numPhi)
						{
                            my_printError("WARNING: Gene # % (%) does not have any phi values. ", i, gene->getId());
                            my_printError("Please check your file to make sure every gene has a phi value. Filling empty genes ");
                            my_printError("with Missing Value Flag for calculations.\n");
                           	gene->observedSynthesisRateValues.resize(numPhi, -1);

						}

						// Finally increment numGenesWithPhi based on stored observedSynthesisRateValues
						for (unsigned j = 0; j < numPhi; j++)
						{
							if (gene->observedSynthesisRateValues[j] != -1)
								numGenesWithPhi[j]++;
						}
					}
				}
			}
		}
		input.close();
	}
}


void Genome::removeUnobservedGenes()
{
	std::vector<Gene> tmp;
	for (unsigned i = 0; i < getGenomeSize(); i++)
	{
		Gene *gene = &(getGene(i));
		for (unsigned j = 0; j < gene -> observedSynthesisRateValues.size(); j++)
		{
			if (gene -> observedSynthesisRateValues[j] != -1)
			{
				tmp.push_back(*gene);
				break;
			}
		}
	}
	genes = tmp;
}

unsigned Genome::getSumRFP()
{
    return totalRFPCount;
}
//------------------------------------//
//---------- Gene Functions ----------//
//------------------------------------//


/* Gene constructor (RCPP EXPOSED)
 * Arguments: gene to add, boolean if it was simulated.
 * Depending on whether a gene was simulated, appends to
 * genes or simulated genes.
*/
void Genome::addGene(Gene& gene, bool simulated)
{
	simulated ? simulatedGenes.push_back(gene) : genes.push_back(gene);
}

/* getGenes (RCPP EXPOSED)
 * Arguments: boolean if simulated genes should be returned.
 * Returns depending on the argument the genes or
 * simulated genes vector.
*/
std::vector <Gene> Genome::getGenes(bool simulated)
{
	return !simulated ? genes : simulatedGenes;
}


/* getNumGenesWithPhiForIndex (RCPP EXPOSED)
 * Arguments: index number.
 * Returns the number of genes with the given Phi
 * expressed as an index.
*/
unsigned Genome::getNumGenesWithPhiForIndex(unsigned index)
{
	return numGenesWithPhi[index];
}


/* getGene (unsigned index) (RCPP EXPOSED)
 * Arguments: index number, simulated
 * Returns the gene from the requested set at index
*/
Gene& Genome::getGene(unsigned index, bool simulated)
{
	return simulated ? simulatedGenes[index] : genes[index];
}


/* getGene (string id) (RCPP EXPOSED)
 * Arguments: id, simulated
 * Returns the gene from the requested set with the id given
*/
Gene& Genome::getGene(std::string id, bool simulated)
{
	Gene tempGene;
	unsigned geneIndex;

	for (geneIndex = 0; geneIndex < getGenomeSize(); geneIndex++)
	{
		tempGene = simulated ? simulatedGenes[geneIndex] : genes[geneIndex];

		if (tempGene.getId() == id) break;
	}
	return simulated ? simulatedGenes[geneIndex] : genes[geneIndex];
}





//-------------------------------------//
//---------- Other Functions ----------//
//-------------------------------------//


/* getGenomeSize (RCPP EXPOSED)
 * Arguments: boolean if requesting
 * size of simulated genes.
 * Returns the size of requested genes structure
*/
unsigned Genome::getGenomeSize(bool simulated)
{
	return simulated ? (unsigned)simulatedGenes.size() : (unsigned)genes.size();
}


/* clear (RCPP EXPOSED)
 * Arguments: None.
 * clears all data structure containing
 * gene information.
*/
void Genome::clear()
{
	genes.clear();
	simulatedGenes.clear();
	numGenesWithPhi.clear();
	RFPCountColumnNames.clear();
}


/* getGenomeForGeneIndices (RCPP EXPOSED)
 * Arguments: vector of indices, boolean if simulated.
 * Returns a genome of genes at indices.
*/
Genome Genome::getGenomeForGeneIndices(std::vector <unsigned> indices, bool simulated)
{
	Genome genome;

	for (unsigned i = 0; i < indices.size(); i++)
	{
		if (indices[i] > getGenomeSize(simulated))
		{
			my_printError("Error in Genome::getGenomeForGeneIndices. An index specified is out of bounds for the genome!\n");
			my_printError("The index % is greater than the size of the genome (%).\n", indices[i], getGenomeSize());
			my_printError("Returning empty Genome.\n");
			genome.clear();
			return genome;
		}
		else
		{
			simulated ? genome.addGene(simulatedGenes[indices[i]], true) : genome.addGene(genes[indices[i]], false);
		}
	}

	return genome;
}


/* getCodonCountsPerGene (RCPP EXPOSED)
 * Arguments: a string which is the
 * codon sequence concerning the user.
 * Returns the number of times the sequence occurs.
*/
std::vector<unsigned> Genome::getCodonCountsPerGene(std::string codon)
{
	std::vector<unsigned> codonCounts(genes.size());
	unsigned codonIndex = SequenceSummary::codonToIndex(codon);
	for (unsigned i = 0u; i < genes.size(); i++)
	{
		Gene gene = genes[i];
		SequenceSummary *sequenceSummary = gene.getSequenceSummary();
		codonCounts[i] = sequenceSummary->getCodonCountForCodon(codonIndex);
	}
	return codonCounts;
}


/* getRFPCountColumnName (NOT EXPOSED)
 * Arguments: None
 * Returns the vector of RFPCountColumnNames.
*/
std::vector <std::string> Genome::getRFPCountColumnNames()
{
	return RFPCountColumnNames;
}


/* addRFPCountColumnName (NOT EXPOSED)
 * Arguments: string categoryName
 * Adds the name of an RFP category to the RFPCountColumnNames vector.
*/
void Genome::addRFPCountColumnName(std::string categoryName)
{
	RFPCountColumnNames.push_back(categoryName);
}




//---------------------------------------//
//---------- Testing Functions ----------//
//---------------------------------------//


std::vector <unsigned> Genome::getNumGenesWithPhi()
{
	return numGenesWithPhi;
}


void Genome::setNumGenesWithPhi(std::vector <unsigned> newVector)
{
	this->numGenesWithPhi = newVector;
}





// -----------------------------------------------------------------------------------------------------//
// ---------------------------------------- R SECTION --------------------------------------------------//
// -----------------------------------------------------------------------------------------------------//



#ifndef STANDALONE


bool Genome::checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound)
{
	bool check = false;
	if (lowerbound <= index && index <= upperbound)
	{
		check = true;
	}
	else
	{
		my_printError("ERROR: Index % is out of bounds. Index must be between % & %\n", index, lowerbound, upperbound);
	}
	return check;
}

//NOTE: This function does the check and performs the function itself because of memory issues.
Gene& Genome::getGeneByIndex(unsigned index, bool simulated)
{
	if (simulated)
	{
		bool checker = checkIndex(index, 1, (unsigned)simulatedGenes.size());

		if (!checker)
		{
			my_printError("Warning: Invalid index given for simulated genes, returning simulated gene 1.\n");
		}

		return checker ? simulatedGenes[index - 1] : simulatedGenes[0];
	}
	else
	{
		bool checker = checkIndex(index, 1, (unsigned)genes.size());

		if (!checker)
		{
			my_printError("Warning: Invalid index given for genes, returning gene 1.\n");
		}

		return checker ? genes[index - 1] : genes[0];
	}
}


Gene& Genome::getGeneById(std::string ID, bool simulated)
{
	return getGene(ID, simulated);
}


Genome Genome::getGenomeForGeneIndicesR(std::vector <unsigned> indices, bool simulated)
{
	Genome genome;

	for (unsigned i = 0; i < indices.size(); i++)
	{
		if (indices[i] < 1 || indices[i] > getGenomeSize(simulated))
		{
			my_printError("Error in Genome::getGenomeForGeneIndices. An index specified is out of bounds for the genome!");
			my_printError("Returning empty Genome.");
			genome.clear();
			return genome;
		}
		else
		{
			simulated ? genome.addGene(simulatedGenes[indices[i]-1], true) : genome.addGene(genes[indices[i]-1], false);
		}
	}

	return genome;
}





//---------------------------------//
//---------- RCPP Module ----------//
//---------------------------------//


RCPP_EXPOSED_CLASS(Gene)
RCPP_EXPOSED_CLASS(Genome)

RCPP_MODULE(Genome_mod)
{
	class_<Genome>("Genome")
		//Constructors & Destructors:
		.constructor("empty constructor")


		//File I/O Functions:
		.method("readFasta", &Genome::readFasta, "reads a genome into the object")
		.method("writeFasta", &Genome::writeFasta, "writes the genome to a fasta file")
		.method("readRFPData", &Genome::readRFPData, "reads RFPData to be used in PA(NSE) models")
		.method("readSimRFPData", &Genome::readSimulatedGenomeFromPAModel, "reads already simulated RFPData to be used in PA(NSE) models")
		.method("writeRFPData", &Genome::writeRFPData, "writes RFPData used in PA(NSE) models")
		.method("readObservedPhiValues", &Genome::readObservedPhiValues)
		.method("removeUnobservedGenes", &Genome::removeUnobservedGenes)

		//Gene Functions:
		.method("addGene", &Genome::addGene) //TEST THAT ONLY!
		.method("getGenes", &Genome::getGenes) //TEST THAT ONLY!
		.method("getNumGenesWithPhi", &Genome::getNumGenesWithPhi) //TEST THAT ONLY!


		//Other Functions:
		.method("checkIndex", &Genome::checkIndex) //TEST THAT ONLY!
		.method("getGenomeSize", &Genome::getGenomeSize, "returns how many genes are in the genome")
		.method("clear", &Genome::clear, "clears the genome")
		.method("getCodonCountsPerGene", &Genome::getCodonCountsPerGene,
			"returns a vector of codon counts for a given gene")



		//R Section:
		.method("getGeneByIndex", &Genome::getGeneByIndex, "returns a gene for a given index")
		.method("getGeneById", &Genome::getGeneById) //TEST THAT ONLY!
		.method("getGenomeForGeneIndices", &Genome::getGenomeForGeneIndicesR,
			"returns a new genome based on the ones requested in the given vector")
		;
}
#endif
