#ifndef GENOME_H
#define GENOME_H


#include "Gene.h"


#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <cmath>
#include <functional>

#ifndef STANDALONE
#include <Rcpp.h>
#endif

class Model;
class Genome
{
	private:

		std::vector <Gene> genes;
		std::vector <Gene> simulatedGenes;
		std::vector <unsigned> numGenesWithPhi; //Number of phi sets is vector size, value is number of genes
												//with a phi value for that set. Values should currently be equal.
        std::vector<std::string> RFPCountColumnNames;
        unsigned prev_genome_size;
        unsigned totalRFPCount;


  	public:

	public:

		//Constructors & Destructors:
		explicit Genome();
		Genome& operator=(const Genome& other);
		bool operator==(const Genome& other) const;
		virtual ~Genome();


		//File I/O Functions:
		void readFasta(std::string filename, bool append = false);
		void writeFasta(std::string filename, bool simulated = false);
		void readRFPData(std::string filename, bool append = false, bool positional = false);
		void writeRFPData(std::string filename, bool simulated = false);
		void readObservedPhiValues(std::string filename, bool byId = true);
		void removeUnobservedGenes();
        void readSimulatedGenomeFromPAModel(std::string filename);


		//Gene Functions:
		void addGene(Gene& gene, bool simulated = false);
		std::vector <Gene> getGenes(bool simulated = false);
		unsigned getNumGenesWithPhiForIndex(unsigned index);
		Gene& getGene(unsigned index, bool simulated = false);
		Gene& getGene(std::string id, bool simulated = false);


		//Other Functions:
		unsigned getGenomeSize(bool simulated = false);
		void clear();
		Genome getGenomeForGeneIndices(std::vector <unsigned> indices, bool simulated = false);
		std::vector <unsigned> getCodonCountsPerGene(std::string codon);
        std::vector <std::string> getRFPCountColumnNames();
		void addRFPCountColumnName(std::string categoryName);
		unsigned getSumRFP();

		//Testing Functions:
		std::vector <unsigned> getNumGenesWithPhi();
		void setNumGenesWithPhi(std::vector <unsigned> newVector);


		//R Section:

#ifndef STANDALONE

		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);
		Gene& getGeneByIndex(unsigned index, bool simulated = false);
		Gene& getGeneById(std::string ID, bool simulated = false);
		Genome getGenomeForGeneIndicesR(std::vector <unsigned> indices, bool simulated = false);

#endif //STANDALONE

	protected:
};

#endif // GENOME_H
