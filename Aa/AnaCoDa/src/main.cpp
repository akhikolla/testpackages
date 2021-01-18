#include "include/MCMCAlgorithm.h"
#include "include/Testing.h"
#include "include/SequenceSummary.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

//#define JEREMY
//#define STANDALONE

#ifdef CEDRIC
int main()
{
	std::cout << "Initializing MCMCAlgorithm object---------------" << std::endl;
	int samples = 2000;
	int thinning = 1000;
	int useSamples = 100;
	std::cout << "\t# Samples: " << samples << "\n";
	std::cout << "\tThinning: " << thinning << "\n";
	std::cout << "\t# Samples used: " << useSamples << "\n";
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thinning, 100, false, true, true);
	mcmc.setRestartFileSettings(std::string("test"), 100, true);

	//mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
	std::cout << "Done!-------------------------------\n\n\n";


	std::cout << "initialize Genome object--------------------------" << std::endl;
	bool withPhi = true;

	Genome genome;
	genome.readFasta("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/tests/testthat/UnitTestingData/testMCMCROCFiles/simulatedAllUniqueR.fasta");
	//genome.readFasta("F:/GitHub/RibModelDev/data/twoMixtures/simulatedAllUniqueR.fasta");
	//genome.readFasta("E:/RibosomeModel/RibModelDev/data/twoMixtures/simulatedAllUniqueR_unevenMixtures.fasta");
	if(withPhi)
	{
		genome.readObservedPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/tests/testthat/UnitTestingData/testMCMCROCFiles/simulatedAllUniqueR_phi_withPhiSet.csv", false);
		//genome.readObservedPhiValues("E:/RibosomeModel/RibModelDev/data/twoMixtures/simulatedAllUniqueR_phi_unevenMixtures.csv", false);
	}

	std::cout << "Done!-------------------------------\n\n\n";
	std::cout << "Initializing shared parameter variables---------------\n";

	std::cout << "Done!-------------------------------\n\n\n";
	std::cout << "Initializing shared parameter variables---------------\n";
	std::vector<unsigned> geneAssignment(genome.getGenomeSize());

	unsigned numMixtures = 1;
	std::vector<double> sphi_init(numMixtures, 1);

	/* For 2 mixture */
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		//geneAssignment[i] = ( ((double)rand() / (double)RAND_MAX) < 0.5 ? 0u : 1u );
		geneAssignment[i] = 0u;

	}
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
	std::cout << "Done!------------------------\n\n\n";

	std::cout << "initialize ROCParameter object" << std::endl;
	std::string mixDef = ROCParameter::allUnique;
	ROCParameter parameter(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

	for (unsigned i = 0u; i < numMixtures; i++)
	{
		unsigned selectionCategry = parameter.getSelectionCategory(i);
		std::cout << "Sphi_init for selection category " << selectionCategry << ": " << sphi_init[selectionCategry] << std::endl;
	}
	std::cout << "\t# mixtures: " << numMixtures << "\n";
	std::cout << "\tmixture definition: " << mixDef << "\n";

	std::vector<std::string> files(1);
	//files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_mutation0.csv");
	//files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_mutation1.csv");
	files[0] = std::string("/home/clandere/CodonUsageBias/R_roctest/bioinf_test_data/second_run/genome1/sim_oneMix_mutation_id_1_csp_-1_1_sphi_0.5.csv");
	//files[1] = std::string("/home/clandere/CodonUsageBias/RibosomeModel/RibModelDev/data/twoMixtures/simulated_mutation1.csv");
	parameter.initMutationCategories(files, parameter.getNumMutationCategories());
	files.resize(1);
	//files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection0.csv");
	//files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection1.csv");
	files[0] = std::string("/home/clandere/CodonUsageBias/R_roctest/bioinf_test_data/second_run/genome1/sim_oneMix_selection_id_1_csp_-1_1_sphi_0.5.csv");
	//files[1] = std::string("/home/clandere/CodonUsageBias/RibosomeModel/RibModelDev/data/twoMixtures/simulated_selection1.csv");
	parameter.initSelectionCategories(files, parameter.getNumSelectionCategories());

	parameter.InitializeSynthesisRate(genome, sphi_init[0]);
	//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
	//parameter.InitializeSynthesisRate(phiVals);
	std::cout << "done initialize ROCParameter object" << std::endl;


	std::cout << "Initializing ROCModel object\n";

	ROCModel model(withPhi);
	model.setParameter(parameter);


	std::cout << "starting MCMC for ROC" << std::endl;
	mcmc.run(genome, model, 1, 0);
	std::cout << std::endl << "Finished MCMC for ROC" << std::endl;


	std::cout << std::endl << "Exiting" << std::endl;
}

#endif // CEDRIC

#ifdef GABE
int main()
{


	if (1)
	{
		//testSequenceSummary();
		//testGene();
		testGenome("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/UnitTestingData");

		/*
		Genome genome;
		genome.readRFPData("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
		{
			geneAssignment[i] = 0u;
		}
		unsigned numMixtures = 1;
		std::vector<double> sphi_init(numMixtures, 2);
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		PAParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, "allUnique");

		std::vector<std::string> files;
		files.push_back("/Users/roxasoath1/Desktop/TONEWTON/RFPAlphaValues.csv");
		tmp.initMutationSelectionCategories(files, 1, PAParameter::alp);
		files[0] = "/Users/roxasoath1/Desktop/TONEWTON/RFPLambdaPrimeValues.csv";
		tmp.initMutationSelectionCategories(files, 1, PAParameter::lmPri);
		std::vector<double> phi = tmp.readPhiValues("/Users/roxasoath1/Desktop/TONEWTON/RFPPsiValues.csv");
		tmp.InitializeSynthesisRate(phi);


		PAModel model;

		model.setParameter(tmp);

		std::cout <<"init done\n";
		model.simulateGenome(genome);
		std::cout <<"writing file\n";
		genome.writeRFPData("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/rfp/simulatedRFPData.csv", true);
*/
		exit(1);
	}
	std::string modelToRun = "RFP"; //can also be ROC or FONSE
	bool withPhi = false;
	bool fromRestart = true;
	unsigned numMixtures = 1;


	std::cout << "Initializing MCMCAlgorithm object---------------" << std::endl;
	int samples = 10;
	int thinning = 10;
	int useSamples = 100;
	std::cout << "\t# Samples: " << samples << "\n";
	std::cout << "\tThinning: " << thinning << "\n";
	std::cout << "\t # Samples used: " << useSamples << "\n";
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thinning, 10, true, true, true);
	mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
	std::cout << "Done!-------------------------------\n\n\n";




	if (modelToRun == "ROC")
	{
		std::cout << "Initializing Genome object--------------------------" << std::endl;
		Genome genome;
		genome.readFasta("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/twoMixtures/simulatedAllUniqueR.fasta");
		if (withPhi)
		{
			genome.readObservedPhiValues("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulatedAllUniqueR_phi.csv", false);
		}
		std::cout << "Done!-------------------------------\n\n\n";



		std::cout << "Initializing shared parameter variables---------------\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout << "Done!------------------------\n\n\n";



		std::cout << "Initializing ROCParameter object--------------------\n" << std::endl;
		ROCParameter parameter;

		if (fromRestart)
		{
			ROCParameter tmp("/Users/roxasoath1/Desktop/RibModelFramework/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = ROCParameter::allUnique;
			ROCParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " << sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			std::vector<std::string> files(2);
			files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_mutation0.csv");
			files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_mutation1.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection0.csv");
			files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection1.csv");
			tmp.initSelectionCategories(files, tmp.getNumSelectionCategories());

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
		}
		std::cout << "Done!--------------------------------\n\n\n" << std::endl;



		std::cout << "Initializing ROCModel object--------------------------\n";

		ROCModel model;
		model.setParameter(parameter);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;


		std::cout << "Running MCMC.............\n" << std::endl;
		mcmc.run(genome, model, 1, 0);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;
	} //END OF ROC
	else if (modelToRun == "RFP")
	{
		std::cout << "Initializing Genome object--------------------------" << std::endl;
		Genome genome;
		genome.readRFPData("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
		std::cout << "Done!-------------------------------\n\n\n";



		std::cout << "Initializing shared parameter variables---------------\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout << "Done!------------------------\n\n\n";



		std::cout << "Initializing PAParameter object--------------------\n" << std::endl;
		PAParameter parameter;

		if (fromRestart)
		{
			PAParameter tmp("/Users/roxasoath1/Desktop/RibModelFramework/10_restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = Parameter::allUnique;
			PAParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++) {
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " <<
				sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			parameter = tmp;
		}
		std::cout << "Done!--------------------------------\n\n\n" << std::endl;



		std::cout << "Initializing PAModel object--------------------------\n";

		PAModel model;
		model.setParameter(parameter);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;


		std::cout << "Running MCMC.............\n" << std::endl;
		mcmc.run(genome, model, 1, 0);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;

	} //END OF RFP
	else if (modelToRun == "FONSE")
	{
		std::cout << "initialize Genome object--------------------------" << std::endl;
		Genome genome;
		genome.readFasta("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/FONSE/genome_2000.fasta");
		std::cout << "Done!-------------------------------\n\n\n";



		std::cout << "Initializing shared parameter variables---------------\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout << "Done!------------------------\n\n\n";



		FONSEParameter parameter;
		std::cout << "initialize Parameter object" << std::endl;
		if (fromRestart)
		{
			FONSEParameter tmp("/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = ROCParameter::allUnique;
			FONSEParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++) {
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " <<
				sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			std::vector<std::string> files(1);
			files[0] = std::string(
					"/Users/roxasoath1/Desktop/RibModelDevScripts/RibModelDev/data/FONSE/genome_2000.mutation.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
			parameter = tmp;
			std::cout << "done initialize Parameter object" << std::endl;
		}


		std::cout << "Initializing Model object\n";

		FONSEModel model;
		model.setParameter(parameter);


		std::cout << "starting MCMC for ROC" << std::endl;
		mcmc.run(genome, model, 4, 0);
		std::cout << std::endl << "Finished MCMC for ROC" << std::endl;

	}
}

#endif // GABE

#ifdef JEREMY
int main()
{
	unsigned index;
	bool fromRestart = false;
	std::string modelToRun = "PA";
	bool withPhi = false;


	std::cout << "Initializing MCMCAlgorithm object---------------" << std::endl;
	int samples = 50000;
	int thinning = 100;
	int useSamples = 100;
	unsigned numMixtures = 1;
	std::cout << "\t# Samples: " << samples << "\n";
	std::cout << "\tThinning: " << thinning << "\n";
	std::cout << "\t # Samples used: " << useSamples << "\n";
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thinning, 10, true, true, false);
	mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
	std::cout << "Done!-------------------------------\n\n\n";

	if (modelToRun == "FONSE") {
		std::cout << "initialize Genome object--------------------------" << std::endl;
		Genome genome;
		std::cout << "Reading fasta file\n";
		genome.readFasta("C:/Users/Alan/Documents/GitHub/RibModelDev/data/FONSE/nse2000.fasta");
		std::cout << "Done!-------------------------------\n\n\n";


		std::cout << "Initializing shared parameter variables---------------\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());

		std::vector<double> sphi_init(numMixtures, 1);

		/* For 1 mixture */
		for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
		{
			geneAssignment[i] = 0u;
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout << "Done!------------------------\n\n\n";

		//ROCParameter parameter;
		FONSEParameter parameter;
		std::cout << "initialize Parameter object" << std::endl;
		std::string mixDef = Parameter::allUnique;
		if (fromRestart)
		{
			FONSEParameter tmp("C:/Users/Alan/Documents/GitHub/RibModelDev/DevRScripts/10restartfile.rst");
			parameter = tmp;
		}
		else
		{
			//ROCParameter tmp(sphi_init, nu/mMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);
			FONSEParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " << sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			//std::vector<std::string> files(1);
			//files[0] = std::string("C:/Users/Jeremy/Documents/GitHub/RibModelDev/data/FONSE/Scereviciae.mut.csv");
			//tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			//tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			std::vector<double> phiVals = parameter.readPhiValues("C:/Users/Alan/Documents/GitHub/RibModelDev/data/FONSE/nse2000.phi.csv");
			tmp.InitializeSynthesisRate(phiVals);
			parameter = tmp;
		}
		std::cout << "done initialize Parameter object" << std::endl;


		std::cout << "Initializing Model object\n";

		bool withPhi = true;
		FONSEModel model;
		//ROCModel model(withPhi);
		model.setParameter(parameter);


		std::cout << "starting MCMC for FONSE" << std::endl;
		mcmc.run(genome, model, 1, 0);
		std::cout << std::endl << "Finished MCMC for FONSE" << std::endl;
	}
	else if (modelToRun == "RFP")
	{
		std::cout << "Initializing Genome object--------------------------" << std::endl;
		Genome genome;
		genome.readRFPData("/home/nax/Work/biolab/Logs/Input/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
		std::cout << "Done!-------------------------------\n\n\n";



		std::cout << "Initializing shared parameter variables---------------\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout << "Done!------------------------\n\n\n";



		std::cout << "Initializing PAParameter object--------------------\n" << std::endl;
		PAParameter parameter;

		if (fromRestart)
		{
			PAParameter tmp("/Users/roxasoath1/Desktop/RibModelFramework/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = Parameter::allUnique;
			PAParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++) {
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " <<
					sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			parameter = tmp;
		}
		std::cout << "Done!--------------------------------\n\n\n" << std::endl;



		std::cout << "Initializing PAModel object--------------------------\n";

		PAModel model;
		model.setParameter(parameter);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;


		std::cout << "Running MCMC.............\n" << std::endl;
		mcmc.run(genome, model, 1, 0);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;
	} //END OF RFP
	if (modelToRun == "ROC")
	{
		std::cout << "Initializing Genome object--------------------------" << std::endl;
		Genome genome;
		genome.readFasta("C:/Users/Alan/Documents/GitHub/RibModelDev/data/realGenomes/Skluyveri.fasta");
		if (withPhi)
		{
			genome.readObservedPhiValues("/Users/roxasoath1/Desktop/RibModelFramework/ribModel/data/simulatedAllUniqueR_phi.csv", false);
		}
		std::cout << "Done!-------------------------------\n\n\n";



		std::cout << "Initializing shared parameter variables---------------\n";
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		std::cout << "Done!------------------------\n\n\n";



		std::cout << "Initializing ROCParameter object--------------------\n" << std::endl;
		ROCParameter parameter;

		if (fromRestart)
		{
			ROCParameter tmp("/Users/roxasoath1/Desktop/RibModelFramework/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = ROCParameter::allUnique;
			ROCParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				std::cout << "Sphi_init for selection category " << selectionCategry << ": " << sphi_init[selectionCategry] << std::endl;
			}
			std::cout << "\t# mixtures: " << numMixtures << "\n";
			std::cout << "\tmixture definition: " << mixDef << "\n";

			/*std::vector<std::string> files(2);
			files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_mutation0.csv");
			files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_mutation1.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			files[0] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection0.csv");
			files[1] = std::string("F:/GitHub/RibModelDev/data/twoMixtures/simulated_selection1.csv");
			tmp.initSelectionCategories(files, tmp.getNumSelectionCategories());

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			*/
			//std::vector<double> phiVals = parameter.readPhiValues("/home/clandere/CodonUsageBias/RibosomeModel/RibModelFramework/ribModel/data/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
			parameter = tmp;
		}
		std::cout << "Done!--------------------------------\n\n\n" << std::endl;



		std::cout << "Initializing ROCModel object--------------------------\n";

		ROCModel model;
		model.setParameter(parameter);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;


		std::cout << "Running MCMC.............\n" << std::endl;
		std::cout << numMixtures << std::endl;
		mcmc.run(genome, model, 1, 0);
		std::cout << "Done!----------------------------------\n\n\n" << std::endl;
	} //END OF ROC
}

#endif // JEREMY

#ifdef HOLLIS
int main()
{
	std::string pathBegin = "/Users/hollisbui/";

	unsigned numMixtures = 1;
	std::vector<double> sphi_init(numMixtures, 2);
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;

	/*
	// SIMULATE GENOME: PA

	Genome genome;
	genome.readRFPData(pathBegin + "TODO");
	std::vector<unsigned> geneAssignment(genome.getGenomeSize());
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
		geneAssignment[i] = 0u;

	PAParameter parameter(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, "allUnique");

	std::vector<std::string> files;
	files.push_back(pathBegin + "miscGilchrist/runMe/HollisTestingData/RFPAlphaValues.csv");
	parameter.initMutationSelectionCategories(files, 1, PAParameter::alp);
	files[0] = pathBegin + "miscGilchrist/runMe/HollisTestingData/RFPLambdaPrimeValues.csv";
	parameter.initMutationSelectionCategories(files, 1, PAParameter::lmPri);

	std::vector<double> phi = parameter.readPhiValues(pathBegin + "miscGilchrist/runMe/HollisTestingData/RFPPhiValues.csv");
	//std::vector<double> phi = tmp.readPhiValues("/Users/roxasoath1/Desktop/TONEWTON/RFPPsiValues.csv");
	parameter.InitializeSynthesisRate(phi);

	PAModel model;

	model.setParameter(parameter);

	model.simulateGenome(genome);
	genome.writeRFPData(pathBegin + "miscGilchrist/runMe/HollisTestingOut/hbuiSimGenome5.30.17.csv", true);
	exit(1);
	 */

	// UNIT TESTING
	//testUtility();
	//testSequenceSummary();
	//testGene();
	//testGenome(pathBegin + "RibModelFramework/tests/testthat/UnitTestingData");
	//testCovarianceMatrix();
	//testParameter();
	//testMCMCAlgorithm();
	//exit(0);


	std::string modelToRun = "PA"; //can be PA, ROC or FONSE
	bool withPhi = false;
	bool fromRestart = false;


	my_print("Initializing MCMCAlgorithm object---------------\n");
	unsigned samples = 1000;
	unsigned thinning = 10;
	int useSamples = 100;
	my_print("\t# Samples: %\n", samples);
	my_print("\tThinning: %\n", thinning);
	my_print("\t # Samples used: %\n", useSamples);
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thinning, 10, true, true, true);
	//mcmc.setRestartFileSettings(pathBegin + "RestartFile.txt", 20, true);
	my_print("Done!-------------------------------\n\n\n");
	my_print("Initializing Genome object--------------------------\n");
	Genome genome;

	if (modelToRun == "ROC")
	{
		genome.readFasta(pathBegin + "RibModelDev/data/twoMixtures/simulatedAllUniqueR.fasta");
		if (withPhi)
		{
			genome.readObservedPhiValues(pathBegin + "RibModelFramework/ribModel/data/simulatedAllUniqueR_phi.csv", false);
		}
		my_print("Done!-------------------------------\n\n\n");


		my_print("Initializing shared parameter variables---------------\n");
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		my_print("Done!------------------------\n\n\n");


		my_print("Initializing ROCParameter object--------------------\n\n");
		ROCParameter parameter;

		if (fromRestart)
		{
			ROCParameter tmp(pathBegin + "RibModelFramework/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = ROCParameter::allUnique;
			ROCParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategory = tmp.getSelectionCategory(i);
				my_print("Sphi_init for selection category %: %\n", selectionCategory, sphi_init[selectionCategory]);
			}
			my_print("\t# mixtures: %\n", numMixtures);
			my_print("\tmixture definition: %\n", mixDef);

			std::vector<std::string> files(2);
			files[0] = std::string(pathBegin + "RibModelDev/data/twoMixtures/simulated_mutation0.csv");
			files[1] = std::string(pathBegin + "RibModelDev/data/twoMixtures/simulated_mutation1.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			files[0] = std::string(pathBegin + "RibModelDev/data/twoMixtures/simulated_selection0.csv");
			files[1] = std::string(pathBegin + "RibModelDev/data/twoMixtures/simulated_selection1.csv");
			tmp.initSelectionCategories(files, tmp.getNumSelectionCategories());

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			//std::vector<double> phiVals =
			// parameter.readPhiValues(pathBegin + "RibModelDev/data/realGenomes/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
			parameter = tmp;
		}
		my_print("Done!--------------------------------\n\n\n");


		my_print("Initializing ROCModel object--------------------------\n");
		ROCModel model;
		model.setParameter(parameter);
		my_print("Done!----------------------------------\n\n\n");


		my_print("Running MCMC.............\n\n");
		mcmc.run(genome, model, 1, 0);
		my_print("Done!----------------------------------\n\n\n");
	} //END OF ROC
	else if (modelToRun == "PA")
	{
		genome.readRFPData(pathBegin + "labnotebooks/Hollis.Bui/RunnableScriptsAndData/June2017PAModelPop/PopPADataOneGene.csv");
		my_print("Done!-------------------------------\n\n\n");


		my_print("Initializing shared parameter variables---------------\n");
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		my_print("Done!------------------------\n\n\n");


		my_print("Initializing PAParameter object--------------------\n\n");
		PAParameter parameter;
		parameter.writeBasicRestartFile(pathBegin + "RibModelFramework/HollisRestartFile.txt");

		if (fromRestart)
		{
			PAParameter tmp(pathBegin + "RibModelFramework/10_restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = Parameter::allUnique;
			PAParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategory = tmp.getSelectionCategory(i);
				my_print("Sphi_init for selection category %: %\n", selectionCategory, sphi_init[selectionCategory]);
			}
			my_print("\t# mixtures: %\n", numMixtures);
			my_print("\tmixture definition: %\n", mixDef);

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			parameter = tmp;
			parameter.writeEntireRestartFile(pathBegin + "RibModelFramework/HollisRestartFile2.txt");
		}
		my_print("Done!--------------------------------\n\n\n");


		my_print("Initializing PAModel object--------------------------\n");
		PAModel model;
		model.setParameter(parameter);
		my_print("Done!----------------------------------\n\n\n");


		my_print("Running MCMC.............\n\n");
		mcmc.run(genome, model, 1, 0);
		my_print("Done!----------------------------------\n\n\n");

		std::vector<double> alphaList (61,0);
		std::vector<double> lambdaPrimeList (61,0);
		std::vector<double> waitingTimes (61,0);
		std::vector<double> waitRates (61,0);
		std::vector<std::string> codonList = SequenceSummary::codons();
		unsigned cat = 0u;

		for (unsigned i = 0u; i < 61; i++)
		{
			std::string codon = codonList[i];

			alphaList[i] = parameter.getCodonSpecificPosteriorMean(cat, samples / 2, codon, 0, false);
			//alphaTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 0, FALSE)
			//alpha.ci[i,] <- quantile(alphaTrace[(samples * 0.5):samples], probs = c(0.025,0.975))

			lambdaPrimeList[i] = parameter.getCodonSpecificPosteriorMean(cat, samples / 2, codon, 1, false);
			//lambdaPrimeTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 1, FALSE)
			//lambdaPrime.ci[i,] <- quantile(lambdaPrimeTrace[(samples * 0.5):samples], probs = c(0.025,0.975))

			waitingTimes[i] = alphaList[i] / lambdaPrimeList[i];
			waitRates[i] = (1.0/waitingTimes[i]);

			my_print("For codon %, alpha is %, lambda prime is %, and waitingTime is %.\n", codonList[i], alphaList[i],
					 lambdaPrimeList[i], waitingTimes[i]);
		}
	} //END OF PA
	else if (modelToRun == "FONSE")
	{
		genome.readFasta(pathBegin + "RibModelDev/data/singleMixture/genome_2000.fasta");
		my_print("Done!-------------------------------\n\n\n");


		my_print("Initializing shared parameter variables---------------\n");
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		my_print("Done!------------------------\n\n\n");


		FONSEParameter parameter;
		my_print("initialize Parameter object\n");
		if (fromRestart)
		{
			FONSEParameter tmp(pathBegin + "RibModelDev/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = ROCParameter::allUnique;
			FONSEParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategory = tmp.getSelectionCategory(i);
				my_print("Sphi_init for selection category %: %\n", selectionCategory, sphi_init[selectionCategory]);
			}
			my_print("\t# mixtures: %\n", numMixtures);
			my_print("\tmixture definition: %\n", mixDef);

			std::vector<std::string> files(1);
			files[0] = std::string(pathBegin + "RibModelDev/data/singleMixture/genome_2000.mutation.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			//std::vector<double> phiVals = parameter.readPhiValues(
			// pathBegin + "RibModelDev/data/realGenomes/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
			parameter = tmp;
			my_print("done initialize Parameter object\n");
		}


		my_print("Initializing Model object\n");
		FONSEModel model;
		model.setParameter(parameter);
		my_print("Done!------------------------\n\n\n");


		my_print("starting MCMC for ROC\n");
		mcmc.run(genome, model, 4, 0);
		my_print("\nFinished MCMC for ROC\n");

	}
}

#endif // Hollis

#ifdef DENIZHAN
int main()
{
    std::vector <double> alphas;
    std::vector <double> lambdas;
    std::vector <std::string> cspFiles;
    std::string pathBegin = "/home/nax/Work/biolab/TestingIn/";
    Genome genome;
	unsigned numMixtures = 1;
	std::vector<double> sphi_init(numMixtures, 2);
	std::vector<std::vector<unsigned> > mixtureDefinitionMatrix;
	std::vector<unsigned> geneAssignment;
	genome.readRFPData("/home/nax/Work/biolab/Logs/Main/simRFP2.csv", false);
    geneAssignment.resize(genome.getGenomeSize());
    my_print("%\n", genome.getGenomeSize());
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		geneAssignment[i] = 0u;
	}

    PANSEParameter parameter(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, "allUnique");

    cspFiles.push_back("/home/nax/Work/biolab/Dev/data/rfp/RFPAlphaValues.csv");
    //parameter.readAlphaValues(cspFiles[0]);
    parameter.initMutationSelectionCategories(cspFiles, 1, parameter.alp);

    cspFiles[0] = ("/home/nax/Work/biolab/Dev/data/rfp/RFPLambdaPrimeValues.csv");
    //parameter.readLambdaValues(cspFiles[1]);
    parameter.initMutationSelectionCategories(cspFiles, 1, parameter.lmPri);

    alphas = parameter.oneMixAlpha();
    lambdas = parameter.oneMixLambda();

    for(int i = 0; i < alphas.size(); i++){
        my_print("%,%\n", SequenceSummary::indexToCodon(i), alphas[i]);
    }
    exit(0);
/*

	unsigned numMixtures = 1;
	std::vector<double> sphi_init(numMixtures, 2);
	std::vector<std::vector<unsigned> > mixtureDefinitionMatrix;

	// SIMULATE GENOME: RFP
	Genome genome;
    //if(testEqualityGenome(genome, genome)){
      //  my_print("So far so good\n");
    //}
    //exit(1);

	genome.readRFPData(pathBegin + "rfp_file_20positions_20genes.csv", false);
    exit(0);
	genome.readFasta(pathBegin + "RibModelDev/data/singleMixture/genome_2000.fasta", false);

	std::vector<unsigned> geneAssignment(genome.getGenomeSize());
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		geneAssignment[i] = 0u;
	}

	PAParameter parameter(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, "allUnique");
	PAModel model;

    //ROCParameter parameter(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, "allUnique");
	//ROCModel model;

	model.setParameter(parameter);

	model.simulateGenome(genome);
	genome.writeRFPData(pathBegin + "labbooks/Denizhan.Pak/Log_Files/sim_genomes/PASim.csv", true);
	genome.writeRFPData(pathBegin + "labbooks/Denizhan.Pak/Log_Files/sim_genomes/PANotSim.csv", false);
	exit(1);*/


	// UNIT TESTING
	//testUtility();
	//testSequenceSummary();
	//testGene();
	//testGenome(pathBegin + "RibModelFramework/tests/testthat/UnitTestingData");
	//testCovarianceMatrix();
	//testParameter();
	//testParameterWithFile(pathBegin + "HollisFile.txt");
	//testPAParameter();
	//testMCMCAlgorithm();
	//exit(0);


	std::string modelToRun = "PA"; //can be RFP, ROC or FONSE
	bool withPhi = false;
    getStdCspForIndex(unsigned i);
	bool fromRestart = false;


	my_print("Initializing MCMCAlgorithm object---------------\n");
	unsigned samples = 100;
	unsigned thinning = 100;
	int useSamples = 1000;
	my_print("\t# Samples: %\n", samples);
	my_print("\tThinning: %\n", thinning);
	my_print("\t # Samples used: %\n", useSamples);
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thinning, 10, true, true, true);
	mcmc.setRestartFileSettings(pathBegin + "RestartFile.txt", 20, true);
	my_print("Done!-------------------------------\n\n\n");


	if (modelToRun == "ROC")
	{
		my_print("Initializing Genome object--------------------------\n");
		Genome genome;
		genome.readFasta(pathBegin + "RibModelDev/data/twoMixtures/simulatedAllUniqueR.fasta");
		if (withPhi)
		{
			genome.readObservedPhiValues(pathBegin + "RibModelFramework/ribModel/data/simulatedAllUniqueR_phi.csv", false);
		}
		my_print("Done!-------------------------------\n\n\n");


		my_print("Initializing shared parameter variables---------------\n");
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		my_print("Done!------------------------\n\n\n");


		my_print("Initializing ROCParameter object--------------------\n\n");
		ROCParameter parameter;

		if (fromRestart)
		{
			ROCParameter tmp(pathBegin + "RibModelFramework/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = ROCParameter::allUnique;
			ROCParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				my_print("Sphi_init for selection category %: %\n", selectionCategry, sphi_init[selectionCategry]);
			}
			my_print("\t# mixtures: %\n", numMixtures);
			my_print("\tmixture definition: %\n", mixDef);

			std::vector<std::string> files(2);
			files[0] = std::string(pathBegin + "RibModelDev/data/twoMixtures/simulated_mutation0.csv");
			files[1] = std::string(pathBegin + "RibModelDev/data/twoMixtures/simulated_mutation1.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			files[0] = std::string(pathBegin + "RibModelDev/data/twoMixtures/simulated_selection0.csv");
			files[1] = std::string(pathBegin + "RibModelDev/data/twoMixtures/simulated_selection1.csv");
			tmp.initSelectionCategories(files, tmp.getNumSelectionCategories());

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			//std::vector<double> phiVals =
			// parameter.readPhiValues(pathBegin + "RibModelDev/data/realGenomes/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
			parameter = tmp;
		}
		my_print("Done!--------------------------------\n\n\n");


		my_print("Initializing ROCModel object--------------------------\n");
		ROCModel model;
		model.setParameter(parameter);
		my_print("Done!----------------------------------\n\n\n");


		my_print("Running MCMC.............\n\n");
		mcmc.run(genome, model, 1, 0);
		my_print("Done!----------------------------------\n\n\n");
	} //END OF ROC
	else if (modelToRun == "PANSE")
	{

        my_print("Initializing Genome object--------------------------\n");
		Genome genome;
		//genome.readRFPData(pathBegin + "RibModelDev/data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv");
		genome.readRFPData(pathBegin + "PopPAData.csv");
		my_print("Done!-------------------------------\n\n\n");


		my_print("Initializing shared parameter variables---------------\n");
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		my_print("Done!------------------------\n\n\n");


		my_print("Initializing PANSEParameter object--------------------\n\n");
		//PANSEParameter parameter;
		//parameter.writeBasicRestartFile("/Users/hollisbui/HollisFile.txt");

		if (fromRestart)
		{
			PANSEParameter tmp(pathBegin + "RibModelFramework/10_restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = Parameter::allUnique;
			PANSEParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategry = tmp.getSelectionCategory(i);
				my_print("Sphi_init for selection category %: %\n", selectionCategry, sphi_init[selectionCategry]);
			}
			my_print("\t# mixtures: %\n", numMixtures);
			my_print("\tmixture definition: %\n", mixDef);

			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			parameter = tmp;
			//parameter.writeEntireRestartFile("/Users/hollisbui/HollisFile2.txt");
		}
		my_print("Done!--------------------------------\n\n\n");


		my_print("Initializing PAModel object--------------------------\n");
		PANSEModel model;
		model.setParameter(parameter);
		my_print("Done!----------------------------------\n\n\n");


		my_print("Running MCMC.............\n\n");
		mcmc.run(genome, model, 1, 0);
		my_print("Done!----------------------------------\n\n\n");

	} //END OF PANSE
	else if (modelToRun == "FONSE")
	{
		my_print("initialize Genome object--------------------------\n");
		Genome genome;
		genome.readFasta(pathBegin + "RibModelDev/data/singleMixture/genome_2000.fasta");
		my_print("Done!-------------------------------\n\n\n");


		my_print("Initializing shared parameter variables---------------\n");
		std::vector<unsigned> geneAssignment(genome.getGenomeSize());
		std::vector<double> sphi_init(numMixtures, 1);

		if (numMixtures == 1)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				geneAssignment[i] = 0u;
			}
		}
		else if (numMixtures == 3)
		{
			for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
			{
				if (i < 961) geneAssignment[i] = 0u;
				else if (i < 1418) geneAssignment[i] = 1u;
				else geneAssignment[i] = 0u;
			}
		}
		std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
		my_print("Done!------------------------\n\n\n");


		FONSEParameter parameter;
		my_print("initialize Parameter object\n");
		if (fromRestart)
		{
			FONSEParameter tmp(pathBegin + "RibModelDev/DevRscripts/10restartFile.rst");
			parameter = tmp;
		}
		else
		{
			std::string mixDef = ROCParameter::allUnique;
			FONSEParameter tmp(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);

			for (unsigned i = 0u; i < numMixtures; i++)
			{
				unsigned selectionCategory = tmp.getSelectionCategory(i);
				my_print("Sphi_init for selection category %: %\n", selectionCategory, sphi_init[selectionCategory]);
			}
			my_print("\t# mixtures: %\n", numMixtures);
			my_print("\tmixture definition: %\n", mixDef);

			std::vector<std::string> files(1);
			files[0] = std::string(pathBegin + "RibModelDev/data/singleMixture/genome_2000.mutation.csv");
			tmp.initMutationCategories(files, tmp.getNumMutationCategories());
			tmp.InitializeSynthesisRate(genome, sphi_init[0]);
			//std::vector<double> phiVals = parameter.readPhiValues(
			// pathBegin + "RibModelDev/data/realGenomes/Skluyveri_ChrA_ChrCleft_phi_est.csv");
			//parameter.InitializeSynthesisRate(phiVals);
			parameter = tmp;
			my_print("done initialize Parameter object\n");
		}


		my_print("Initializing Model object\n");
		FONSEModel model;
		model.setParameter(parameter);
		my_print("Done!------------------------\n\n\n");


		my_print("starting MCMC for ROC\n");
		mcmc.run(genome, model, 4, 0);
		my_print("\nFinished MCMC for ROC\n");

	}
}

#endif // Denizhan

#ifdef ALEX
int main()
{
	srand(1500);
    std::vector <double> alphas;
    std::vector <double> lambdas;
    std::vector <std::string> cspFiles;
    std::string pathBegin = "/home/acope3/Panse_project/Logs/Input";
    Genome genome;
	unsigned numMixtures = 1;
	std::vector<double> sphi_init(numMixtures, 1);
	std::vector<unsigned> geneAssignment;
	genome.readRFPData("/home/acope3/Panse_project/Logs/Input/Stochastic_Simulation/Alex_simulations_10_24_remove_stops.csv", false);
    geneAssignment.resize(genome.getGenomeSize());
    my_print("%\n", genome.getGenomeSize());
	for (unsigned i = 0u; i < genome.getGenomeSize(); i++)
	{
		geneAssignment[i] = 0u;
	}

	std::vector<double> phi;
	std::size_t pos;
	std::ifstream currentFile;
	std::string tmpString;
	my_print("Initializing gene expression...\n");
	currentFile.open("/home/acope3/Panse_project/Input/PopData/orderedRandGeneIDPhiMean.csv");
	currentFile >> tmpString;
	while (currentFile >> tmpString)
	{
		pos = tmpString.find(',');
		if (pos != std::string::npos)
		{
			std::string val = tmpString.substr(pos + 1, std::string::npos);
			phi.push_back(std::atof(val.c_str()));
		}
	}
	my_print("Initializing CSP\n");
	std::vector<std::vector<unsigned>> mixtureDefinitionMatrix;
	std::string mixDef = PANSEParameter::allUnique;
	PANSEParameter parameter(sphi_init, numMixtures, geneAssignment, mixtureDefinitionMatrix, true, mixDef);
    cspFiles.push_back("/home/acope3/Panse_project/Input/PopData/JeremyRFPAlphaValues.csv");
    parameter.initMutationSelectionCategories(cspFiles, 1, parameter.alp);
    cspFiles[0] = ("/home/acope3/Panse_project/Input/PopData/JeremyRFPLambdaPrimeValues.csv");
    parameter.initMutationSelectionCategories(cspFiles, 1, parameter.lmPri);
    cspFiles[0] = ("/home/acope3/Panse_project/Logs/Input/Stochastic_Simulation/simNSEMay_4.csv");
    parameter.initMutationSelectionCategories(cspFiles, 1, parameter.nse);
    parameter.InitializeSynthesisRate(phi);
	PANSEModel model;
	model.setParameter(parameter);

	my_print("Initializing MCMCAlgorithm object---------------\n");
	unsigned samples = 50;
	unsigned thinning = 2;

	my_print("\t# Samples: %\n", samples);
	my_print("\tThinning: %\n", thinning);
	MCMCAlgorithm mcmc = MCMCAlgorithm(samples, thinning, 10, false, true, false);
	mcmc.setRestartFileSettings("RestartFile.txt", 20, true);
	my_print("Done!-------------------------------\n\n\n");

	my_print("Done!--------------------------------\n\n\n");


	my_print("Initializing PAModel object--------------------------\n");
	model.setParameter(parameter);
	my_print("Done!----------------------------------\n\n\n");


	my_print("Running MCMC.............\n\n");
	mcmc.run(genome, model, 1, 0);
	my_print("Done!----------------------------------\n\n\n");

	std::vector<std::string> codons = parameter.getGroupList();
	std::vector<double> alpha,lmprime,nse;
	double tmp;
	for (unsigned i=0; i < codons.size();i++)
	{
		tmp = parameter.getCodonSpecificPosteriorMean(0, samples, codons[i],0, true, false);
		alpha.push_back(tmp);
		tmp = parameter.getCodonSpecificPosteriorMean(0, samples, codons[i],1, true, false);
		lmprime.push_back(tmp);
		tmp = parameter.getCodonSpecificPosteriorMean(0, samples, codons[i],2, true, false);
		nse.push_back(tmp);
	}
	std::ofstream myFile("cpp_runs_10_24.csv");
	myFile << "Codon,Alpha,LambdaPrime,NSE\n";
	for (unsigned i=0;i < codons.size();i++)
	{
		myFile << codons[i] <<","<<alpha[i]<<","<<lmprime[i]<<","<<nse[i]<<"\n";
	}
	std::cout << mcmc.getLogPosteriorMean(samples);


}

#endif //Alex

