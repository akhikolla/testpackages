/*
Author: Robert Hijmans
Date: April 2016
License: GPL (>=3)
*/


struct QueftsSoil {
	// Base (unfertilzied) soil supply of nutrients for standard crop with growth duration of 120 days
	double N_base_supply, P_base_supply, K_base_supply;
	// Recovery fractions of applied fertilizer :
	double N_recovery, P_recovery, K_recovery;
	// fraction uptake from soil supply as function of length of season, standard season is 120 days
	std::vector<double> UptakeAdjust;
};



struct QueftsCrop {
	// minimum and maximum concentration of NPK in vegetative organs and in storage organs (kg/kg)
    double NminStore, NminVeg, NmaxStore, NmaxVeg;
    double PminStore, PminVeg, PmaxStore, PmaxVeg;
    double KminStore, KminVeg, KmaxStore, KmaxVeg;
	
	// maximum amount of vegetative organs at zero yield of storage organs
	double Yzero;
	// fraction of crop's nitrogen uptake supplied by biological fixation
	double Nfix;
};



struct QueftsModel {
	
// INPUT
	QueftsSoil soil;
	QueftsCrop crop;

	// Crop biomass (water-limited production, by organ; dry-matter, kg/ha)
	double leaf_att, stem_att, store_att;
	double SeasonLength=120; // days	
	// fertilizer supplied
	double N_fertilizer=0, P_fertilizer=0, K_fertilizer=0;

	
// OUTPUT
	// nutrient supply from soil and fertilzer
	double N_supply, P_supply, K_supply;
	// nutrient uptake
	double UN, UP, UK;
	// nutrient limited yield of leaves, stems, and storage organ
	double leaf_lim, stem_lim, store_lim;
	// fertilizer required to reach attainable yield
	double N_gap, P_gap, K_gap;
	
// MODEL 
	void run();
	std::vector<double> output() {
		return std::vector<double> {N_gap, P_gap, K_gap, N_supply, P_supply, K_supply, leaf_lim, stem_lim, store_lim};
	}

	std::vector<double> runbatch(std::vector<double> Ns, std::vector<double> Ps, std::vector<double> Ks, std::vector<double> Ya);
		
};

