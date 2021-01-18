/* Calculate nutrient requirements and nutrient limited yields according to the QUEFTS model 
After Janssen, B.H., Guiking, F.C.T., Van der Eijk, D., Smaling, E.M.A., Wolf, J., Van Reuler,
H., 1990. A system for quantitative evaluation of the fertility of tropical soils
(QUEFTS). Geoderma 46, 299–318.

Robert Hijmans
April 2016
*/

using namespace std;
#include <vector>
#include <algorithm>
#include <cmath>
#include "SimUtil.h"
#include "quefts.h"


// QUEFTS functions
double uptake(double Si, double iZero, double iYratA, double iYratD, double Sj, double jZero, double jYratA, double jYratD) {
	double up;
	if (Si <= (iZero + jYratA * (Sj - jZero) / iYratD))  {
		up = Si;
	} else if (Si < (iZero + (Sj - jZero) * (2 * jYratD / iYratA - jYratA / iYratD)))  { 
		up = Si -  (0.25 * pow(Si - iZero - jYratA * (Sj - jZero) / iYratD, 2)) / ((jYratD / iYratA - jYratA / iYratD) * (Sj - jZero));
	} else {
		up = iZero + jYratD * (Sj - jZero) / iYratA;
	} 
	return(up);
}


double yield(double iYA, double jYA, double iUPT, double izero, double iYratD, double iYratA, double iYD, double jYD, double zYD, double YldO) {
	double EY;
	double LYD = std::min({iYD, jYD, zYD, YldO});
	
	if (iYD <= jYA)  {
		EY = iYD;
		if (iYD > LYD) {
			EY = LYD;
		}
		return(EY);
	} else if ((iYA >= LYD) | (iYD <= LYD))  {
		if ((iYA < LYD) & (iYD == LYD)) {
			LYD = std::min({jYD, zYD, YldO});
		} else {
			EY = LYD;
		return(EY);
		} // end if;
	} // end if;
	EY = jYA + 2*(LYD-jYA) * (iUPT-izero - jYA/iYratD) / (LYD/iYratA - jYA / iYratD) - (LYD - jYA) * pow(iUPT-izero - jYA/iYratD, 2) / pow(LYD/iYratA - jYA/iYratD, 2);
	if (jYA > LYD) {
		EY = LYD;
	}
	return(EY);
}


// nutrient requirements at balanced concentrations
std::vector<double> requirements(double minVeg, double maxVeg, double minStore, double maxStore, double supply, double recovery, double leaf_att, double stem_att, double store_att, double fix) {

// fertilizer requirements
	double limVeg = (minVeg + maxVeg) / 2;
    double limStore = (maxStore + minStore) / 2;
    double req = fix ? 0 : (limVeg * (leaf_att + stem_att) + limStore * store_att);
    double reqFert = std::max(0.0, (req - supply) / recovery);

// nutrient requirements at accumulated and diluted nutrient concentration
    double reqAcc = (maxVeg * (leaf_att + stem_att) + maxStore * store_att) * (1 - fix);
    double reqDil = (minVeg * (leaf_att + stem_att) + minStore * store_att) * (1 - fix);

	std::vector<double> out = {reqAcc, reqDil, reqFert};
	return(out);
}


void QueftsModel::run() {

	double Nzero = 0, Pzero = 0, Kzero = 0;

	// which plant organs to consider
    int organs;
	double Yo;
    if (store_att > 200)  {
        Yo = store_att;
        organs = 3;
		Nzero = crop.Yzero * (1 - crop.Nfix) * (crop.NmaxVeg + 2 * crop.NminVeg) / 3;
		Pzero = crop.Yzero * (crop.PmaxVeg + 2 * crop.PminVeg) / 3;
		Kzero = crop.Yzero * (crop.KmaxVeg + 2 * crop.KminVeg) / 3;
    } else if (stem_att <= 200)  {
        Yo = std::max(5.0, leaf_att);
        organs = 1;
    } else {
        Yo = stem_att;
        organs = 2;
    } 

// nutrient supply of soil during growing season 
	double relgrowseason = approx(soil.UptakeAdjust, SeasonLength);
    N_supply = soil.N_base_supply * relgrowseason + N_fertilizer * soil.N_recovery;
    P_supply = soil.P_base_supply * relgrowseason + P_fertilizer * soil.P_recovery;
    K_supply = soil.K_base_supply * relgrowseason + K_fertilizer * soil.K_recovery;
  
	std::vector<double> Nreq = requirements(crop.NminVeg, crop.NmaxVeg, crop.NminStore, crop.NmaxStore, N_supply, soil.N_recovery, leaf_att, stem_att, store_att, crop.Nfix);
	std::vector<double> Preq = requirements(crop.PminVeg, crop.PmaxVeg, crop.PminStore, crop.PmaxStore, P_supply, soil.P_recovery, leaf_att, stem_att, store_att, 0.);
	std::vector<double> Kreq = requirements(crop.KminVeg, crop.KmaxVeg, crop.KminStore, crop.KmaxStore, K_supply, soil.K_recovery, leaf_att, stem_att, store_att, 0.);
	N_gap = Nreq[2];
	P_gap = Preq[2];
	K_gap = Kreq[2];

// yield/nutrient uptake ratios at accumulated and diluted N, P, and K concentrations in the plant parts
    double NYratA = Yo / std::max(0.1, Nreq[0] - Nzero);
    double NYratD = Yo / std::max(0.1, Nreq[1] - Nzero);
    double PYratA = Yo / std::max(0.1, Preq[0] - Pzero);
    double PYratD = Yo / std::max(0.1, Preq[1] - Pzero);
    double KYratA = Yo / std::max(0.1, Kreq[0] - Kzero);
    double KYratD = Yo / std::max(0.1, Kreq[1] - Kzero);


// N uptake
    double UNP = uptake(N_supply, Nzero, NYratA, NYratD, P_supply, Pzero, PYratA, PYratD);
    double UNK = uptake(N_supply, Nzero, NYratA, NYratD, K_supply, Kzero, KYratA, KYratD);
    UN  = std::max(0.0, std::min(UNP, UNK));

// P uptake
    double UPN = uptake(P_supply, Pzero, PYratA, PYratD, N_supply, Nzero, NYratA, NYratD);
    double UPK = uptake(P_supply, Pzero, PYratA, PYratD, K_supply, Kzero, KYratA, KYratD);
    UP  = std::max(0.0, std::min(UPN, UPK));
	
// K uptake
    double UKN = uptake(K_supply, Kzero, KYratA, KYratD, N_supply, Nzero, NYratA, NYratD);
    double UKP = uptake(K_supply, Kzero, KYratA, KYratD, P_supply, Pzero, PYratA, PYratD);
    UK  = std::max(0.0, std::min(UKN, UKP));
	
//  N, P, and K limited yields at accumulated and diluted nutrient concentration
    double YNA = NYratA * std::max(0.0, UN - Nzero);
    double YND = NYratD * std::max(0.0, UN - Nzero);
    double YPA = PYratA * std::max(0.0, UP - Pzero);
    double YPD = PYratD * std::max(0.0, UP - Pzero);
    double YKA = KYratA * std::max(0.0, UK - Kzero);
    double YKD = KYratD * std::max(0.0, UK - Kzero);

// N and P limited yields
	double EYNP = yield(YNA, YPA, UN, Nzero, NYratD, NYratA, YND, YPD, YKD, Yo);
	double EYPN = yield(YPA, YNA, UP, Pzero, PYratD, PYratA, YPD, YND, YKD, Yo);

// N and K limited yields
	double EYNK = yield(YNA, YKA, UN, Nzero, NYratD, NYratA, YND, YKD, YPD, Yo);
	double EYKN = yield(YKA, YNA, UK, Kzero, KYratD, KYratA, YKD, YND, YPD, Yo);

// P and K limited yields
	double EYPK = yield(YPA, YKA, UP, Pzero, PYratD, PYratA, YPD, YKD, YND, Yo);
	double EYKP = yield(YKA, YPA, UK, Kzero, KYratD, KYratA, YKD, YPD, YND, Yo);  

// N, P, and K limited yields
    double limYield = std::min({YND, YPD, YKD, Yo});
    double yield = std::min(limYield, (EYNP + EYPN + EYNK + EYKN + EYPK + EYKP) / 6);
	
	  
//  nutrient-limited yield
    if (organs == 3)  {
        store_lim = yield;
        leaf_lim = leaf_att / (leaf_att + stem_att) * 
				((stem_att + leaf_att - crop.Yzero) * yield / store_att + crop.Yzero);
        if (yield <= 0) {
			leaf_lim = leaf_att / (leaf_att + stem_att) * 
					std::min({
						crop.Yzero,	
						UN / (0.333 * (crop.NminVeg + 2 * crop.NminVeg) * (1. - crop.Nfix)),
						UP / (0.333 * (crop.PminVeg + 2 * crop.PminVeg)), 
						UK / (0.333 * (crop.KminVeg + 2 * crop.KminVeg))});
		}
		stem_lim = (stem_att / leaf_att) * leaf_lim;
	} else if (organs != 2) {
		leaf_lim = yield;
		stem_lim = (stem_att / leaf_att) * leaf_lim;
		store_lim = (store_att / leaf_att) * leaf_lim;
	} else {
		stem_lim = yield;
		leaf_lim = (leaf_att / stem_att) * stem_lim;
		store_lim = (store_att / stem_att) * stem_lim;
	}  
	
}



std::vector<double> QueftsModel::runbatch(std::vector<double> Ns, std::vector<double> Ps, std::vector<double> Ks, std::vector<double> Ya) {
	
	size_t n = Ns.size();
	std::vector<double> out(n, NAN); 
	for (size_t i=0; i<n; i++) {
		if (isnan(Ns[i])) continue;	
		soil.N_base_supply = Ns[i];
		soil.P_base_supply = Ps[i];
		soil.K_base_supply = Ks[i];
		store_att = Ya[i];
		leaf_att = 0.45 * store_att;
		stem_att = 0.55 * store_att;
		run();
		out[i] = store_lim;
	}
	return out;
}	


