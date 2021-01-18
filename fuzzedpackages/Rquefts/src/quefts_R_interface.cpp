/* 
Rcpp interface to QUEFTS

Author: Robert Hijmans
Date: April 2016

License: GPL (>=3)
*/

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#include <vector>
#include "R_interface_util.h"
#include "quefts.h"


// [[Rcpp::export(name = ".quefts")]]
NumericVector quefts(List soil, List crop, List fertilizer) {
  
  // parameters
	QueftsModel q;
	q.soil.N_base_supply = doubleFromList(soil, "N_base_supply");
	q.soil.P_base_supply = doubleFromList(soil, "P_base_supply");
	q.soil.K_base_supply = doubleFromList(soil, "K_base_supply");
	q.soil.N_recovery = doubleFromList(soil, "N_recovery");
	q.soil.P_recovery = doubleFromList(soil, "P_recovery");
	q.soil.K_recovery = doubleFromList(soil, "K_recovery");
	q.soil.UptakeAdjust = vecFromList(soil, "UptakeAdjust");

	q.N_fertilizer = doubleFromList(fertilizer, "N");
	q.P_fertilizer = doubleFromList(fertilizer, "P");
	q.K_fertilizer = doubleFromList(fertilizer, "K");
	
	
	q.crop.NminStore = doubleFromList(crop, "NminStore");
	q.crop.NminVeg = doubleFromList(crop, "NminVeg");
	q.crop.NmaxStore = doubleFromList(crop, "NmaxStore");
	q.crop.NmaxVeg = doubleFromList(crop, "NmaxVeg");
	
	q.crop.PminStore = doubleFromList(crop, "PminStore");
	q.crop.PminVeg = doubleFromList(crop, "PminVeg");
	q.crop.PmaxStore = doubleFromList(crop, "PmaxStore");
	q.crop.PmaxVeg = doubleFromList(crop, "PmaxVeg");
	
	q.crop.KminStore = doubleFromList(crop, "KminStore");
	q.crop.KminVeg = doubleFromList(crop, "KminVeg");
	q.crop.KmaxStore = doubleFromList(crop, "KmaxStore");
	q.crop.KmaxVeg = doubleFromList(crop, "KmaxVeg");
	
	q.crop.Yzero = doubleFromList(crop, "Yzero");
	q.crop.Nfix = doubleFromList(crop, "Nfix");

	q.leaf_att = doubleFromList(crop, "leaf_att");
	q.stem_att = doubleFromList(crop, "stem_att");
	q.store_att = doubleFromList(crop, "store_att");
	q.SeasonLength = doubleFromList(crop, "SeasonLength");
	
// run model  
	q.run();

// prepare output	
	NumericVector result = NumericVector::create(q.N_supply, q.P_supply, q.K_supply, q.UN, q.UP, q.UK, q.leaf_lim, q.stem_lim, q.store_lim, q.N_gap, q.P_gap, q.K_gap);
	result.attr("names") = StringVector::create("N_supply", "P_supply", "K_supply", "N_uptake", "P_uptake", "K_uptake","leaf_yield", "stem_yield", "storage_yield", "N_gap", "P_gap", "K_gap");
	return(result);
}

