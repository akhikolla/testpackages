#include <Rcpp.h>
#include "quefts.h"

using namespace Rcpp;


NumericVector runout(QueftsModel* q) {
	q->run();
	NumericVector out = NumericVector::create(
		_["leaf_lim"] = q->leaf_lim,
		_["stem_lim"] = q->stem_lim,
		_["store_lim"] = q->store_lim,
		_["N_supply"] = q->N_supply,
		_["P_supply"] = q->P_supply,
		_["K_supply"] = q->K_supply,
		_["N_uptake"] = q->UN,
		_["P_uptake"] = q->UP,
		_["K_uptake"] = q->UK,
		_["N_gap"] = q->N_gap,
		_["P_gap"] = q->P_gap,
		_["K_gap"] = q->K_gap
	);
	return(out);
}


RCPP_EXPOSED_CLASS(QueftsModel)
RCPP_EXPOSED_CLASS(QueftsCrop)
RCPP_EXPOSED_CLASS(QueftsSoil)
	
RCPP_MODULE(QUEFTS){
    using namespace Rcpp;

    class_<QueftsSoil>("QueftsSoil")
		.field("N_base_supply", &QueftsSoil::N_base_supply, "N_base_supply")
		.field("P_base_supply", &QueftsSoil::P_base_supply, "P_base_supply")
		.field("K_base_supply", &QueftsSoil::K_base_supply, "K_base_supply")
		.field("N_recovery", &QueftsSoil::N_recovery, "N_recovery")
		.field("P_recovery", &QueftsSoil::P_recovery, "P_recovery")
		.field("K_recovery", &QueftsSoil::K_recovery, "K_recovery")
		.field("UptakeAdjust", &QueftsSoil::UptakeAdjust, "UptakeAdjust")
	;

    class_<QueftsCrop>("QueftsCrop")
		.field("NminStore", &QueftsCrop::NminStore, "NminStore")
		.field("NminVeg", &QueftsCrop::NminVeg, "NminVeg")
		.field("NmaxStore", &QueftsCrop::NmaxStore, "NmaxStore")
		.field("NmaxVeg", &QueftsCrop::NmaxVeg, "NmaxVeg")
		.field("PminStore", &QueftsCrop::PminStore, "PminStore")
		.field("PminVeg", &QueftsCrop::PminVeg, "PminVeg")
		.field("PmaxStore", &QueftsCrop::PmaxStore, "PmaxStore")
		.field("PmaxVeg", &QueftsCrop::PmaxVeg, "PmaxVeg")
		.field("KminStore", &QueftsCrop::KminStore, "KminStore")
		.field("KminVeg", &QueftsCrop::KminVeg, "KminVeg")
		.field("KmaxStore", &QueftsCrop::KmaxStore, "KmaxStore")
		.field("KmaxVeg", &QueftsCrop::KmaxVeg, "KmaxVeg")
		.field("Yzero", &QueftsCrop::Yzero, "Yzero")
		.field("Nfix", &QueftsCrop::Nfix, "Nfix") 
	;
	
    class_<QueftsModel>("QueftsModel")
		.constructor()
//		.method("run", &QueftsModel::run, "run") 
//		.method("output", &QueftsModel::output, "output") 
		.method("run", &runout, "run the model")

		.method("runbatch", &QueftsModel::runbatch, "run the model")
		
		.field("crop", &QueftsModel::crop, "crop")
		.field("soil", &QueftsModel::soil, "soil")

		// yield without nutrient limitation
		.field("leaf_att", &QueftsModel::leaf_att, "leaf_att")
		.field("stem_att", &QueftsModel::stem_att, "stem_att")
		.field("store_att", &QueftsModel::store_att, "store_att")
		// season length
		.field("SeasonLength", &QueftsModel::SeasonLength, "SeasonLength")

		// NPK input
		.field("N", &QueftsModel::N_fertilizer, "N_fertilizer")
		.field("P", &QueftsModel::P_fertilizer, "P_fertilizer")
		.field("K", &QueftsModel::K_fertilizer, "K_fertilizer")
		
		// output
/*		.field("N_supply", &QueftsModel::N_supply, "N_supply")
		.field("P_supply", &QueftsModel::P_supply, "P_supply")
		.field("K_supply", &QueftsModel::K_supply, "K_supply")

		// nutrient limited yield
		.field("leaf_lim", &QueftsModel::leaf_lim, "leaf_lim")
		.field("stem_lim", &QueftsModel::stem_lim, "stem_lim")
		.field("store_lim", &QueftsModel::store_lim, "store_lim")

		.field("N_gap", &QueftsModel::N_gap, "N_gap")
		.field("P_gap", &QueftsModel::P_gap, "P_gap")
		.field("K_gap", &QueftsModel::K_gap, "K_gap")
*/		
	;			
}
