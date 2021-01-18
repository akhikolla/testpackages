#include "get_util.h"
#include "chc.h"
#include "de.h"
#include "ssga.h"
#include "pso.h"
#include "sade.h"
#include "jade.h"
#include "sadeaf.h"
#include "jdelsgo.h"
#include "jdebin.h"
#include "jderand.h"
#include "jdemc.h"
#include "solis.h"
#include "cmaes.h"
#include "simplex.h"
#include "solisn.h"
#include "solisn2.h"
#include "cmaeshan.h"
#include "mts1.h"
#include "mts2.h"
#include "debug.h"
#include <iostream>
#include <cassert>
#include <cstdio>
#include <sstream>

double string_to_double( const std::string& s )
{
   std::istringstream i(s);
   double x;
   if (!(i >> x))
     return 0;
   return x;
} 
using namespace realea;
using std::unique_ptr;

void set_InitVerbose(void) {
    enable_print_info();
}

bool getCrossFactor(string name, double *pfactor) {
    int num;
    
    num = sscanf(name.c_str(), "ssgac-%lf", pfactor);
    return (num >= 1);
}

bool find_str(const string &str, const string substr) {
    return (str.find(substr)!=string::npos);
}

ICrossBinary *get_Cross(string crossover) {
    ICrossBinary *cross=NULL;
    char c_alpha[20];
    double alpha,pr;
    int num;

 
    if (find_str(crossover, "pblx")) { 

       if (find_str(crossover, "-") ) {
	    sscanf(crossover.c_str(), "pblx-%s", c_alpha);
	    alpha = string_to_double(c_alpha);
       }
       else {
	    alpha=0.5;
       }

       alpha = fabs(alpha);

       cross = new CrossPBLX(alpha);
    }
    else if (find_str(crossover, "blx")) {
       num = sscanf(crossover.c_str(), "blx-%s", c_alpha);
       alpha = string_to_double(c_alpha);

	if (num == 1) {
	   cross = new CrossBLX(alpha);
	}
    }
    else if (find_str(crossover, "dim")) {
       num = sscanf(crossover.c_str(), "dim-%s", c_alpha);
       alpha = string_to_double(c_alpha);
 //     pr = string_to_double(c_pr);
       pr = 0.8;
       print_info("alpha: %f\npr: %f\n", alpha, pr);
       cross = new CrossDim(alpha,pr);
    }

    return cross;
}

string get_EANames(string sep) {
    string result;
    result = "chc" + sep + "de" + sep +"jde" +sep +"jade" +sep +"sade" +sep +"sadeF" +sep +"sadeaf" +sep +"jdebin" +sep +"jdeexp" +sep;
    result = result + "jdemc" +sep +"jdemcinfo" +sep +"jderand" +sep;
    result = result + "ssga" + sep +"pso";
    return result;
}

IEA *get_EA(string alg, Random &random) {
	string print_str;
	ICrossBinary *cross=NULL;
	string cross_str;
	IEA *ea;

	print_str = "EA:: ";

	string::size_type pos_cross = alg.find("_");

	if (pos_cross != string::npos) {
	   cross_str = alg.substr(pos_cross+1);
	   alg = alg.substr(0, pos_cross);
	}
	else {
	    cross_str = "blx-0.5";
	}

	cross = get_Cross(cross_str);

	if (cross == NULL && (alg != "de" && alg != "pso")) {
	   throw string("Crossover " +alg +":" + cross_str +" is unknown");
	}
	if (alg == "ssga") {
	   SSGA *ssga = new SSGA(&random);
	   ssga->setCross(cross);
	   ssga->setMutation(new MutationBGA());
	   ssga->setSelect(new SelectNAM(3));
	   ssga->setReplacement(new ReplaceWorst());
           print_str += "SSGA with " +cross_str +" (by default NAM-3 and RW)";
	   ea = ssga;	
	}
	else if (alg == "ssga2") {
	   SSGA *ssga = new SSGA(&random);
	   ssga->setCross(cross);
	   ssga->setMutation(new MutationBGA());
	   ssga->setSelect(new SelectTournament(3));
	   ssga->setReplacement(new ReplaceDC());
	   print_str += "SSGA with " +cross_str +" with Tournament-3 and ReplaceDC";
	   ea = ssga;	
	}
	else if (alg == "chc") {
	   CHC *chc = new CHC(&random);
	   chc->setCross(cross);
	   print_str += "CHC with " + cross_str;
	   ea = chc;
	}
	else if (alg == "jade") {
	   JADE *jade= new JADE(&random);
	   print_str += "JADE with F=0.5";
	   ea = jade;
	}

	else if (alg == "pso") {
	  delete cross;
	  cross = NULL;
	  PSO *pso = new PSO(&random);
	  ea = pso;
	}
	else if (alg == "de") {
	   double f, cr;
	   DE *de = new DE(&random);

	   if (!cross_str.empty()) {
	    sscanf(cross_str.c_str(), "%lf,%lf", &f, &cr);
	   }
	   else {
	      f = 0.5;
	      cr = 0.9;
	   }
	   de->setF(f);
	   de->setCR(cr);
	   char val[50];
           sprintf(val, "DE with F=%.1f and CR=%.1f", f, cr);
	   print_str += val;
	   delete cross;
	   ea = de;
	}
        else if (alg == "jde") {
           JDE *de = new JDE(&random, 0);
	   char val[50];
	   sprintf(val, "JDE : not dynamic");
	   print_str += val;
	   ea = de;
	}
        else if (alg == "jdebin") {
           JDEBin *de = new JDEBin(&random, 0);
	   char val[50];
	   sprintf(val, "JDEBin : not dynamic");
	   print_str += val;
	   ea = de;
	}
        else if (alg == "jdeexp") {
           JDEBin *de = new JDEBin(&random, 0);
	   char val[50];
	   de->setStrategy("jDEexp");
	   sprintf(val, "JDEExp : not dynamic");
	   print_str += val;
	   ea = de;
	}

	else if (alg == "jdemc") {
           JDEMC *de = new JDEMC(&random, 0);
	   char val[50];
	   sprintf(val, "JDEMC : not dynamic");
	   print_str += val;
	   ea = de;
	}
	else if (alg == "jderand") {
           JDERand *de = new JDERand(&random, 0);
	   char val[50];
	   sprintf(val, "JDERand : not dynamic");
	   print_str += val;
	   ea = de;
	}
	else if (alg == "jdemcinfo") {
           JDEMC *de = new JDEMC(&random, 0);
	   de->setDebug();
	   char val[50];
	   sprintf(val, "JDEMC : not dynamic");
	   print_str += val;
	   ea = de;
	}
        else if (alg == "sade") {
           SADE *de = new SADE(&random);
           char val[50];
	   sprintf(val, "SaDE: averageF: 0.5");
	   print_str += val;
           ea = de;
        }
        else if (alg == "sadeaf") {
           SADEAF *sadeaf = new SADEAF(&random);
           sadeaf->setAverageF(0.5);
           char val[50];
	   sprintf(val, "SaDEAF: averageF: 0.5");
	   print_str += val;
           ea = sadeaf;
        }
	else if (alg.find("sadeF")==0) {
           SADE *de = new SADE(&random);
           double averageF = atof(alg.substr(5).c_str());
 	   de->setAverageF(averageF);
           char val[50];
           sprintf(val, "SaDE\taverageF: %f", averageF);
	   print_str += val;
           ea = de;
        }
	else {
	    delete cross;
	    throw string("EA '" +alg +"' is unknown");
	}
	
	print_info("%s\n", print_str.c_str());

	return ea;
}

string get_LSNames(string sep) {
    return "cmaes" + sep + "sw" + sep + "swn" + sep +"diswf" +sep +"diswa" + sep + "simplex" + sep +"mtsls1" +sep + "mts1" + sep +"mts2";
}


ILocalSearch *get_LS(string arg_ls, DomainRealPtr domain, Random *random) {
    ILocalSearch *ls;
    string print_str = "LS: ";
	
	if (arg_ls == "sw") {
	   SolisWets *sw = new SolisWets();
	   sw->setDelta(0.2);
	   print_str += "Solis Wets\nSW::sigma : Sigma 0.2";
	   ls = sw;
	}
	else if (arg_ls.find("swn")==0) {
	   SWNDim *sw = new SWNDim();
	   string c_strategy = arg_ls.substr(3);
	   int strategy = atoi(c_strategy.c_str());
	   sw->setStrategy(strategy);
	   sw->setDelta(0.2);
	   print_str += "Solis Wets NDim\tSW::sigma : Sigma 0.2\tStrategy: " +c_strategy;
	   ls = sw;
	}
	else if (arg_ls.find("ssw")==0) {
	   SWN2Dim *sw = new SWN2Dim();
	   int strategy = 3;
	   sw->setStrategy(strategy);
	   sw->setDelta(1e-15,0.4);
	   print_str += "Solis Wets NDim\tSW::sigma : Sigma 0.2";
	   ls = sw;
	}
       else if (arg_ls == "mts1") {
	    ls = new MTSLS1(0.4, 1e-15);
	    print_str += "MTS1 Local Search\nMTS::maxsigma : 0.4\nMTS::minsigma : 1E-15";
	}
        else if (arg_ls == "mtsls1") {
            ls = new MTSLS1(0.4, 1e-15);
            print_str += "MTS1 Local Search\nMTS::maxsigma : 0.4\nMTS::minsigma : 1E-15";
        }
	else if (arg_ls == "mts2") {
	    ls = new MTSLS2(0.4, 1e-15);
	    print_str += "MTS2 Local Search\nMTS::maxsigma : 0.4\nMTS::minsigma : 1E-15";
	}
	else if (arg_ls == "simplex") {
	    ls = new SimplexDim();
	    print_str += "SimplexDim";
	}
	else if (arg_ls.find("cmaes")!=string::npos) {
 	   CMAESHansen *cmaes = new CMAESHansen("cmaesinit.par");
 	   cmaes->searchNeighborhood(0.5);
	   string strategy = "";

	   if (arg_ls == "cmaesno") 
	      cmaes->disableBoundChecking();
           else if (arg_ls == "cmaesalways" || arg_ls == "cmaes")
	      cmaes->enableBoundChecking();
           else if (arg_ls == "cmaesmyrandom") {
               cmaes->enableBoundChecking();
               cmaes->setMyRandom();
           }
	   else {
		throw string("localsearch '" +arg_ls +"' is unknown");
	   }

	   ls = cmaes;
           print_str += "CMAESHansen: " +arg_ls +"\nCMAES::Neighborhood: 0.5";
	}
	else {
		throw string("localsearch '" +arg_ls +"' is unknown");
	}

     print_debug("%s\n", print_str.c_str());

     return ls;
}

void set_Effort(Hybrid *hybrid, string effort) {
    double ratio;

    if (effort != "") {
       ratio = string_to_double(effort);
       assert(ratio > 0 && ratio < 1);
    }
    else {
       ratio = 0.5;
    }

     print_debug("LS::Effort: %f\n", ratio);
     hybrid->setEffortRatio(ratio);
}

void set_MaxEval(IEA *ea, string maxeval) {
   if (maxeval != "") {
      int imaxeval = atoi(maxeval.c_str());
      assert(imaxeval > 0);
      unsigned maxeval = imaxeval; 

      print_debug("EA::MaxEval: %u\n", maxeval);
      ea->setMaxEval(maxeval);

   }

}
