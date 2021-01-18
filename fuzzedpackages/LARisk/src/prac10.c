#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


double doseft(char* dist, double* info);						// function to calculate dose
double latencyft(char* site);							      		// function to calculate latency
double weightft(char* site);    									  // function to return weight
void data_call_1();							          					// function to retrieve dbaseline data
void data_call_2();											          	// function to retrieve dincidence data
void data_call_3(double * data_baseline);						// function to retrieve NEW dbaseline data
void data_call_4(double * data_incidence);				                	// function to retrieve NEW dincidence data
double survft(int age, char* sex);					          	  	      	// function to return surv
double srateft(int e_age, int a_age, char* sex);		                // function to return srate
double incidenceft(int age, char* sex, char* site);                	// function to return incidence
void brft(int* sex1, int* site1, int *age, int* changedata, double* dbaseline, double* dincidence , double* out);   				// function to calculate brft baserisk
double vmvft(double* vec, double(*mat)[5], int column);			        // function to calculate vec * mat * vex
double DDREFft(char * DDREF_op, char * exposure_rate, double dose); // function to calculate DDREF
double approxft(double a, double b, double c);						          // function to calculate approx
double r_trift(double a, double b, double c);           						// function to generate random number from triangular distribution

// structure declaration for data entry
struct baseline
{
  double Prob_d_m;
  double Prob_d_f;
};
struct incidence
{
  double Rate_m;
  double Rate_f;
};

// Global Variables Declaration
struct baseline dbaseline[101];						// Global variable declaration for storing dbaseline data
struct incidence dincidence[1919];				// Global variable declaration for storing dincidence data

char * weightmat_site[19] = { "stomach","colon","liver","lung","breast","ovary","uterus","prostate","bladder","brain/cns","thyroid","remainder","oral","oesophagus","rectum","gallbladder","pancreas","kidney","leukemia" };
double weightmat_value[19] = { 0.7,0.7,0.7,0.3,0,0.7,0.7,0.7,0.7,1,1,0.7,0.7,0.7,0.7,1,0.7,0.7,0.7 };


void larft(int * dosedist1, double * doseinfo, int * sex1, int * site1,
           int * exposure_rate1, int * age, double * exposeage, int * DDREF_op1, int * sim ,double * ci,
           int * weight_site, double * weight_value, int * weight_length,
           int * changedata, double * dbaseline, double * dincidence ,double * out)
{

  if(*changedata==2)
  {
    data_call_1();			// function to retrieve dbaseline data > Run before survft and srateft
    data_call_2();
  }
  else
  {
    data_call_3(dbaseline);
    data_call_4(dincidence);
  }



  // Variable Declaration
  char * sex, * exposure_rate, * site , * DDREF_op, * dosedist;
  int sum_el, len;
  double latency, estar=0, DDREF, weight, z_value, dose;

  //Re-declaration after taking over
  if (*sex1 == 1)
    sex = "male";
  else
    sex = "female";

  switch (*site1)
  {
  case 1:
    site = "stomach";
    break;
  case 2:
    site = "colon";
    break;
  case 3:
    site = "liver";
    break;
  case 4:
    site = "lung";
    break;
  case 5:
    site = "breast";
    break;
  case 6:
    site = "ovary";
    break;
  case 7:
    site = "uterus";
    break;
  case 8:
    site = "prostate";
    break;
  case 9:
    site = "bladder";
    break;
  case 10:
    site = "brain/cns";
    break;
  case 11:
    site = "thyroid";
    break;
  case 12:
    site = "remainder";
    break;
  case 13:
    site = "oral";
    break;
  case 14:
    site = "oesophagus";
    break;
  case 15:
    site = "rectum";
    break;
  case 16:
    site = "gallbladder";
    break;
  case 17:
    site = "pancreas";
    break;
  case 18:
    site = "kidney";
    break;
  case 19:
    site = "leukemia";
    break;
  default:
    site = "";

  }

  if (*exposure_rate1 == 1)
    exposure_rate = "acute";
  else
    exposure_rate = "chronic";


  switch (*dosedist1)
  {
  case 1:
    dosedist = "fixedvalue";
    break;
  case 2:
    dosedist = "lognormal";
    break;
  case 3:
    dosedist = "normal";
    break;
  case 4:
    dosedist = "triangular";
    break;
  case 5:
    dosedist = "logtriangular";
    break;
  case 6:
    dosedist = "uniform";
    break;
  case 7:
    dosedist = "loguniform";
    break;
  default:
    dosedist = "fixedvalue";
  }

  if (*DDREF_op1 == 1)
    DDREF_op = "T";
  else
    DDREF_op = "F";


  // change weight
  if (1 <= *weight_length )
  {
    for (int i = 0; i < *weight_length; i++)
    {

      for (int j= 0; j < 19; j++)
      {
        if (weight_site[i]==j+1)
        {
          weightmat_value[j] = weight_value[i];
          break;
        }
      }
    }

  }



  // matrix for cleanup result
  /*
  double LARMAT_leukemia[10][3];		//Assignment by the number of data
  double F_LARMAT_leukemia[10][3];	//Assignment by the number of data
  double SUMMAT[10][300];				//Dynamic allocation by number of data and secondary array
  double F_SUMMAT[10][300];			//Dynamic allocation by number of data and secondary array
  */

  //z_value < -qnorm();				from normal dist.
  GetRNGstate();
  z_value = qnorm(0.5+(*ci/2.0), 0.0, 1.0, 1, 0);
  PutRNGstate();




  if(!strcmp(site, "leukemia"))
  {
    // decrlaration
    /*
    double * LAR_total = (double *)malloc(sizeof(double)* sim);				// Dynamic allocation : use stdlib.h
    double * var_log_ERR = (double *)malloc(sizeof(double)* sim);			// Dynamic allocation
    double * var_log_EAR = (double *)malloc(sizeof(double)* sim);			// Dynamic allocation
    double * LAR_total_SD = (double *)malloc(sizeof(double)* sim);			// Dynamic allocation
    double * F_LAR_total = (double *)malloc(sizeof(double)* sim);			// Dynamic allocation
    double * F_var_log_ERR = (double *)malloc(sizeof(double)* sim);			// Dynamic allocation
    double * F_var_log_EAR = (double *)malloc(sizeof(double)* sim);			// Dynamic allocation
    double * F_LAR_total_SD = (double *)malloc(sizeof(double)* sim);		// Dynamic allocation
    */


    double LAR_total[*sim];				// Dynamic allocation : use stdlib.h
    double var_log_ERR[*sim];			// Dynamic allocation
    double var_log_EAR[*sim];			// Dynamic allocation
    double LAR_total_SD[*sim];			// Dynamic allocation
    double F_LAR_total[*sim];			// Dynamic allocation
    double F_var_log_ERR[*sim];			// Dynamic allocation
    double F_var_log_EAR[*sim];			// Dynamic allocation
    double F_LAR_total_SD[*sim];		// Dynamic allocation


    double leukemia_tmp[3] = { 0,0,0 }, F_leukemia_tmp[3] = { 0,0,0 };

    // decrlaration sigma
    double sigma_ERR_m[5][5] = {	{0.349923, -0.371562, 0.0279017, 0.0277847, 0.0206183},	{-0.371562, 0.542042, 0.0101649, 0.00150142, -0.00642578},	{0.0279017, 0.0101649, 0.0373066, 0.0169665, 0.0134684}, {0.0277847, 0.00150142, 0.0169665, 0.110026, 0.0578845}, {0.0206183, -0.00642578, 0.0134684, 0.0578845, 0.0696666}	};
    double sigma_EAR_m[5][5] = {	{0.677437, -0.545608, 0.0296734, 0, -0.00472918},	{-0.545608, 0.542886, 0.00169235, 0, 0.00737659},	{0.0296734, 0.00169235, 0.0217357, 0, 0.0103813},	{0, 0, 0, 0, 0},	{-0.00472918, 0.00737659, 0.0103813, 0, 0.0158421}		};
    double sigma_ERR_f[5][5] = {	{0.420773, -0.399592, 0.0270153, 0.0244283, 0.0188047},	{-0.399592, 0.542042, 0.0101649, 0.00150142, -0.00642578},	{0.0270153, 0.0101649, 0.0373066, 0.0169665, 0.0134684},	{0.0244283, 0.00150142, 0.0169665, 0.110026, 0.0578845},	{0.0188047, -0.00642578, 0.0134684, 0.0578845, 0.0696666}	};
    double sigma_EAR_f[5][5] = {	{0.216888, -0.306358, 0.0150695, 0, -0.00201906},	{-0.306358, 0.542886, 0.00169235, 0, 0.00737659},	{0.0150695, 0.00169235, 0.0217357, 0, 0.0103813},	{0, 0, 0, 0, 0},		{-0.00201906, 0.00737659, 0.0103813, 0, 0.0158421}		};

    // leukemia beta
    double gamma_ERR = -0.4, delta_ERR = -0.48, phi_ERR = 0.42, theta_ERR = 0.87;
    double gamma_EAR = 0.29, delta_EAR = 0, phi_EAR = 0.56, theta_EAR = 0.88;
    double beta_ERR, beta_EAR;

    if (!strcmp(sex, "male"))
    {
      beta_ERR = 1.1;
      beta_EAR = 1.62;
    }
    else
    {
      beta_ERR = 1.2;
      beta_EAR = 0.93;
    }

    if (!strcmp(exposure_rate,"chronic")) {
      theta_ERR = 0;
      theta_EAR = 0;
    }

    for (int j = 0; j < *sim; j++)
    {
      //dose
      dose = doseft(dosedist, doseinfo);

      //latency
      latency=latencyft(site);

      //exposeage + latency
      sum_el=(int)round(*exposeage + latency);

      len = (100 - sum_el) + 1;				// to calculate length, +1



      double * LAR_ERR = (double *)malloc(sizeof(double)* len);			// dynamic allocation : stdlib.h??????
      double * LAR_EAR = (double *)malloc(sizeof(double)* len);			// dynamic allocation
      double * denominator_ERR = (double *)malloc(sizeof(double)* len);	// dynamic allocation
      double * numerator_ERR = (double *)malloc(sizeof(double)* len);		// dynamic allocation
      double * denominator_EAR = (double *)malloc(sizeof(double)* len);	// dynamic allocation
      double * numerator_EAR = (double *)malloc(sizeof(double)* len);		// dynamic allocation
      double S_LAR_ERR = 0, S_LAR_EAR = 0, S_denominator_ERR = 0, S_numerator_ERR = 0, S_denominator_EAR = 0, S_numerator_EAR = 0;	// for sum

      double incrate, srate,log_cal;
      for (int i = sum_el; i < 101; i++)										// Repeat by the entered size
      {
        if (*exposeage < (double)30)
          estar = ((*exposeage - (double)30) / (double)10);
        else
          estar = (double)0;

        incrate = incidenceft(i, sex, site);
        srate = srateft((int)*exposeage, i, sex);
        if(srate <= (double)0)
          srate = pow(0.1, 100);
        log_cal = log((i - (int)*exposeage) / (double)25);

        LAR_ERR[i - sum_el] = incrate * beta_ERR*dose*((double)1 + theta_ERR * dose)*exp((gamma_ERR*estar) + (delta_ERR*log_cal) + (phi_ERR*estar*log_cal))*srate;
        LAR_EAR[i - sum_el] = 10 * beta_EAR*dose*((double)1 + theta_EAR * dose)*exp((gamma_EAR*estar) + (delta_EAR*log_cal) + (phi_EAR*estar*log_cal))*srate;

        numerator_ERR[i - sum_el] = exp((delta_ERR + phi_ERR * estar)*log_cal)*incrate*srate*log_cal;
        denominator_ERR[i - sum_el] = exp((delta_ERR + phi_ERR * estar)*log_cal)*incrate*srate;

        numerator_EAR[i - sum_el] = exp((delta_EAR + phi_EAR * estar)*log_cal)*srate*log_cal;
        denominator_EAR[i - sum_el] = exp((delta_EAR + phi_EAR * estar)*log_cal)*srate;

        S_LAR_ERR += LAR_ERR[i - sum_el];
        S_LAR_EAR += LAR_EAR[i - sum_el];
        S_numerator_ERR += numerator_ERR[i - sum_el];
        S_denominator_ERR += denominator_ERR[i - sum_el];
        S_numerator_EAR += numerator_EAR[i - sum_el];
        S_denominator_EAR += denominator_EAR[i - sum_el];

      }

      if(sum_el >= 100)
      {
        S_LAR_ERR = pow(0.1,100);
        S_LAR_EAR = pow(0.1,100);

        S_numerator_EAR = pow(0.1,100);
        S_denominator_EAR = pow(0.1,100);
      }

      DDREF = 1;

      weight = weightft(site);

      LAR_total[j] = pow(S_LAR_ERR, weight) * pow(S_LAR_EAR, ((double)1 - weight)) / DDREF;

      double F_S_LAR_ERR = 0, F_S_LAR_EAR = 0;					// for F_SUM
      double F_S_numerator_ERR = 0, F_S_denominator_ERR = 0;
      double F_S_numerator_EAR = 0, F_S_denominator_EAR = 0;

      if (*age <= sum_el)
        F_LAR_total[j] = LAR_total[j];
      else
      {
        for (int i = (*age - (int)sum_el); i < len ; i++)		// dose not exist function for sum
        {
          F_S_LAR_ERR += LAR_ERR[i];
          F_S_LAR_EAR += LAR_EAR[i];

          F_S_numerator_ERR += numerator_ERR[i];
          F_S_denominator_ERR += denominator_ERR[i];
          F_S_numerator_EAR += numerator_ERR[i];
          F_S_denominator_EAR += denominator_ERR[i];
        }

        if(F_S_LAR_ERR <= (double)0)
        {
          F_S_LAR_ERR = pow(0.1,100);
          F_S_LAR_EAR = pow(0.1,100);

          F_S_numerator_EAR = pow(0.1,100);
          F_S_denominator_EAR = pow(0.1,100);
        }

        F_LAR_total[j] = pow(F_S_LAR_ERR, weight) * pow(F_S_LAR_EAR, ((double)1 - weight)) / DDREF;

      }

      // var (log(LAR_ERR)) calculation
      double sum_n_d_ERR = 0, F_sum_n_d_ERR = 0;

      sum_n_d_ERR = S_numerator_ERR / S_denominator_ERR;
      F_sum_n_d_ERR = F_S_numerator_ERR / F_S_denominator_ERR;

      double atrans_ERR[5] = { (double)1 / beta_ERR ,dose / ((double)1 + theta_ERR * dose), estar, sum_n_d_ERR, sum_n_d_ERR*estar };		// vetor

      // var (log(LAR_EAR)) calculation
      double sum_n_d_EAR = 0, F_sum_n_d_EAR = 0;

      sum_n_d_EAR = S_numerator_EAR / S_denominator_EAR;
      F_sum_n_d_EAR = F_S_numerator_EAR / F_S_denominator_EAR;

      double  atrans_EAR[5] = { (double)1 / beta_EAR ,dose / ((double)1 + theta_EAR * dose), estar, sum_n_d_EAR, sum_n_d_EAR*estar };		// vetor


      if (!strcmp(sex, "male"))
      {
        var_log_ERR[j] = vmvft(atrans_ERR, sigma_ERR_m,sizeof(sigma_ERR_m)/sizeof(sigma_ERR_m[0]));
        var_log_EAR[j] = vmvft(atrans_EAR, sigma_EAR_m, sizeof(sigma_EAR_m) / sizeof(sigma_EAR_m[0]));
      }
      else
      {
        var_log_ERR[j] = vmvft(atrans_ERR, sigma_ERR_f, sizeof(sigma_ERR_f) / sizeof(sigma_ERR_f[0]));
        var_log_EAR[j] = vmvft(atrans_EAR, sigma_EAR_f, sizeof(sigma_EAR_f) / sizeof(sigma_EAR_f[0]));
      }

      LAR_total_SD[j] = sqrt( var_log_EAR[j] + ( ( pow(( log(S_LAR_ERR / S_LAR_EAR)),2)) * weight*((double)1 - weight) ) );

      // Start here

      if (*age <= sum_el)
      {
        F_LAR_total_SD[j] = LAR_total_SD[j];
      }

      else
      {
        double F_atrans_ERR[5] = { (double)1 / beta_ERR, dose / ((double)1 + theta_ERR * dose), estar, (F_sum_n_d_ERR), (F_sum_n_d_ERR)*estar};
        double F_atrans_EAR[5] = { (double)1 / beta_EAR, dose / ((double)1 + theta_EAR * dose), estar, (F_sum_n_d_EAR), (F_sum_n_d_EAR)*estar};
        if (!strcmp(sex, "male"))
        {
          F_var_log_ERR[j] = vmvft(F_atrans_ERR, sigma_ERR_m, sizeof(sigma_ERR_m) / sizeof(sigma_ERR_m[0]));
          F_var_log_EAR[j] = vmvft(F_atrans_EAR, sigma_EAR_m, sizeof(sigma_EAR_m) / sizeof(sigma_EAR_m[0]));
        }
        else
        {
          F_var_log_ERR[j] = vmvft(F_atrans_ERR, sigma_ERR_f, sizeof(sigma_ERR_f) / sizeof(sigma_ERR_f[0]));
          F_var_log_EAR[j] = vmvft(F_atrans_EAR, sigma_EAR_f, sizeof(sigma_EAR_f) / sizeof(sigma_EAR_f[0]));
        }
        F_LAR_total_SD[j] = sqrt( F_var_log_EAR[j] + ( ( pow(log(F_S_LAR_ERR / F_S_LAR_EAR),2) )*weight*((double)1 - weight)   )   );

      }

      if(var_log_ERR[j]<0) var_log_ERR[j] = var_log_ERR[j];
      if(F_var_log_ERR[j]<0) F_var_log_ERR[j] = F_var_log_ERR[j];



      free(LAR_ERR);
      free(LAR_EAR);
      free(denominator_ERR);
      free(numerator_ERR);
      free(denominator_EAR);
      free(numerator_EAR);

      //  To put a value in LARMAT_leukemia

      leukemia_tmp[0] += exp(log(LAR_total[j]) - z_value * LAR_total_SD[j]);
      leukemia_tmp[1] += LAR_total[j];
      leukemia_tmp[2] += exp(log(LAR_total[j]) + z_value * LAR_total_SD[j]);

      F_leukemia_tmp[0] += exp(log(F_LAR_total[j]) - z_value * F_LAR_total_SD[j]);
      F_leukemia_tmp[1] += F_LAR_total[j];
      F_leukemia_tmp[2] += exp(log(F_LAR_total[j]) + z_value * F_LAR_total_SD[j]);


      //SUMMAT[0][j] = LAR_total[j];		// k
      //F_SUMMAT[0][j] = F_LAR_total[j];	// k
    }



    for(int i=0;i<(2*(*sim)+6);i++)
    {
      if(i<*sim)
        out[i]=LAR_total[i];

      else if((*sim <= i) & (i<2*(*sim)))
        out[i]=F_LAR_total[i-(*sim)];

      else if((2*(*sim) <= i) & (i < (2*(*sim)+3)))
        out[i]=leukemia_tmp[i-2*(*sim)] / *sim;

      else if(((2*(*sim)+3) <= i) & (i < (2*(*sim)+6)))
        out[i]=F_leukemia_tmp[i-(2*(*sim)+3)] / *sim;
      //out[i]=1;

    }



    //LARMAT_leukemia[0][i] = leukemia_tmp[i] / *sim;		// k
    //F_LARMAT_leukemia[0][i] = leukemia_tmp[i] / *sim;		// k

    /*
    free(LAR_total);
    free(var_log_ERR);
    free(var_log_EAR);
    free(LAR_total_SD);
    free(F_LAR_total);
    free(F_var_log_ERR);
    free(F_var_log_EAR);
    free(F_LAR_total_SD);*/

  }
  else    // end leukemia
  {

    // declaration

    //double * LAR_total = (double *)malloc(sizeof(double)* sim);					//dynamic allocation : stdlib.h??????
    //double * F_LAR_total = (double *)malloc(sizeof(double)* sim);				//dynamic allocation

    double LAR_total[*sim];
    double F_LAR_total[*sim];

    // declaration beta
    double gamma_ERR = (-0.3), eta_ERR = (-1.4), gamma_EAR = (-0.41), eta_EAR = 2.8;
    if (!strcmp(site, "breast"))
      gamma_EAR = (-0.5);
    else if (!strcmp(site, "oral"))
      eta_EAR = 0.5;
    else if (!strcmp(site, "liver"))
      eta_EAR = 4.1;
    else if (!strcmp(site, "lung"))
      eta_EAR = 5.2;
    else if (!strcmp(site, "bladder"))
      eta_EAR = 6.0;
    else if (!strcmp(site, "thyroid"))
    {
      gamma_ERR = (-0.83);
      eta_ERR = 0.0;
    }

    //betadist declaration :
    // row : "stomach","colon","liver","lung","breast","ovary","uterus","prostate","bladder","brain/cns","thyroid","remainder"
    // column : Mean_ERR_m, Sd_ERR_m, Mean_ERR_f, Sd_ERR_f, Mean_EAR_m, Sd_EAR_m, Mean_EAR_f, Sd_EAR_f

    //double betadist[12][8] = {
    //	{-1.5606477,0.3293327,-0.73396918,0.21848782,1.5892352,0.3042856,1.58923521,0.2103887},
    //	{-0.4620355,0.2779496,-0.84397007,0.41324215,1.1631508,0.2895357,0.47000363,0.3536465},
    //	{-1.1394343,0.3536465,-1.13943428,0.58739416,0.7884574,0.4523131,0,0.4674953},
    //	{-1.1394343,0.3929707,0.33647224,0.20505427,0.8329091,0.3862571,1.22377543,0.1929403},
    //	{0,0,0,0,0,0,2.30258509,0.1804418},
    //	{0,0,-0.96758403,0.67322891,0,0,-0.35667494,0.5998406},
    //	{0,0,0.055,0.08418367,0,0,1.2,0.7142857},
    //	{0.12,0.2908163,0,0,0.11,0.4540816,0,0},
    //	{-0.6931472,0.5232833,0.50077529,0.44830562,0.1823216,0.567506,-0.28768207,0.4425003},
    //	{-0.3424903,0.4183019,-1.42711636,0.42166404,0,0,0,0},
    //	{-0.6348783,0.6783827,0.04879016,0.67192404,0,0,0,0},
    //	{-0.1392621,0.3375603,-0.22314355,0.45055679,1.0043016,0.2894181,0.05826891,0.3922258}
    //};

    double beta_ERR = 0.0, beta_EAR = 1.0;
    double astar = 0.0;

    for (int j = 0; j < *sim; j++)
    {
      dose = doseft(dosedist, doseinfo);

      //	"stomach","colon","liver","lung","breast","ovary","uterus","prostate","bladder","brain/cns","thyroid","remainder"
      //	"oral", "oesophagus","rectum","gallbladder", "pancreas","kidney"

      if (!strcmp(site, "stomach")) {
        if (!strcmp(sex, "male")) {
          GetRNGstate();
          beta_ERR = exp(rnorm(-1.5606477,0.3293327));
          beta_EAR = exp(rnorm(1.5892352,0.3042856));
          PutRNGstate();
        }
        else {
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.73396918,0.21848782));
          beta_EAR = exp(rnorm(1.58923521,0.2103887));
          PutRNGstate();
        }
      }
      else if (!strcmp(site, "colon")) {
        if (!strcmp(sex, "male")) {
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.4620355,0.2779496));
          beta_EAR = exp(rnorm(1.1631508,0.2895357));
          PutRNGstate();
        }
        else {
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.84397007,0.41324215));
          beta_EAR = exp(rnorm(0.47000363,0.3536465));
          PutRNGstate();
        }
      }
      else if (!strcmp(site, "liver")) {
        if (!strcmp(sex, "male")) {
          GetRNGstate();
          beta_ERR = exp(rnorm(-1.1394343,0.3536465));
          beta_EAR = exp(rnorm(0.7884574,0.4523131));
          PutRNGstate();
        }
        else {
          GetRNGstate();
          beta_ERR = exp(rnorm(-1.13943428,0.58739416));
          beta_EAR = exp(rnorm(0,0.4674953));
          PutRNGstate();
        }
      }
      else if (!strcmp(site, "lung")) {
        if (!strcmp(sex, "male")) {
          GetRNGstate();
          beta_ERR = exp(rnorm(-1.1394343,0.3929707));
          beta_EAR = exp(rnorm(0.8329091,0.3862571));
          PutRNGstate();
        }
        else {
          GetRNGstate();
          beta_ERR = exp(rnorm(0.33647224,0.20505427));
          beta_EAR = exp(rnorm(1.22377543,0.1929403));
          PutRNGstate();
        }
      }
      else if (!strcmp(site, "breast")) {
        beta_ERR = (double)0;
        GetRNGstate();
        beta_EAR = exp(rnorm(2.30258509,0.1804418));
        PutRNGstate();
      }
      else if (!strcmp(site, "ovary")) {
        GetRNGstate();
        beta_ERR = exp(rnorm(-0.96758403,0.67322891));
        beta_EAR = exp(rnorm(-0.35667494,0.5998406));
        PutRNGstate();
      }
      else if (!strcmp(site, "uterus")) {
        GetRNGstate();
        beta_ERR = rnorm(0.055, 0.08418367);
        beta_EAR = rnorm(1.2, 0.7142857);
        PutRNGstate();
      }
      else if (!strcmp(site, "prostate")) {
        GetRNGstate();
        beta_ERR = rnorm(0.12, 0.2908163);
        beta_EAR = rnorm(0.11, 0.4540816);
        PutRNGstate();
      }
      else if (!strcmp(site, "bladder")) {
        if (!strcmp(sex, "male")) {
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.6931472,0.5232833));
          beta_EAR = exp(rnorm(0.1823216,0.567506));
          PutRNGstate();
        }
        else {
          GetRNGstate();
          beta_ERR = exp(rnorm(0.50077529,0.44830562	));
          beta_EAR = exp(rnorm(-0.28768207,0.4425003	));
          PutRNGstate();
        }
      }
      else if (!strcmp(site, "brain/cns")) {
        if (!strcmp(sex, "male")) {
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.3424903,0.4183019));
          beta_EAR = (double)0;
          PutRNGstate();
        }
        else {
          GetRNGstate();
          beta_ERR = exp(rnorm(-1.42711636,0.42166404));
          beta_EAR = (double)0;
          PutRNGstate();
        }
      }
      else if (!strcmp(site, "thyroid")) {
        if (!strcmp(sex, "male")) {
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.6348783,0.6783827));
          beta_EAR = (double)0;
          PutRNGstate();
        }
        else {
          GetRNGstate();
          beta_ERR = exp(rnorm(0.04879016,0.67192404));
          beta_EAR = (double)0;
          PutRNGstate();
        }
      }
      else if (!strcmp(site, "remainder")) {
        if (!strcmp(sex, "male")) {
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.1392621,0.3375603));
          beta_EAR = exp(rnorm(1.0043016,0.2894181));
          PutRNGstate();
        }
        else {
          GetRNGstate();
          beta_ERR = exp(rnorm(-0.22314355,0.45055679));
          beta_EAR = exp(rnorm(0.05826891,0.3922258));
          PutRNGstate();
        }
      }
      else if (!strcmp(site, "oral")) {
        if (!strcmp(sex, "male")) {
          beta_ERR = approxft(-0.0047, 0.23, 0.66);
          beta_EAR = approxft(0.08, 0.44, 1.09);
        }
        else {
          beta_ERR = approxft(0.13, 0.53, 1.24);
          beta_EAR = approxft(0.06, 0.29, 0.66);
        }

      }
      else if (!strcmp(site, "oesophagus")) {
        if (!strcmp(sex, "male")) {
          beta_ERR = approxft(0.093, 0.51, 1.13);
          beta_EAR = approxft(0.11, 0.88, 2.10);
        }
        else {
          beta_ERR = approxft(-0.10, 0.82, 3.09);
          beta_EAR = approxft(-0.08, 0.14, 0.63);
        }
      }
      else if (!strcmp(site, "rectum")) {
        beta_ERR = approxft(-0.038, 0.12, 0.38);
        beta_EAR = approxft(-0.104, 0.34, 1.09);
      }
      else if (!strcmp(site, "gallbladder")) {
        beta_ERR = approxft(-0.39, -0.018, 0.29);
        beta_EAR = (double)0;
      }
      else if (!strcmp(site, "pancreas")) {
        beta_ERR = approxft(-0.0055, 0.36, 0.88);
        beta_EAR = approxft(0.0904, 0.49, 1.08);
      }
      else if (!strcmp(site, "kidney")) {
        beta_ERR = approxft(-0.012, 0.34, 1.00);
        beta_EAR = approxft(0.082, 0.31, 0.68);
      }

      //latency
      latency = latencyft(site);

      //exposeage + latency
      sum_el = (int)round(*exposeage + latency);

      len = (100 - sum_el) + 1;				// To calculate length, +1

      double * LAR_ERR = (double *)malloc(sizeof(double)* len);			// dynamic allocation : use stdlib.h
      double * LAR_EAR = (double *)malloc(sizeof(double)* len);			// dynamic allocation

      double srate = 0, S_LAR_ERR = 0, S_LAR_EAR = 0;
      for (int i = sum_el; i < 101; i++)										// Repeat by the entered size
      {
        // making breast eta_EAR
        if (!strcmp(site, "breast") && i <= 50) {
          eta_EAR = (double)3.5;
        }
        else {
          if (!strcmp(site,"breast"))
            eta_EAR = (double)1;
        }

        // making astar,estar
        if ((int)*exposeage < 30)
          estar = ((*exposeage - (double)30) / (double)10);
        else
          estar = (double)0;

        astar = i / (double)60;

        if (!strcmp(site, "breast")) {
          estar = ((*exposeage - (double)25) / (double)10);
          astar = i / (double)50;
        }
        if (!strcmp(site, "thyroid")) {
          estar = ((*exposeage - (double)30) / (double)10);
        }

        srate = srateft((int)*exposeage, i, sex);
        LAR_ERR[i - sum_el] = incidenceft(i, sex, site)    * (beta_ERR*dose) * exp(gamma_ERR*(estar)) *(pow((astar), eta_ERR))   *srate;
        LAR_EAR[i - sum_el] = (double)10 * beta_EAR * dose * exp(gamma_EAR*(estar))*pow((astar), eta_EAR) * srate;

        S_LAR_ERR += LAR_ERR[i - sum_el];
        S_LAR_EAR += LAR_EAR[i - sum_el];
      }

      //DDREF
      DDREF = DDREFft(DDREF_op, exposure_rate, dose);

      //weight
      weight = weightft(site);

      // LAR_total
      LAR_total[j] = (S_LAR_ERR * (weight) + S_LAR_EAR * ((double)1 - weight)) / DDREF;


      double F_S_LAR_ERR = 0, F_S_LAR_EAR = 0;

      if (*age <= sum_el)
        F_LAR_total[j] = LAR_total[j];
      else
      {
        for (int i = (*age - (int)sum_el); i < len; i++)
        {
          F_S_LAR_ERR += LAR_ERR[i];
          F_S_LAR_EAR += LAR_EAR[i];
        }
        F_LAR_total[j] = (F_S_LAR_ERR * (weight) + F_S_LAR_EAR * ((double)1 - weight)) / DDREF;
      }

      //SUMMAT[0][j] = LAR_total[j];		// k
      //F_SUMMAT[0][j] = F_LAR_total[j];	// k


      free(LAR_ERR);
      free(LAR_EAR);
    }

    for(int i=0;i<(2*(*sim));i++)
    {
      if(i<*sim)
        out[i]=LAR_total[i];
      else
        out[i]=F_LAR_total[i-(*sim)];
    }

  }


}



// function to calculate dose
double doseft(char* dist, double* info)
{
  char *dist1[7] = { "fixedvalue","lognormal","normal","triangular","logtriangular","uniform","loguniform" };
  double doseresult;

  if (!strcmp(dist, dist1[0]))
    doseresult = info[0] / 1000.0;

  else if (!strcmp(dist, dist1[1])) {
    GetRNGstate();
    doseresult = rlnorm(log(info[0]), log(info[1])) / 1000.0;
    PutRNGstate();
  }

  else if (!strcmp(dist, dist1[2])) {
    GetRNGstate();
    doseresult = rnorm(info[0], info[1]) / 1000.0;
    PutRNGstate();
  }
  else if (!strcmp(dist, dist1[3]))
    doseresult = r_trift(info[0], info[1],info[2]) / 1000.0;

  else if (!strcmp(dist, dist1[4]))
    doseresult = r_trift(info[0], exp((log(info[0])+log(info[2]))/2),info[2]) / 1000.0;

  else if (!strcmp(dist, dist1[5])) {
    GetRNGstate();
    doseresult = runif(info[0],info[1]) / 1000.0;
    PutRNGstate();
  }
  else
  {
    GetRNGstate();
    doseresult = exp(runif(log(info[0]),log(info[1]))) / 1000.0;
    PutRNGstate();
  }


  return doseresult;
}

// function to calculate latency
double latencyft(char* site)
{
  char *site1[2] = { "leukemia","thyroid" };
  double latencyresult;

  if (!strcmp(site, site1[0])) {

    while (1) {
      GetRNGstate();
      latencyresult = r_trift(2.0, 2.25, 2.5) - 0.401*log((1.0 / runif(0.0,1.0)) - 1.0);
      PutRNGstate();
      if ((0.5 <= latencyresult) & (latencyresult <= 4.1)) break;
    }

  }
  else if (!strcmp(site, site1[1])) {
    while (1) {
      GetRNGstate();
      latencyresult = r_trift(3.0, 5.0, 7.0) - 0.544*log((1.0 / runif(0.0, 1.0)) - 1.0);
      PutRNGstate();
      if ((2.5 <= latencyresult) & (latencyresult <= 7.6)) break;
    }
  }
  else {
    while (1) {
      GetRNGstate();
      latencyresult = r_trift(5.0, 7.5, 10.0) - 0.76*log((1.0 / runif(0.0, 1.0)) - 1.0);
      PutRNGstate();
      if ((4.0 <= latencyresult) & (latencyresult <= 11.0)) break;
    }
  }


  return latencyresult;
}

// function to determines the weight
double weightft(char* site)
{

  double weightresult = 0;

  for (int i = 0; i < 19; i++)
  {
    if (!strcmp(site, weightmat_site[i]))
    {
      weightresult = weightmat_value[i];
      break;
    }

  }
  return(weightresult);
}









// function to retrieve dbaseline data
void data_call_1()
{
  double aa[202] = { 0.00369,0.00032,0.00025,0.00018,0.00015,0.00013,0.00013,0.00013,0.00013,0.00012,0.00012,0.00012,0.00013,0.00016,0.0002,0.00026,0.00033,0.0004,0.00044,0.00047,0.0005,0.00053,0.00057,0.00061,0.00065,0.00069,0.00073,0.00076,0.00078,0.00081,0.00085,0.00088,0.00091,0.00094,0.001,0.00109,0.00118,0.00126,0.00136,0.00151,0.00173,0.00197,0.00221,0.00245,0.00268,0.00292,0.00318,0.0035,0.00385,0.00422,0.00457,0.00491,0.00532,0.00577,0.0062,0.00664,0.00715,0.00776,0.00849,0.00919,0.00986,0.01044,0.01113,0.01194,0.01294,0.01428,0.01591,0.01777,0.01992,0.02229,0.02501,0.02799,0.03131,0.03486,0.03861,0.04281,0.04777,0.05349,0.05989,0.06668,0.07376,0.08105,0.08957,0.09911,0.10958,0.12023,0.13161,0.14374,0.15663,0.17028,0.1847,0.19987,0.21579,0.23244,0.2498,0.26784,0.28652,0.3058,0.32563,0.34594,1,0.00275,0.0003,0.00022,0.00015,0.00011,0.00009,0.00008,0.00008,0.00008,0.00008,0.00008,0.00009,0.0001,0.00012,0.00014,0.00015,0.00016,0.00019,0.00023,0.00026,0.00029,0.00031,0.00034,0.00037,0.00041,0.00043,0.00045,0.00047,0.00049,0.00052,0.00053,0.00053,0.00053,0.00055,0.00058,0.00062,0.00067,0.00072,0.00076,0.0008,0.00083,0.00087,0.00092,0.00099,0.00108,0.00118,0.00128,0.00137,0.00146,0.00156,0.00166,0.00177,0.00189,0.00201,0.00211,0.00221,0.00235,0.00254,0.00278,0.00304,0.00336,0.00372,0.0041,0.00449,0.0049,0.00548,0.00624,0.00715,0.00814,0.00915,0.0103,0.01166,0.0134,0.01553,0.01793,0.02068,0.02374,0.02744,0.03179,0.0366,0.04193,0.04779,0.05494,0.06257,0.07052,0.07959,0.0895,0.10029,0.11198,0.12458,0.13811,0.15257,0.16793,0.18419,0.2013,0.21921,0.23787,0.2572,0.27711,0.2975,1 };

  for (int i = 0; i < 202; i++)
  {
    if (i < 101)
    {
      dbaseline[i].Prob_d_m = aa[i];
    }
    else
    {
      dbaseline[i-101].Prob_d_f = aa[i];

    }

  }

}

// function to retrieve dincidence data
void data_call_2()
{
  double aa[3838] = { 0,0.210238621,0.210238621,0.210238621,0.210238621,0.22534585,0.22534585,0.22534585,0.22534585,0.22534585,0.348141634,0.348141634,0.348141634,0.348141634,0.348141634,0.534510682,0.534510682,0.534510682,0.534510682,0.534510682,0.366120659,0.366120659,0.366120659,0.366120659,0.366120659,0.66600749,0.66600749,0.66600749,0.66600749,0.66600749,1.562054767,1.562054767,1.562054767,1.562054767,1.562054767,2.212185781,2.212185781,2.212185781,2.212185781,2.212185781,3.664313165,3.664313165,3.664313165,3.664313165,3.664313165,8.013510049,8.013510049,8.013510049,8.013510049,8.013510049,12.43040608,12.43040608,12.43040608,12.43040608,12.43040608,19.53185754,19.53185754,19.53185754,19.53185754,19.53185754,23.47095166,23.47095166,23.47095166,23.47095166,23.47095166,33.84922053,33.84922053,33.84922053,33.84922053,33.84922053,41.67671192,41.67671192,41.67671192,41.67671192,41.67671192,41.83925359,41.83925359,41.83925359,41.83925359,41.83925359,39.88632398,39.88632398,39.88632398,39.88632398,39.88632398,38.50279145,38.50279145,38.50279145,38.50279145,38.50279145,38.50279145,38.50279145,38.50279145,38.50279145,38.50279145,38.50279145,38.50279145,38.50279145,38.50279145,38.50279145,38.50279145,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0,0,0,0,0,0.9,0.9,0.9,0.9,0.9,2.7,2.7,2.7,2.7,2.7,8.1,8.1,8.1,8.1,8.1,18.2,18.2,18.2,18.2,18.2,31.1,31.1,31.1,31.1,31.1,48,48,48,48,48,61.5,61.5,61.5,61.5,61.5,66.1,66.1,66.1,66.1,66.1,63.9,63.9,63.9,63.9,63.9,48.4,48.4,48.4,48.4,48.4,48.4,48.4,48.4,48.4,48.4,48.4,48.4,48.4,48.4,48.4,48.4,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,0.3,0.3,0.3,0.3,0.3,2.4,2.4,2.4,2.4,2.4,7,7,7,7,7,16.5,16.5,16.5,16.5,16.5,39.9,39.9,39.9,39.9,39.9,68.9,68.9,68.9,68.9,68.9,122.6,122.6,122.6,122.6,122.6,191.1,191.1,191.1,191.1,191.1,274.2,274.2,274.2,274.2,274.2,373,373,373,373,373,455.9,455.9,455.9,455.9,455.9,530.7,530.7,530.7,530.7,530.7,543.7,543.7,543.7,543.7,543.7,438.9,438.9,438.9,438.9,438.9,438.9,438.9,438.9,438.9,438.9,438.9,438.9,438.9,438.9,438.9,438.9,
                      0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.3,0.3,0.3,0.3,0.3,0.7,0.7,0.7,0.7,0.7,1.3,1.3,1.3,1.3,1.3,2.7,2.7,2.7,2.7,2.7,5.6,5.6,5.6,5.6,5.6,11.4,11.4,11.4,11.4,11.4,22.8,22.8,22.8,22.8,22.8,47.9,47.9,47.9,47.9,47.9,75.7,75.7,75.7,75.7,75.7,120.5,120.5,120.5,120.5,120.5,171.2,171.2,171.2,171.2,171.2,210.7,210.7,210.7,210.7,210.7,237,237,237,237,237,231.7,231.7,231.7,231.7,231.7,253,253,253,253,253,253,253,253,253,253,253,253,253,253,253,253,
                      0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.4,0.4,0.4,0.4,0.4,1.6,1.6,1.6,1.6,1.6,3.9,3.9,3.9,3.9,3.9,8.9,8.9,8.9,8.9,8.9,13.4,13.4,13.4,13.4,13.4,23.9,23.9,23.9,23.9,23.9,45.8,45.8,45.8,45.8,45.8,76.8,76.8,76.8,76.8,76.8,105.9,105.9,105.9,105.9,105.9,125.8,125.8,125.8,125.8,125.8,164.6,164.6,164.6,164.6,164.6,166.3,166.3,166.3,166.3,166.3,171.3,171.3,171.3,171.3,171.3,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.3,0.3,0.3,0.3,0.3,1,1,1,1,1,1.4,1.4,1.4,1.4,1.4,3.7,3.7,3.7,3.7,3.7,9.1,9.1,9.1,9.1,9.1,17.2,17.2,17.2,17.2,17.2,33.8,33.8,33.8,33.8,33.8,49.7,49.7,49.7,49.7,49.7,81.5,81.5,81.5,81.5,81.5,97.7,97.7,97.7,97.7,97.7,123.2,123.2,123.2,123.2,123.2,122.1,122.1,122.1,122.1,122.1,122.1,122.1,122.1,122.1,122.1,122.1,122.1,122.1,122.1,122.1,122.1,
                      0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.9,0.9,0.9,0.9,0.9,2.3,2.3,2.3,2.3,2.3,6.5,6.5,6.5,6.5,6.5,11.4,11.4,11.4,11.4,11.4,19.1,19.1,19.1,19.1,19.1,32.4,32.4,32.4,32.4,32.4,50.7,50.7,50.7,50.7,50.7,73.9,73.9,73.9,73.9,73.9,86.6,86.6,86.6,86.6,86.6,87.4,87.4,87.4,87.4,87.4,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,
                      0.4,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,0.2,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,1.4,1.4,1.4,1.4,1.4,3.4,3.4,3.4,3.4,3.4,9.1,9.1,9.1,9.1,9.1,26.5,26.5,26.5,26.5,26.5,58.2,58.2,58.2,58.2,58.2,96,96,96,96,96,136.5,136.5,136.5,136.5,136.5,158.4,158.4,158.4,158.4,158.4,178.5,178.5,178.5,178.5,178.5,215.2,215.2,215.2,215.2,215.2,220.5,220.5,220.5,220.5,220.5,245.8,245.8,245.8,245.8,245.8,238.7,238.7,238.7,238.7,238.7,238.7,238.7,238.7,238.7,238.7,238.7,238.7,238.7,238.7,238.7,238.7,
                      0,0.1,0.1,0.1,0.1,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.9,0.9,0.9,0.9,0.9,1.9,1.9,1.9,1.9,1.9,3.7,3.7,3.7,3.7,3.7,7.7,7.7,7.7,7.7,7.7,18.6,18.6,18.6,18.6,18.6,44.9,44.9,44.9,44.9,44.9,98.4,98.4,98.4,98.4,98.4,179.7,179.7,179.7,179.7,179.7,315.5,315.5,315.5,315.5,315.5,472.6,472.6,472.6,472.6,472.6,667.8,667.8,667.8,667.8,667.8,670.4,670.4,670.4,670.4,670.4,650.1,650.1,650.1,650.1,650.1,650.1,650.1,650.1,650.1,650.1,650.1,650.1,650.1,650.1,650.1,650.1,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.4,0.4,0.4,0.4,0.4,0.9,0.9,0.9,0.9,0.9,0.4,0.4,0.4,0.4,0.4,0.6,0.6,0.6,0.6,0.6,1.8,1.8,1.8,1.8,1.8,2.2,2.2,2.2,2.2,2.2,0.6,0.6,0.6,0.6,0.6,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0.4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,1,1,1,1,1,3.1,3.1,3.1,3.1,3.1,12.6,12.6,12.6,12.6,12.6,38.3,38.3,38.3,38.3,38.3,98.7,98.7,98.7,98.7,98.7,213.8,213.8,213.8,213.8,213.8,298,298,298,298,298,374.1,374.1,374.1,374.1,374.1,356,356,356,356,356,377.3,377.3,377.3,377.3,377.3,377.3,377.3,377.3,377.3,377.3,377.3,377.3,377.3,377.3,377.3,377.3,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,0.6,0.6,0.6,0.6,1.6,1.6,1.6,1.6,1.6,2.7,2.7,2.7,2.7,2.7,4.4,4.4,4.4,4.4,4.4,9.9,9.9,9.9,9.9,9.9,19.3,19.3,19.3,19.3,19.3,31,31,31,31,31,53.9,53.9,53.9,53.9,53.9,82.9,82.9,82.9,82.9,82.9,111.5,111.5,111.5,111.5,111.5,149,149,149,149,149,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,134.2,
                      3.1,1.6,1.6,1.6,1.6,0.1,0.1,0.1,0.1,0.1,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.4,0.4,0.4,0.4,0.4,1,1,1,1,1,3,3,3,3,3,4.8,4.8,4.8,4.8,4.8,8.1,8.1,8.1,8.1,8.1,11.6,11.6,11.6,11.6,11.6,17.7,17.7,17.7,17.7,17.7,25.7,25.7,25.7,25.7,25.7,30.8,30.8,30.8,30.8,30.8,37.6,37.6,37.6,37.6,37.6,45.2,45.2,45.2,45.2,45.2,43.5,43.5,43.5,43.5,43.5,38.1,38.1,38.1,38.1,38.1,37.4,37.4,37.4,37.4,37.4,37.4,37.4,37.4,37.4,37.4,37.4,37.4,37.4,37.4,37.4,37.4,
                      5.3,2.2,2.2,2.2,2.2,2.6,2.6,2.6,2.6,2.6,2.4,2.4,2.4,2.4,2.4,2.2,2.2,2.2,2.2,2.2,1.4,1.4,1.4,1.4,1.4,1.6,1.6,1.6,1.6,1.6,2,2,2,2,2,2.9,2.9,2.9,2.9,2.9,2.8,2.8,2.8,2.8,2.8,3.6,3.6,3.6,3.6,3.6,4.8,4.8,4.8,4.8,4.8,6.7,6.7,6.7,6.7,6.7,6.5,6.5,6.5,6.5,6.5,9.8,9.8,9.8,9.8,9.8,10.8,10.8,10.8,10.8,10.8,15.4,15.4,15.4,15.4,15.4,17,17,17,17,17,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,
                      0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.3,0.3,0.3,0.3,0.3,1.6,1.6,1.6,1.6,1.6,4.1,4.1,4.1,4.1,4.1,12.9,12.9,12.9,12.9,12.9,25.1,25.1,25.1,25.1,25.1,37.1,37.1,37.1,37.1,37.1,43.3,43.3,43.3,43.3,43.3,43.7,43.7,43.7,43.7,43.7,49.7,49.7,49.7,49.7,49.7,45.9,45.9,45.9,45.9,45.9,43.7,43.7,43.7,43.7,43.7,39.1,39.1,39.1,39.1,39.1,28.7,28.7,28.7,28.7,28.7,18.9,18.9,18.9,18.9,18.9,12.9,12.9,12.9,12.9,12.9,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
                      20.87682672,8.619783454,8.619783454,8.619783454,8.619783454,3.981110009,3.981110009,3.981110009,3.981110009,3.981110009,5.918407786,5.918407786,5.918407786,5.918407786,5.918407786,9.033230529,9.033230529,9.033230529,9.033230529,9.033230529,8.908936029,8.908936029,8.908936029,8.908936029,8.908936029,10.09257504,10.09257504,10.09257504,10.09257504,10.09257504,13.1011045,13.1011045,13.1011045,13.1011045,13.1011045,15.88349391,15.88349391,15.88349391,15.88349391,15.88349391,18.27687908,18.27687908,18.27687908,18.27687908,18.27687908,29.59534961,29.59534961,29.59534961,29.59534961,29.59534961,49.6712988,49.6712988,49.6712988,49.6712988,49.6712988,73.13595547,73.13595547,73.13595547,73.13595547,73.13595547,108.1908011,108.1908011,108.1908011,108.1908011,108.1908011,161.8380241,161.8380241,161.8380241,161.8380241,161.8380241,218.2703708,218.2703708,218.2703708,218.2703708,218.2703708,277.2187964,277.2187964,277.2187964,277.2187964,277.2187964,302.0802478,302.0802478,302.0802478,302.0802478,302.0802478,322.3233684,322.3233684,322.3233684,322.3233684,322.3233684,322.3233684,322.3233684,322.3233684,322.3233684,322.3233684,322.3233684,322.3233684,322.3233684,322.3233684,322.3233684,322.3233684,
                      8.439568249,7.042993798,7.042993798,7.042993798,7.042993798,4.056225292,4.056225292,4.056225292,4.056225292,4.056225292,4.293746825,4.293746825,4.293746825,4.293746825,4.293746825,3.634672639,3.634672639,3.634672639,3.634672639,3.634672639,2.86794516,2.86794516,2.86794516,2.86794516,2.86794516,2.612798615,2.612798615,2.612798615,2.612798615,2.612798615,3.325664988,3.325664988,3.325664988,3.325664988,3.325664988,3.141303809,3.141303809,3.141303809,3.141303809,3.141303809,3.083385468,3.083385468,3.083385468,3.083385468,3.083385468,5.691413387,5.691413387,5.691413387,5.691413387,5.691413387,6.190040276,6.190040276,6.190040276,6.190040276,6.190040276,9.404227707,9.404227707,9.404227707,9.404227707,9.404227707,12.24977955,12.24977955,12.24977955,12.24977955,12.24977955,18.12129988,18.12129988,18.12129988,18.12129988,18.12129988,22.20729905,22.20729905,22.20729905,22.20729905,22.20729905,26.72313616,26.72313616,26.72313616,26.72313616,26.72313616,24.0491071,24.0491071,24.0491071,24.0491071,24.0491071,29.70215341,29.70215341,29.70215341,29.70215341,29.70215341,29.70215341,29.70215341,29.70215341,29.70215341,29.70215341,29.70215341,29.70215341,29.70215341,29.70215341,29.70215341,29.70215341,
                      0,0.111932196,0.111932196,0.111932196,0.111932196,0.163128478,0.163128478,0.163128478,0.163128478,0.163128478,0.19142041,0.19142041,0.19142041,0.19142041,0.19142041,0.363973977,0.363973977,0.363973977,0.363973977,0.363973977,0.537273338,0.537273338,0.537273338,0.537273338,0.537273338,0.974144843,0.974144843,0.974144843,0.974144843,0.974144843,1.467065429,1.467065429,1.467065429,1.467065429,1.467065429,1.888153609,1.888153609,1.888153609,1.888153609,1.888153609,2.354838628,2.354838628,2.354838628,2.354838628,2.354838628,4.111931502,4.111931502,4.111931502,4.111931502,4.111931502,4.800010621,4.800010621,4.800010621,4.800010621,4.800010621,6.147395238,6.147395238,6.147395238,6.147395238,6.147395238,5.897809723,5.897809723,5.897809723,5.897809723,5.897809723,8.375229771,8.375229771,8.375229771,8.375229771,8.375229771,10.57284348,10.57284348,10.57284348,10.57284348,10.57284348,10.96764862,10.96764862,10.96764862,10.96764862,10.96764862,10.96730307,10.96730307,10.96730307,10.96730307,10.96730307,10.64887797,10.64887797,10.64887797,10.64887797,10.64887797,10.64887797,10.64887797,10.64887797,10.64887797,10.64887797,10.64887797,10.64887797,10.64887797,10.64887797,10.64887797,10.64887797,
                      0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.4,0.4,0.4,0.4,0.4,0.5,0.5,0.5,0.5,0.5,1.5,1.5,1.5,1.5,1.5,1.6,1.6,1.6,1.6,1.6,1.9,1.9,1.9,1.9,1.9,2.8,2.8,2.8,2.8,2.8,5.4,5.4,5.4,5.4,5.4,5.2,5.2,5.2,5.2,5.2,9.1,9.1,9.1,9.1,9.1,9.1,9.1,9.1,9.1,9.1,9.1,9.1,9.1,9.1,9.1,9.1,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,1.3,1.3,1.3,1.3,1.3,4.5,4.5,4.5,4.5,4.5,10.4,10.4,10.4,10.4,10.4,18.4,18.4,18.4,18.4,18.4,29.7,29.7,29.7,29.7,29.7,34.2,34.2,34.2,34.2,34.2,48.5,48.5,48.5,48.5,48.5,63.1,63.1,63.1,63.1,63.1,83.6,83.6,83.6,83.6,83.6,126.2,126.2,126.2,126.2,126.2,164.4,164.4,164.4,164.4,164.4,199.5,199.5,199.5,199.5,199.5,210.5,210.5,210.5,210.5,210.5,187.9,187.9,187.9,187.9,187.9,187.9,187.9,187.9,187.9,187.9,187.9,187.9,187.9,187.9,187.9,187.9,
                      0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.7,0.7,0.7,0.7,0.7,2.5,2.5,2.5,2.5,2.5,5.2,5.2,5.2,5.2,5.2,9.1,9.1,9.1,9.1,9.1,17.1,17.1,17.1,17.1,17.1,29.8,29.8,29.8,29.8,29.8,44.8,44.8,44.8,44.8,44.8,61.7,61.7,61.7,61.7,61.7,82.6,82.6,82.6,82.6,82.6,106.2,106.2,106.2,106.2,106.2,131.8,131.8,131.8,131.8,131.8,134,134,134,134,134,123.2,123.2,123.2,123.2,123.2,123.2,123.2,123.2,123.2,123.2,123.2,123.2,123.2,123.2,123.2,123.2,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.3,0.3,0.3,0.3,0.3,1.2,1.2,1.2,1.2,1.2,3.1,3.1,3.1,3.1,3.1,5.1,5.1,5.1,5.1,5.1,9.2,9.2,9.2,9.2,9.2,14.4,14.4,14.4,14.4,14.4,25.1,25.1,25.1,25.1,25.1,31.2,31.2,31.2,31.2,31.2,41.8,41.8,41.8,41.8,41.8,56.8,56.8,56.8,56.8,56.8,72,72,72,72,72,81.4,81.4,81.4,81.4,81.4,94.5,94.5,94.5,94.5,94.5,79.1,79.1,79.1,79.1,79.1,79.1,79.1,79.1,79.1,79.1,79.1,79.1,79.1,79.1,79.1,79.1,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.6,1.5,1.5,1.5,1.5,1.5,3.7,3.7,3.7,3.7,3.7,4.5,4.5,4.5,4.5,4.5,12.4,12.4,12.4,12.4,12.4,21.3,21.3,21.3,21.3,21.3,31.7,31.7,31.7,31.7,31.7,49.8,49.8,49.8,49.8,49.8,63.1,63.1,63.1,63.1,63.1,89.3,89.3,89.3,89.3,89.3,98.5,98.5,98.5,98.5,98.5,98.5,98.5,98.5,98.5,98.5,98.5,98.5,98.5,98.5,98.5,98.5,
                      0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.1,0.1,0.1,0.1,0.1,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.9,0.9,0.9,0.9,0.9,1.2,1.2,1.2,1.2,1.2,3.1,3.1,3.1,3.1,3.1,5.9,5.9,5.9,5.9,5.9,11.2,11.2,11.2,11.2,11.2,17.2,17.2,17.2,17.2,17.2,32,32,32,32,32,45,45,45,45,45,60.1,60.1,60.1,60.1,60.1,76.5,76.5,76.5,76.5,76.5,71.9,71.9,71.9,71.9,71.9,71.9,71.9,71.9,71.9,71.9,71.9,71.9,71.9,71.9,71.9,71.9,
                      1.4,0.6,0.6,0.6,0.6,0.2,0.2,0.2,0.2,0.2,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.3,0.3,0.3,0.3,0.3,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,1.9,1.9,1.9,1.9,1.9,5.2,5.2,5.2,5.2,5.2,10.6,10.6,10.6,10.6,10.6,18.1,18.1,18.1,18.1,18.1,29.7,29.7,29.7,29.7,29.7,42.9,42.9,42.9,42.9,42.9,61.5,61.5,61.5,61.5,61.5,75.1,75.1,75.1,75.1,75.1,96,96,96,96,96,104.2,104.2,104.2,104.2,104.2,89,89,89,89,89,89,89,89,89,89,89,89,89,89,89,89,
                      0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.3,0.3,0.3,0.3,0.3,0.7,0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.9,3.4,3.4,3.4,3.4,3.4,7,7,7,7,7,14.4,14.4,14.4,14.4,14.4,24.1,24.1,24.1,24.1,24.1,40.2,40.2,40.2,40.2,40.2,58.4,58.4,58.4,58.4,58.4,73.7,73.7,73.7,73.7,73.7,110.3,110.3,110.3,110.3,110.3,160.2,160.2,160.2,160.2,160.2,192.2,192.2,192.2,192.2,192.2,197.4,197.4,197.4,197.4,197.4,197.4,197.4,197.4,197.4,197.4,197.4,197.4,197.4,197.4,197.4,197.4,
                      0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,1.9,1.9,1.9,1.9,1.9,8.5,8.5,8.5,8.5,8.5,27.2,27.2,27.2,27.2,27.2,62.1,62.1,62.1,62.1,62.1,104.2,104.2,104.2,104.2,104.2,145.1,145.1,145.1,145.1,145.1,131.2,131.2,131.2,131.2,131.2,111.2,111.2,111.2,111.2,111.2,102.1,102.1,102.1,102.1,102.1,83.7,83.7,83.7,83.7,83.7,65.3,65.3,65.3,65.3,65.3,49.4,49.4,49.4,49.4,49.4,33.7,33.7,33.7,33.7,33.7,25.5,25.5,25.5,25.5,25.5,25.5,25.5,25.5,25.5,25.5,25.5,25.5,25.5,25.5,25.5,25.5,
                      0,0,0,0,0,0.5,0.5,0.5,0.5,0.5,1.7,1.7,1.7,1.7,1.7,3,3,3,3,3,2.2,2.2,2.2,2.2,2.2,2.5,2.5,2.5,2.5,2.5,4.2,4.2,4.2,4.2,4.2,5.1,5.1,5.1,5.1,5.1,9,9,9,9,9,12.7,12.7,12.7,12.7,12.7,15,15,15,15,15,16.4,16.4,16.4,16.4,16.4,17.2,17.2,17.2,17.2,17.2,15.5,15.5,15.5,15.5,15.5,16.1,16.1,16.1,16.1,16.1,16.5,16.5,16.5,16.5,16.5,14.1,14.1,14.1,14.1,14.1,19.8,19.8,19.8,19.8,19.8,19.8,19.8,19.8,19.8,19.8,19.8,19.8,19.8,19.8,19.8,19.8,
                      0,0,0,0,0,0,0,0,0,0,0.127613606,0.127613606,0.127613606,0.127613606,0.127613606,0.181986989,0.181986989,0.181986989,0.181986989,0.181986989,2.149093351,2.149093351,2.149093351,2.149093351,2.149093351,8.117873691,8.117873691,8.117873691,8.117873691,8.117873691,16.39969569,16.39969569,16.39969569,16.39969569,16.39969569,21.18416244,21.18416244,21.18416244,21.18416244,21.18416244,32.41366111,32.41366111,32.41366111,32.41366111,32.41366111,41.07205144,41.07205144,41.07205144,41.07205144,41.07205144,45.44690907,45.44690907,45.44690907,45.44690907,45.44690907,46.46287098,46.46287098,46.46287098,46.46287098,46.46287098,45.48462347,45.48462347,45.48462347,45.48462347,45.48462347,42.16830803,42.16830803,42.16830803,42.16830803,42.16830803,40.13133063,40.13133063,40.13133063,40.13133063,40.13133063,42.4400316,42.4400316,42.4400316,42.4400316,42.4400316,44.65259109,44.65259109,44.65259109,44.65259109,44.65259109,41.0742436,41.0742436,41.0742436,41.0742436,41.0742436,41.0742436,41.0742436,41.0742436,41.0742436,41.0742436,41.0742436,41.0742436,41.0742436,41.0742436,41.0742436,41.0742436,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.3,0.3,0.3,0.3,0.3,0.7,0.7,0.7,0.7,0.7,1.8,1.8,1.8,1.8,1.8,2.9,2.9,2.9,2.9,2.9,4.3,4.3,4.3,4.3,4.3,8,8,8,8,8,13.8,13.8,13.8,13.8,13.8,19.4,19.4,19.4,19.4,19.4,31.1,31.1,31.1,31.1,31.1,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,
                      3.3,0.7,0.7,0.7,0.7,0.2,0.2,0.2,0.2,0.2,0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,1.3,1.3,1.3,1.3,1.3,2.2,2.2,2.2,2.2,2.2,2.5,2.5,2.5,2.5,2.5,5,5,5,5,5,7.5,7.5,7.5,7.5,7.5,8.4,8.4,8.4,8.4,8.4,11.3,11.3,11.3,11.3,11.3,13.4,13.4,13.4,13.4,13.4,14.6,14.6,14.6,14.6,14.6,17.2,17.2,17.2,17.2,17.2,10.7,10.7,10.7,10.7,10.7,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
                      1.9,2,2,2,2,2.4,2.4,2.4,2.4,2.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.1,1.1,1.1,1.1,1.1,1.7,1.7,1.7,1.7,1.7,1.5,1.5,1.5,1.5,1.5,2.3,2.3,2.3,2.3,2.3,2.3,2.3,2.3,2.3,2.3,2.6,2.6,2.6,2.6,2.6,4.1,4.1,4.1,4.1,4.1,4.9,4.9,4.9,4.9,4.9,6,6,6,6,6,6,6,6,6,6,9.4,9.4,9.4,9.4,9.4,11.1,11.1,11.1,11.1,11.1,11.5,11.5,11.5,11.5,11.5,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
                      0,0,0,0,0,0.5,0.5,0.5,0.5,0.5,1.5,1.5,1.5,1.5,1.5,7.5,7.5,7.5,7.5,7.5,26.8,26.8,26.8,26.8,26.8,70.1,70.1,70.1,70.1,70.1,124.5,124.5,124.5,124.5,124.5,164.1,164.1,164.1,164.1,164.1,198.8,198.8,198.8,198.8,198.8,236.5,236.5,236.5,236.5,236.5,259,259,259,259,259,227.7,227.7,227.7,227.7,227.7,200.5,200.5,200.5,200.5,200.5,147.2,147.2,147.2,147.2,147.2,89,89,89,89,89,52.5,52.5,52.5,52.5,52.5,32.9,32.9,32.9,32.9,32.9,12.9,12.9,12.9,12.9,12.9,12.9,12.9,12.9,12.9,12.9,12.9,12.9,12.9,12.9,12.9,12.9,
                      17.50468369,7.499457129,7.499457129,7.499457129,7.499457129,2.446927169,2.446927169,2.446927169,2.446927169,2.446927169,3.828408193,3.828408193,3.828408193,3.828408193,3.828408193,4.731661702,4.731661702,4.731661702,4.731661702,4.731661702,6.111484218,6.111484218,6.111484218,6.111484218,6.111484218,9.20025685,9.20025685,9.20025685,9.20025685,9.20025685,10.74101475,10.74101475,10.74101475,10.74101475,10.74101475,12.01970956,12.01970956,12.01970956,12.01970956,12.01970956,16.43769709,16.43769709,16.43769709,16.43769709,16.43769709,26.37307791,26.37307791,26.37307791,26.37307791,26.37307791,35.08092869,35.08092869,35.08092869,35.08092869,35.08092869,45.6050949,45.6050949,45.6050949,45.6050949,45.6050949,67.2886473,67.2886473,67.2886473,67.2886473,67.2886473,85.11570721,85.11570721,85.11570721,85.11570721,85.11570721,128.0109867,128.0109867,128.0109867,128.0109867,128.0109867,167.2169035,167.2169035,167.2169035,167.2169035,167.2169035,193.7556877,193.7556877,193.7556877,193.7556877,193.7556877,225.5280227,225.5280227,225.5280227,225.5280227,225.5280227,225.5280227,225.5280227,225.5280227,225.5280227,225.5280227,225.5280227,225.5280227,225.5280227,225.5280227,225.5280227,225.5280227,
                      3.784796473,6.492067365,6.492067365,6.492067365,6.492067365,2.936312603,2.936312603,2.936312603,2.936312603,2.936312603,2.552272128,2.552272128,2.552272128,2.552272128,2.552272128,1.880532215,1.880532215,1.880532215,1.880532215,1.880532215,2.283411686,2.283411686,2.283411686,2.283411686,2.283411686,2.164766318,2.164766318,2.164766318,2.164766318,2.164766318,2.095807756,2.095807756,2.095807756,2.095807756,2.095807756,3.730254691,3.730254691,3.730254691,3.730254691,3.730254691,3.416824676,3.416824676,3.416824676,3.416824676,3.416824676,4.726358048,4.726358048,4.726358048,4.726358048,4.726358048,4.748946679,4.748946679,4.748946679,4.748946679,4.748946679,6.576283278,6.576283278,6.576283278,6.576283278,6.576283278,6.255252737,6.255252737,6.255252737,6.255252737,6.255252737,8.862161734,8.862161734,8.862161734,8.862161734,8.862161734,12.50551379,12.50551379,12.50551379,12.50551379,12.50551379,13.82877434,13.82877434,13.82877434,13.82877434,13.82877434,16.97320714,16.97320714,16.97320714,16.97320714,16.97320714,11.40951211,11.40951211,11.40951211,11.40951211,11.40951211,11.40951211,11.40951211,11.40951211,11.40951211,11.40951211,11.40951211,11.40951211,11.40951211,11.40951211,11.40951211,11.40951211

  };

  for (int i = 0; i < 3838; i++)
  {
    if(i<1919)
      dincidence[i].Rate_m = aa[i];
    else
      dincidence[i-1919].Rate_f = aa[i];

  }
}

// function to retrieve NEW dbaseline data
void data_call_3(double * data_baseline)
{
  for (int i = 0; i < 202; i++)
  {
    if (i < 101)
    {
      dbaseline[i].Prob_d_m = data_baseline[i];
    }
    else
    {
      dbaseline[i - 101].Prob_d_f = data_baseline[i];

    }

  }

}
// function to retrieve NEW dincidence
void data_call_4(double * data_incidence)
{
  for (int i = 0; i < 3838; i++)
  {
    if (i < 1919)
      dincidence[i].Rate_m = data_incidence[i];
    else
      dincidence[i - 1919].Rate_f = data_incidence[i];

  }
}



// function to calculate surv
double survft(int age, char* sex)
{
  double product = 1.0;

  if (!strcmp(sex, "male")) {
    for (int i = 0; i < age + 1; i++)
    {
      product *= 1.0 - dbaseline[i].Prob_d_m;
    }
  }
  else {
    for (int i = 0; i < age + 1; i++)
    {
      product *= 1.0 - dbaseline[i].Prob_d_f;
    }

  }
  return product;
}


// function to calculate srate
double srateft(int e_age, int a_age, char* sex)
{
  double product = 1.0;

  if (!strcmp(sex, "male")) {
    for (int i = e_age + 1; i < a_age + 1; i++)
    {
      product *= 1.0 - dbaseline[i].Prob_d_m;
    }
  }
  else {
    for (int i = e_age + 1; i < a_age + 1; i++)
    {
      product *= 1.0 - dbaseline[i].Prob_d_f;
    }
  }
  return product;
}

// function to calculate incidence
double incidenceft(int age, char* sex, char* site)
{
  double incidenceresult = 1.0;
  int site_tmp = 1;
  char *site1[19] = { "oral","oesophagus","stomach","colon","rectum","gallbladder","pancreas","liver","lung","breast",
                      "ovary","uterus","prostate","bladder","kidney","brain/cns","thyroid","remainder","leukemia" };

  for (int i = 0; i < 19; i++)		// Find the right location for the site
  {
    if (!strcmp(site, site1[i]))
    {
      site_tmp = i;
      break;
    }

  }

  if (!strcmp(sex, "male"))		// male: 3, female: 4
    incidenceresult = dincidence[(site_tmp * 101 + (age))].Rate_m;
  else
    incidenceresult = dincidence[(site_tmp * 101 + (age))].Rate_f;

  return incidenceresult;
}



// function to calculate brft baserisk
void brft(int* sex1, int* site1, int *age, int* changedata, double* dbaseline, double* dincidence , double* out)
{
  if(*changedata==2)
  {
    data_call_1();			// function to retrieve dbaseline data > Run before survft and srateft
    data_call_2();
  }
  else
  {
    data_call_3(dbaseline);
    data_call_4(dincidence);
  }

  char * sex, *site;

  if (*sex1 == 1)
    sex = "male";
  else
    sex = "female";


  switch (*site1)
  {
  case 1:
    site = "stomach";
    break;
  case 2:
    site = "colon";
    break;
  case 3:
    site = "liver";
    break;
  case 4:
    site = "lung";
    break;
  case 5:
    site = "breast";
    break;
  case 6:
    site = "ovary";
    break;
  case 7:
    site = "uterus";
    break;
  case 8:
    site = "prostate";
    break;
  case 9:
    site = "bladder";
    break;
  case 10:
    site = "brain/cns";
    break;
  case 11:
    site = "thyroid";
    break;
  case 12:
    site = "remainder";
    break;
  case 13:
    site = "oral";
    break;
  case 14:
    site = "oesophagus";
    break;
  case 15:
    site = "rectum";
    break;
  case 16:
    site = "gallbladder";
    break;
  case 17:
    site = "pancreas";
    break;
  case 18:
    site = "kidney";
    break;
  case 19:
    site = "leukemia";
    break;
  default:
    site = "";
  }

  if (100 <= *age)
    *age = 100;

  for (int i = *age; i < 101; i++)
  {
    *out += survft(i, sex)*incidenceft(i, sex, site);
  }
}


// function to calculate vec * mat * vex
double vmvft(double* vec, double(*mat)[5], int column)
{
  double prod_tmp1, prod_tmp2 = 0;
  for (int a = 0; a < 5; a++)
  {
    prod_tmp1 = 0;
    for (int b = 0; b < 5; b++)
    {
      prod_tmp1 += vec[b] * mat[a][b];

    }
    prod_tmp2 += prod_tmp1 * vec[a];
  }
  return prod_tmp2;
}

// function to return DDREF
double DDREFft(char * DDREF_op, char * exposure_rate, double dose)
{
  double DDREF_tmp, Dlim;
  if (!strcmp(DDREF_op, "F")) {
    DDREF_tmp = (double)1;
  }
  else {
    if (!strcmp(exposure_rate, "chronic")) {
      GetRNGstate();
      DDREF_tmp = rlnorm(log(1.5), log(1.35));
      PutRNGstate();
    }
    else {

      GetRNGstate();
      Dlim = exp(runif(0.03,0.2));
      PutRNGstate();
      //Dlim = rlunif(1, 0.03, 0.2, base = exp(1))
      if (dose < Dlim) {
        GetRNGstate();
        DDREF_tmp = rlnorm(log(1.5), log(1.35)) ;
        PutRNGstate();
      }
      else {
        DDREF_tmp = (double)1;
      }
    }
  }
  return DDREF_tmp;
}

// function to calculate approx
double approxft(double a, double b, double c) {

  double r_unif = 0.0;
  double result = 0.0;

  GetRNGstate();
  r_unif = runif(0.025, 0.975);
  PutRNGstate();


  if (r_unif <= 0.5) {
    result = ((r_unif - 0.025)*(b - a) / 0.475) + a;
  }
  else {
    result = ((r_unif - 0.5)*(c - b) / 0.475) + b;
  }

  return result;

}

// function to generate random number from triangulardistribution
double r_trift(double a, double b, double c)
{
  double F = (c - a) / (b - a);
  double r_unif;
  double result;

  GetRNGstate();
  r_unif = runif(0.0, 1.0);
  PutRNGstate();

  if ((0 < r_unif) & (r_unif < F)) {
    result = a + sqrt(r_unif*(b - a)*(c - a));
  }
  else {
    result = b - sqrt((1 - r_unif)*(b - a)*(b - a));
  }
  return result;
}



// RegisteringDynamic Symbols
void R_init_LARisk(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
