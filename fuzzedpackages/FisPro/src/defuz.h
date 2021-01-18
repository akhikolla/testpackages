//***********************************************************************
//
//
//                              DEFUZ.H
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC 
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  July 31, 2009 - Contact : fispro@supagro.inra.fr
// File : DEFUZ class functions used by FISPRO, part of library fispro

//**********************************************************************

#define NOTHING             0
#define NO_ACTIVE_RULE      1  
#define AMBIGUITY           2
#define NON_CONNEX_AREA     3  
#define NON_CONNEX_SEGMENT  4
#define EMPTY_OUTPUT        5
#define AMBIGU_THRES        0.1
#define MIN_THRES           0.1
#define EQUALITY_THRES      0.1 

//**********************************************************************

//                           OUT_CRISP

//**********************************************************************


class DEFUZ_Sugeno : public DEFUZ
//*******************************
{
  protected :


    public :

    DEFUZ_Sugeno() : DEFUZ()
    {
      Thres = 0.;      // Not used
    }

  virtual ~DEFUZ_Sugeno() {}

  virtual double EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display);

  virtual void WriteHeader(FILE *p, FISOUT *O) const
    {
      fprintf(p, "     %s", INFERRED);
      fprintf(p, "    %s", ALARM);
    }
};

class DEFUZ_SugenoClassif : public DEFUZ_Sugeno
//*********************************************
{
  protected :
  
    double * Classes;  // labels of the different classes 

  public :

    DEFUZ_SugenoClassif();

  virtual ~DEFUZ_SugenoClassif()
    {
      if(Classes) delete [] Classes;
    }

  //! 'c' is the array of class labels (double), 'n' its size
  void SetClasses(double *c, int n)
    {
      if(Classes) delete [] Classes;

      NbClas = n;
      if(NbClas == 0) return;

      Classes = new double [NbClas];
      for(int i = 0; i < NbClas; i++) Classes[i] = c[i];
    }

  //! 'T' is the array of output observed values in data file, 'n' its size
  //! the Classes array is allocated inside InitUniq
  void InitClasses(double *T, int n)
    {
      if(Classes) delete [] Classes;
      Classes = NULL;
          
      InitUniq(T, n, Classes, NbClas);
    }      

  double * ClasLabel(void) { return Classes; } 

  //! Computes an inferred value using the sugeno defuzzification operator 
  //! and returns the closest class label.
  virtual double EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display);

  virtual void WriteHeader(FILE *p, FISOUT *O) const
    {
      DEFUZ_Sugeno::WriteHeader(p, O);
      fprintf(p, "    %s", CLASINF);
      fprintf(p, "    %s", CLASALARM);
    }
};

//**********************************************************************

//                           OUT_FUZZY

//**********************************************************************

class DEFUZ_SugenoFuzzy : public DEFUZ
//************************************
// Output is not calculated with the rule conclusion, 
// which is the MF number, but with the kernel center 
// of this MF.
{
  protected :
    double * Consequences;

  public :

    DEFUZ_SugenoFuzzy();

  virtual ~DEFUZ_SugenoFuzzy() 
    {
      delete [] Consequences;
    }

  void InitConsequences(FISOUT * O);

  //! Each rule is associated a crisp consequency, computed as the middle of 
  //! the MF kernel, and then the sugeno defuzzification operator is applied. 
  virtual double EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display);

  virtual void WriteHeader(FILE *p, FISOUT *O) const;

};

class DEFUZ_WeArea : public DEFUZ
//*******************************
{
  public :

    DEFUZ_WeArea();
  virtual ~DEFUZ_WeArea() {}
  
  //! Close to centroid defuzzification. As the centroids are computed at the MF level, 
  //! the common areas are taken into account twice, one time with each MF. 
  virtual double EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display);

  virtual void WriteHeader(FILE *p, FISOUT *O) const;
};

class DEFUZ_MeanMax : public DEFUZ
//********************************
{
  public :

  //! 'v' is a threshold. Two degrees whose difference is lower than 'v' will be 
  //! considered as max. The conflict management is depending on the the strategy 's'.  
  DEFUZ_MeanMax(void);

  virtual ~DEFUZ_MeanMax() {}

  //! Returns the mean of the maximum. 
  virtual double EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display);

   virtual void WriteHeader(FILE *p, FISOUT *O) const;
};

class DEFUZ_MaxCrisp : public DEFUZ
//**********************************
{
  protected :

    //! The classes are not used by the inference process but are needed for
    //! learning purposes.
    double * Classes;

  DEFUZ_MaxCrisp& operator=(const DEFUZ_MaxCrisp&);
  DEFUZ_MaxCrisp(const DEFUZ_MaxCrisp&);
  public :
  
    DEFUZ_MaxCrisp()
    : DEFUZ() 
    { 
      Classes = NULL;
      Thres = AMBIGU_THRES;
    }

  virtual ~DEFUZ_MaxCrisp()
    {
      if(Classes) delete [] Classes;
    }

  void SetClasses(double *c, int n)
    {
      if(Classes) delete [] Classes;
      Classes = NULL;

      NbClas = n;
      if(NbClas == 0) return;

      Classes = new double [NbClas];
      for(int i = 0; i < NbClas; i++) Classes[i] = c[i];
    }

  //! 'T' is the array of output observed values in data file, 'n' its size
  void InitClasses(double *T, int n)
    {
      if(Classes) delete [] Classes;
      Classes = NULL;
      
      InitUniq(T, n, Classes, NbClas);
    }      

  double * ClasLabel(void) { return Classes; } 

  //! The MaxCrisp defuzzification operator is used to provide a more precise performance 
  //! index when dealing with crisp outputs, whose label represents non ordered classes. 
  //! Diagnotic applications are a good example. It is well suited for learning processes. 
  //! The inferred output is the crisp label corresponding to the maximum degree of match. 
  //! When the difference between two levels is less than a threshold, an AMBIGUITY alarm is fired. 
  //! In this case the inferred output is the MF corresponding to the maw and, when the max 
  //! are equal the MF with the smaller number.  
  //! The error associated to an item, by the Perf function, depends on the value of the  
  //! classification flag. When it is set to 'no', the error is (1 - max of matching degree) when 
  //! the rule conclusion corresponding to max is the observed output, 1 otherwise. This error 
  //! can be interpreted as a distance between the inferred value and the observed one. 
  //! When the classification flag is set to 'yes', the error is either 0 or 1. 
  //! The value of 1 used in the case of misclassification, 
  //! whatever the inferred value, it aims to consider that there 
  //! is not any difference, or preference, among the other wrong classes.
  double EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display);

  void WriteHeader(FILE *p, FISOUT *O) const;
};

class DEFUZ_ImpFuzzy : public DEFUZ
//*********************************
// Output is calculated by taking the center of the kernel when exists, 
// otherwise the center of the support 

{
  protected :

  public :
 
 DEFUZ_ImpFuzzy() : DEFUZ() {}
  
  virtual ~DEFUZ_ImpFuzzy() {}
  
  virtual double EvalOut(RULE ** TabR, int NbR, FISOUT * O, FILE * fa, FILE *display);

  virtual void WriteHeader(FILE *p, FISOUT *O) const;

};

//************************   DEFUZ.H   ************************



