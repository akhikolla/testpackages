#ifndef RESULT_H
#define RESULT_H

#include <Rcpp.h>
#include <string>
#include <map>
#include <vector>
//Data-class containing all information regarding one simon's two-stage design.
//Also all information regarding (non)-stochastic curtailment can be stored.
class Result
{
public:
    // More than "ri" responses under "n" enrolled patients are needed to proceed with the study
    struct StoppingRule  
    {
        int ri; 
        int n; 
        double pet;
    };

    struct Curtailment
    {
        float cut;  // Chosen cut.
        float en_sc; // Expected sample size (with choosen cut).
        float pet_sc; // Probability of early termination (with choosen cut).
        float type1_errorRate;
        float type2_errorRate;
        // Variables for the 95%-confidence intervall
        float en_lower;
        float en_upper;
        float pet_lower;
        float pet_upper;
        float alpha_lower;
        float alpha_upper;
        float beta_lower;
        float beta_upper;
        std::vector<float*> *stoppingRulesNSC;
        
        Curtailment(const Curtailment& other){
          
          
          cut = other.cut;
          en_sc = other.en_sc;
          pet_sc = other.pet_sc;
          type1_errorRate = other.type1_errorRate;
          type2_errorRate = other.type2_errorRate;
          en_lower = other.en_lower;
          en_upper = other.en_upper;
          pet_lower = other.pet_lower;
          pet_upper = other.pet_upper;
          alpha_lower = other.alpha_lower;
          alpha_upper = other.alpha_upper;
          beta_lower = other.beta_lower;
          beta_upper = other.beta_upper;
          
          stoppingRulesNSC = new std::vector<float*>();
          for(unsigned int i = 0; i < other.stoppingRulesNSC->size(); i++){
            float *entry = new float[3];
            memcpy(entry, other.stoppingRulesNSC->at(i),3*sizeof(float));
            stoppingRulesNSC->push_back(entry);
          }
        }
        
        Curtailment& operator= (const Curtailment& other){
          if(this != &other){
            
            cut = other.cut;
            en_sc = other.en_sc;
            pet_sc = other.pet_sc;
            type1_errorRate = other.type1_errorRate;
            type2_errorRate = other.type2_errorRate;
            en_lower = other.en_lower;
            en_upper = other.en_upper;
            pet_lower = other.pet_lower;
            pet_upper = other.pet_upper;
            alpha_lower = other.alpha_lower;
            alpha_upper = other.alpha_upper;
            beta_lower = other.beta_lower;
            beta_upper = other.beta_upper;
            
            stoppingRulesNSC = new std::vector<float*>();
            for(unsigned int i = 0; i < other.stoppingRulesNSC->size(); i++){
              float *entry = new float[3];
              memcpy(entry, other.stoppingRulesNSC->at(i),3*sizeof(float));
              stoppingRulesNSC->push_back(entry);
            }
          }
          return *this;
        }
        
        Curtailment(){
          stoppingRulesNSC = NULL;
        }
        
        ~Curtailment(){
          if(stoppingRulesNSC != NULL)
            delete stoppingRulesNSC;
        }
    };

    // Default constructor
    Result();
    // Constructor
    // n: Number of patients enrolled in the whole trial.
    // r: Critical value for the whole trial (more than "r" responses needed at the end of the study to reject the null hypothesis).    
    // n1: Bumber of patients enrolled in the first stage.    
    // r1: Critical value for the first stage (more than "r1" responses needed to proceed to the second stage).    
    // alpha: Type I error rate.
    // beta: Type II error rate.
    // petP1: Probability of early termination under the alternativ hypothesis.
    // petP0: Probability of early termination under the null hypothesis.
    // enP1: Expected samplesize under the alternative hypothesis.
    // enP0: Expected samplesize under the null hypothesis.
    // iD: ID of the design.
    // p0: Response probability under the null hypothesis.
    // p1: Response probability under the alternativ hypothesis.
    Result(int n, int r, int n1, int r1, double alpha, double beta, double petP1, double petP0, double enP1, double enP0, int iD, double p0, double p1);
    ~Result();
    // Sets the start and stop value for wich the weight q minimizes the risk p = q * n + (1-q)* enP0.
    void setAdmissible(double start, double stop);
    // Sets the start and stop value for wich the weight q minimizes the risk p = q * n + (1-q)* enP0.
    // Also sets a name for the admissible design (such as "minimax", "optimal" or "admissible" should be used).
    void setAdmissible(double start, double stop, std::string name);
    // Adds a curtailment result to the map "curtailmentResults".
    void addCurtailmentResult(Curtailment curResult);    
    // Returns an R-representation of all informations stored in this class.
    SEXP get_R_Representation();

    void setUseCurtailment(bool useCurtailment);
    void setCut(int cut);
    bool getAdmissible();
    int getN();
    int getR();
    int getN1();
    int getR1();
    double getEnP0();
    double getEnP1();
    double getPetP0();
    std::map<int, Curtailment>* getCurtailmentResults();

private:    
    int iD;
    int n, r, n1, r1;
    double alpha, beta, petP1, petP0, enP1, enP0, admissibleStart, admissibleStop, p0, p1;
    bool admissible;
    std::string *name;
    std::vector<StoppingRule> *stoppingRulesNSC;
    std::map<int, Curtailment> *curtailmentResults;

    // Should (non-)stochastic curtailment be used?
    bool useCurtailment;
    // Specifies the threshold in percent if (non-)stochastic curtailment is used.
    int cut;
};

#endif // RESULT_H
