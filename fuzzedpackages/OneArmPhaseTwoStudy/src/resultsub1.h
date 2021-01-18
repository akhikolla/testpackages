#ifndef RESULTSUB1_H
#define RESULTSUB1_H


#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>

// Data-class containing all information regarding one subset design.
// Also all information regarding (non)-stochastic curtailment can be stored.
class ResultSub1
{
public:

    struct StoppingRule_key
    {
    public:
        StoppingRule_key() : enrolled(0), ep2(0) {}
        StoppingRule_key(int enrolled, int ep2) : enrolled(enrolled), ep2(ep2) {}
        int enrolled, ep2;

        friend bool operator <(const ResultSub1::StoppingRule_key &lhs, const ResultSub1::StoppingRule_key &rhs)
        {
            if(lhs.enrolled != rhs.enrolled){
                return lhs.enrolled < rhs.enrolled;
            }
            else
            {
                return lhs.ep2 < rhs.ep2;
            }
        }
    };
    // More than "t_int" responses for the subset endpoint 
    // or more than "u_int" responses for the superset endpoint
    // under "enrolled_int" enrolled patients are needed to proceed with the study
    struct StoppingRule_SubD1 
    {
    public:
        StoppingRule_SubD1() : t_int(0), u_int(0), enrolled_int(0), cp(0) {}
        StoppingRule_SubD1(int t, int u, int enrolled, long double c_p) : t_int(t), u_int(u), enrolled_int(enrolled), cp(c_p) {}
        int t_int, u_int, enrolled_int;
        long double cp; //Conditional power.
    };

    struct Curtailment_SubD1
    {
    public:
        float cut; // Chosen cut.
        float en_sc; // Expected sample size (with choosen cut).
        float pet_sc; // Probability of early termination (with choosen cut).
        float type1_errorRate;
        float type2_errorRate;

        std::vector<StoppingRule_SubD1> *stoppingRulesNSC;
        
        Curtailment_SubD1(const Curtailment_SubD1& other){
          
          cut = other.cut;
          en_sc = other.en_sc;
          pet_sc = other.pet_sc;
          type1_errorRate = other.type1_errorRate;
          type2_errorRate = other.type2_errorRate;
          
          stoppingRulesNSC = new std::vector<StoppingRule_SubD1>();
          *stoppingRulesNSC = *(other.stoppingRulesNSC);
        }
        
        Curtailment_SubD1& operator= (const Curtailment_SubD1& other){
          if(this != &other){
            
            cut = other.cut;
            en_sc = other.en_sc;
            pet_sc = other.pet_sc;
            type1_errorRate = other.type1_errorRate;
            type2_errorRate = other.type2_errorRate;
            
            stoppingRulesNSC = new std::vector<StoppingRule_SubD1>();
            *stoppingRulesNSC = *(other.stoppingRulesNSC);
          }
          return *this;
        }
        
        Curtailment_SubD1(){
          stoppingRulesNSC = NULL;
        }
        
        ~Curtailment_SubD1(){
          if(stoppingRulesNSC != NULL)
            delete stoppingRulesNSC;
        }
    };

    // Default constructor
    ResultSub1();
    // Constructor
    // n: Number of patients enrolled in the whole trial.
    // r / s: Critical values for the whole trial (more than "r" responses for the subset endpoint 
    //        or more than "s" responses for the superset endpoint are needed at the end of the study to reject the null hypothesis).    
    // n1: Number of patients enrolled in the first stage.    
    // r1: Critical value for the first stage (more than "r1" responses needed to proceed to the second stage).    
    // alpha: Type I error rate.
    // beta: Type II error rate.
    // petP0: Probability of early termination under the null hypothesis.
    // enP0: Expected samplesize under the null hypothesis.
    // iD: ID of the design.
    // pc0: The response probability under the null hypothesis for the subset endpoint.
    // pt0: The response probability under the null hypothesis for the superset endpoint.
    // pc1: The response probability under the alternative hypothesis for the subset endpoint.
    // pt1: The response probability under the alternative hypothesis for the superset endpoint.
    ResultSub1(int n, int r, int s, int n1, int r1, double alpha, double beta, double petP0, double enP0, int iD, double pc0, double pt0, double pc1, double pt1);
    ~ResultSub1();
    
    // Returns true if the design is admissible which is the case if "setAdmissible" was called.
    bool getAdmissible();
    // Sets the start and stop value for wich the weight q minimizes the rist p = q * n + (1-q)* enP0.
    void setAdmissible(double start, double stop);
    // Sets the start and stop value for wich the weight q minimizes the rist p = q * n + (1-q)* enP0.
    // Also sets a name for the admissible design (such as "minimax", "optimal" or "admissible" should be used).
    void setAdmissible(double start, double stop, std::string name);
    int getN();
    void setN(int n);
    int getR();
    void setR(int r);
    int getS();
    void setS(int s);
    int getN1();
    void setN1(int n1);
    int getR1();
    void setR1(int r1);
    double getEnP0();
    void setEnP0(double enP0);
    double getPetP0();
    void setPetP0(double petP0);
    // Adds a curtailment result to the map "curtailmentResults"
    void addCurtailmentResult(Curtailment_SubD1 curResult);
    std::map<int, Curtailment_SubD1>* getCurtailmentResults();
    void setUseCurtailment(bool useCurtailment);
    void setCut(int cut);

    int getID();
    
    // Returns an R-representation of all informations stored in this class.
    SEXP get_R_Representation();
  
private:
    int iD;
    int n, r, n1, r1, s;
    double alpha, beta, petP0, enP0, admissibleStart, admissibleStop;
    double pc0, pt0, pc1, pt1;
    bool admissible;
    // Should (non-)stochastic curtailment be used?
    bool useCurtailment;
    // Specifies the threshold in percent if (non-)stochastic curtailment is used.
    int cut;
    std::string *name;

    std::map<int, Curtailment_SubD1> *curtailmentResults;
};


#endif // RESULTSUB1_H
