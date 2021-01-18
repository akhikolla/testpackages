#include "sub1.h"
#include "sub1Module.h"
using namespace Rcpp;

Sub1Design::Sub1Design(): Design()
{
  multinomialLookupTable = new std::map<MultiKey, long double>();
  alphaLookupTable = new std::map<AlphaBetaKey, long double>();
  betaLookupTable = new std::map<AlphaBetaKey, long double>();
  simon = new SimonDesign();
  allSub1Results = new std::vector<ResultSub1*>();
  n = 0;
  n1 = 0;
  r = 0;
  r1 = 0;
  s = 0;
  maxn = 0;
  alpha = 0;
  beta = 0;
  p0 = 0;
  p1 = 0;
  pc0 = 0;
  pt0 = 0;
  pc1 = 0;
  pt1 = 0;
  enCurrent = 100000000000;
}

Sub1Design::~Sub1Design()
{
  delete multinomialLookupTable;
  delete alphaLookupTable;
  delete betaLookupTable;
  delete simon;
  delete allSub1Results;
}

void Sub1Design::setAlpha(double a){ Design::setAlpha(a); } 

double Sub1Design::getAlpha(){return Design::getAlpha();}

void Sub1Design::setBeta(double b){ Design::setBeta(b); }

void Sub1Design::setP0(double p0){ Design::setP0(p0); }

void Sub1Design::setP1(double p1){ Design::setP1(p1); }


void Sub1Design::setPc0(double pc0)
{
  alphaLookupTable->clear();
  betaLookupTable->clear();
  this->pc0 = pc0; 
}
void Sub1Design::setPt0(double pt0)
{
  alphaLookupTable->clear();
  betaLookupTable->clear();
  this->pt0 = pt0; 
}
void Sub1Design::setPc1(double pc1)
{
  alphaLookupTable->clear();
  betaLookupTable->clear();
  this->pc1 = pc1; 
}
void Sub1Design::setPt1(double pt1)
{
  alphaLookupTable->clear();
  betaLookupTable->clear();
  this->pt1 = pt1; 
}


long double Sub1Design::logFact(int n){ return Design::logFact(n); }

long double Sub1Design::bin(int n, int r, long double p){ return Design::bin(n, r, p); }

long double Sub1Design::binsum(int n, int r, long double p){ return Design::binsum(n, r, p); }

long double Sub1Design::multinomial(int x1, int x2, int n, long double p1, long double p2)
{
    MultiKey key(x1, x2, n, p1, p2);
    std::map<MultiKey, long double>::iterator it = multinomialLookupTable->find(key);
    if(it != multinomialLookupTable->end())
    {
        return it->second;
    }
    else
    {
        long double prob = 0;
        prob = exp(logFact(n) - (logFact(x1) + logFact(x2) + logFact(n - x1 - x2))) * pow(p1, x1) * pow(p2, x2) * pow(1 - p1 - p2, n - x1 - x2);

        multinomialLookupTable->insert(std::pair<MultiKey,long double>(key,prob));

        return prob;
    }

}

// Based on formula 4.39 from "Kunz C (2011) To-Stage Designs for Phase II Trials with One or Two Endpoints. 
// PH.D. thesis, Medical Faculty of Heidelberg, Rupechts-Karls-Universitaet"
long double Sub1Design::multinomialTest(int r1, int r, int s, int n1, int n, long double p1, long double p2)
{
  if(p1 < p2)
  {
    int x1, x2, y1, y2;
    int n2 = n - n1;
    int y1_sum_boader = s < n1 ? s : n1;
    int x2_sum_boader = 0;
    int y2_sum_boader = 0;
    long double result = 0;
    long double tmp1, tmp2;
    
    for(x1 = 0; x1 <= r1; x1++)
    {
      for(y1 = x1; y1 <= y1_sum_boader; y1++)
      {
        tmp1 = multinomial(x1, y1 - x1, n1, p1, p2 - p1);
        
        x2_sum_boader = (r-x1) < (s-y1) ? (r-x1) : (s-y1);
        y2_sum_boader = n2 < (s-y1) ? n2 : (s-y1);
        tmp2 = 0;
        
        for(x2 = 0; x2 <= x2_sum_boader; x2++)
        {
          for(y2 = x2; y2 <= y2_sum_boader; y2++)
          {
            tmp2 += multinomial(x2, y2- x2, n2, p1, p2 - p1);
          }
        }
        
        result += tmp1 * tmp2;
      }
    }
    return result;
  }
  else
    return 0;
  
  
}

long double Sub1Design::aproximateMaxN()
{
  if(pc0 != 0 && pt0 != 0 && pc1 != 0 && pt1 != 0 && alpha != 0 && beta != 0)
  {
    int minn1 = beta != (1-pc1) ? (ceil( log(beta) / log(1-pc1) )) : 2;
    int minn = minn1 + 1;
    
    simon->setP0(pc0);
    simon->setP1(pc1);
    simon->setAlpha(alpha);
    simon->setBeta(beta);
    int startValue = simon->aproximateMaxN();
  
    return aproximateMaxNInternal(alpha, beta, pc0, pt0, pc1, pt1, minn, startValue);  
  }
  return -1;
}

long double Sub1Design::aproximateMaxNInternal(long double alpha, long double beta, long double pc0, long double pt0, long double pc1, long double pt1, int n1, int startValue)
{
  bool stop = false;
  int maxn = startValue > n1 ? startValue : (n1+1);
  int r, s, n2;
  long double a, b;
  
  while(!stop)
  {
    n2 = maxn - n1;
    r = 1;
    while(r <= n2 && !stop)
    {
      s = r;
      while((s <= (maxn -1)) && !stop)
      {
        a = this->calcAlpha(n1, 0, maxn, r, s, pc0, pt0);
        if(a < alpha)
        {
          b = this->calcBeta(n1, 0, maxn, r, s, pc1, pt1);
          if(b < beta)
          {
            stop = true;
          }
        }
        if(!stop)
          s++;  
      }
      if(!stop)
        r++;
    }
    if(!stop)
      maxn++;
  }
  this->maxn = maxn;
  
  return maxn;
}

// Based on formula 4.11 from "Kunz C (2011) To-Stage Designs for Phase II Trials with One or Two Endpoints. 
// PH.D. thesis, Medical Faculty of Heidelberg, Rupechts-Karls-Universitaet"
long double Sub1Design::calcAlpha(int n1, int r1, int n, int r, int s, long double pc0, long double pt0)
{
  if(pc0 < pt0)
  {
    long double alpha = 0;
    long double tmp1 = 0;
    long double tmp2 = 0;
    int sumboarder1 = r < n1 ? r:n1;
    int sumboarder2 = s < n1 ? s:n1;
    int sumboarder3 = 0;
    int sumboarder4 = 0;
    int n2 = n - n1;
    int x1, y1, x2, y2;
    std::map<AlphaBetaKey, long double>::iterator it;
    AlphaBetaKey key;
    alpha = binsum(n1, r1, pc0);
    
    for(x1 = r1 + 1; x1 <= sumboarder1; x1++)
    {
      for(y1 = x1; y1 <= sumboarder2; y1++)
      {
        sumboarder3 = (r - x1) < (s - y1) ? (r-x1):(s-y1);
        sumboarder4 = n2 < (s-y1) ? n2:(s-y1);
        tmp2 = 0;
        key.setValues(sumboarder3, sumboarder4, n2);
        it = alphaLookupTable->find(key);
        if(it != alphaLookupTable->end())
        {
          tmp2 = it->second;
        }
        else
        {
          for(x2 = 0; x2 <= sumboarder3; x2++)
          {
            for(y2 = x2; y2 <= sumboarder4; y2++)
            {
              tmp2 += multinomial(x2, (y2 -x2), n2, pc0, pt0 - pc0);
            }
          }
          alphaLookupTable->insert(std::pair<AlphaBetaKey,long double>(key,tmp2));
        }
        
        tmp1 += multinomial(x1, (y1-x1), n1, pc0, pt0 - pc0) * tmp2;
      }
    }
    
    alpha += tmp1;
    alpha = 1 - alpha;
    
    return alpha;
  }
    
  return 0;
}

// Based on formula 4.12 from "Kunz C (2011) To-Stage Designs for Phase II Trials with One or Two Endpoints. 
// PH.D. thesis, Medical Faculty of Heidelberg, Rupechts-Karls-Universitaet"
long double Sub1Design::calcBeta(int n1, int r1, int n, int r, int s, long double pc1, long double pt1)
{
  if(pc1 < pt1)
  {
    long double beta = 0;
    long double tmp1 = 0;
    long double tmp2 = 0;
    int sumboarder1 = r < n1 ? r:n1;
    int sumboarder2 = s < n1 ? s:n1;
    int sumboarder3 = 0;
    int sumboarder4 = 0;
    int n2 = n - n1;
    int x1, y1, x2, y2;
    std::map<AlphaBetaKey, long double>::iterator it;
    AlphaBetaKey key;
    beta = binsum(n1, r1, pc1);
    
    for(x1 = r1 + 1; x1 <= sumboarder1; x1++)
    {
      for(y1 = x1; y1 <= sumboarder2; y1++)
      {
        sumboarder3 = (r - x1) < (s - y1) ? (r-x1):(s-y1);
        sumboarder4 = n2 < (s-y1) ? n2:(s-y1);
        tmp2 = 0;
        key.setValues(sumboarder3, sumboarder4, n2);
        it = betaLookupTable->find(key);
        
        if(it != betaLookupTable->end())
        {
          tmp2 = it->second;
        }
        else
        {
          for(x2 = 0; x2 <= sumboarder3; x2++)
          {
            for(y2 = x2; y2 <= sumboarder4; y2++)
            {
              tmp2 += multinomial(x2, (y2 -x2), n2, pc1, pt1 - pc1);
            }
          }
          betaLookupTable->insert(std::pair<AlphaBetaKey,long double>(key,tmp2));
        }
        tmp1 += multinomial(x1, (y1-x1), n1, pc1, pt1 - pc1) * tmp2;
      }
    }
    
    beta += tmp1;
    
    return beta;
  }
  
  return 0;
}

long double Sub1Design::calcPet_p0(int r1, int n1, long double pc0)
{
  return binsum(n1, r1, pc0);
}

long double Sub1Design::calcEN_p0(int n1, int n, long double pet_p0)
{
  long double ld_n1 = (long double)n1;
  long double ld_n2 = (long double)n - ld_n1;
  return (ld_n1 + (1-pet_p0)*ld_n2);
}

bool Sub1Design::getDesign(int r1, int n1, int r, int s, int n, bool skipN1)
{
    if(multinomialTest(r1, r, s, n1, n, pc1, pt1) < beta)
    {
      long double a,b;
      a = this->calcAlpha(n1, r1, n, r, s, pc0, pt0);
      if(a < alpha)
      {
        b = this->calcBeta(n1, r1, n, r, s, pc1, pt1);

        if(b < beta) // design found
        {
          if(firstDesignFound)  //don't save the first found design because calculateStudySolutions skips r1 untill the first design was found
          {
              long double pet_p0 = this->calcPet_p0(r1, n1, pc0);
              long double en_p0 = this->calcEN_p0(n1, n, pet_p0);
              
              if(en_p0 <= enCurrent)
              {
                //save found design
                allSub1Results->push_back(new ResultSub1(n, r, s, n1, r1, a, b, pet_p0, en_p0, allSub1Results->size(), pc0, pt0, pc1, pt1));

                //save position of minimax design
                if(!minMaxFound)
                {
                  minMaxFound = true;
                  miniN = n;
                }

                if(n == miniN)
                {
                  miniMaxPos++;
                }

                optimalPos++;
                enCurrent = en_p0;

                return true;
              }
              else if(skipN1)
              {
                  allSub1Results->push_back(new ResultSub1(n, r, s, n1, r1, a, b, pet_p0, en_p0, allSub1Results->size(), pc0, pt0, pc1, pt1));
                  return true;
              }
              else
               return false;
          }
          return true;

        }// if(b < beta...
        else
            return false;
      } // if(a < alpha...
      else
          return false;
    } // if(multinomialTest...
    else
        return false;
}

void Sub1Design::calculateStudySolutions(bool skipS, bool skipR, bool skipN1, int lowerBorder, int upperBorder)
{
  if(pc0 != 0 && pt0 != 0 && pc1 != 0 && pt1 != 0 && alpha != 0 && beta != 0)
  {
    if((maxn == 0) & (upperBorder == 0))
        maxn = this->aproximateMaxN();

    allSub1Results->clear();
    
    miniMaxPos = -1;
    optimalPos = -1;
    
    miniN = maxn;
    minMaxFound = false;
    
    int minn1 = beta != (1-pc1) ? (ceil( log(beta) / log(1-pc1) )) : 2;
    int minn = lowerBorder == 0 ? (minn1 + 1) : lowerBorder;
    int local_maxN = (upperBorder == 0) ? maxn : upperBorder;

    enCurrent = 100000000000;
    int n, n1, r1, r, s;
    int  r1Max = 0;

    // "r1" should not be skipped if "lowerBorder" = "upperBorder" therefore "firstDesignFound" is set to true in this case.
    firstDesignFound = (lowerBorder == upperBorder);

    for(n = minn; n <= local_maxN; n++)
    {
      for(n1 = minn1; n1 <= (n-1) && n1 <= enCurrent; n1++)
      {
        if(!firstDesignFound) // skip "r1" as long as no design was found
            r1Max = 1;
        else
            r1Max = n1;
        for(r1 = 0; r1 < r1Max; r1++)
        {
          if(binsum(n1, r1, pc1) < beta) 
          {
            for(r = n - n1 + r1; r >= (r1 + 1); r--)
            {
              if(simon->calcAlpha(n1, r1, n, r, pc0) < alpha)
              {
                for(s = r; s < n; s++)
                {
                  if(getDesign(r1, n1, r, s, n, skipN1))
                  {
                      if(!firstDesignFound)
                      {
                          firstDesignFound  = true;
                          n = n -2;
                          n1 = (n);
                          r1 = n1;
                          r = 0;
                          s = n;
                      }
                      if(skipS)
                          s = n;
                      if(skipR)
                          r = r1;
                      if(skipN1)
                          n1 = n;
                  }

               } // for(int s = r; ...
              } // if(simon->calcAlpha ...
              else
              {
                r = r1 -1;
              }
            } // for(int r = n - n1 + r1; ...
          } // if(binsum(n1, r1, pc1) ...
          else
          {
            r1 = n1;
          }
        } // for(int r1 = 0; ...
      } // for(int n1 = minn1; ...
    } // for(int n = minn; ...
    if(!skipN1)
      this->setAdmissible(allSub1Results);
  } // if(pc0 != 0 && pt0 != 0 ...
}

long double Sub1Design::get_p_exact(int t, int u, int r1, int n1, int n, long double pc0, long double pt0)
{
    long double sum = 0;
    int n2 = n - n1;
    int x2_start = 0;
    int y2_start = 0;
    long double tmp1, tmp2, tmp3, tmp4;

    int x1, x2, y1, y2;

    for(x1 = r1 + 1; x1 <= n1; x1++)
    {
        x2_start = 0 > (t - x1) ? 0 : (t - x1);
        if(x2_start <= n2)
        {
            for(y1 = x1; y1 <= n1; y1++)
            {
                for(x2 = x2_start; x2 <= n2; x2++)
                {
                    y2_start = x2 > (u - y1) ? x2 : (u - y1);
                    for(y2 = y2_start; y2 <= n2; y2++)
                    {
                        tmp1 = exp(logFact(n1) - logFact(n1-y1) -logFact(y1)) * exp(logFact(y1) - logFact(y1-x1) -logFact(x1))
                                * exp(logFact(n2) - logFact(n2-y2) -logFact(y2)) * exp(logFact(y2) - logFact(y2-x2) -logFact(x2));
                        tmp2 = pow(pc0, (x1+x2));
                        tmp3 = pow((pt0-pc0), (y1 + y2 - x1 - x2));
                        tmp4 = pow((1-pt0), (n - y1 - y2));


                        sum = sum + tmp1*tmp2 * tmp3 * tmp4;                        
                    }
                }
            }
        }
    }

    return sum;

}

// Based on formula 4.37 from "Kunz C (2011) To-Stage Designs for Phase II Trials with One or Two Endpoints. 
// PH.D. thesis, Medical Faculty of Heidelberg, Rupechts-Karls-Universitaet"
long double Sub1Design::get_conditionalPower(int t, int u, int enrolled, int r1, int n1, int r, int s, int n, long double pc1, long double pt1)
{
    long double cp;
    long double sum = 0;

    int r1_t = r1 - t;
    int r_t = r - t;
    int s_u = s - u;
    int n1_enrolled = n1 - enrolled;
    int n_enrolled = n - enrolled;


    if( ((0 <= n1_enrolled) && (n1_enrolled <= r1_t)) || (((n1 < enrolled) && ( enrolled < n)) && ((0 <= n_enrolled) && (n_enrolled <= r_t))
                                                            && ((0 <= n_enrolled) && (n_enrolled <= s_u)))) // unabel to reach study goal
    {
        cp = 0;
        return cp;
    }
    else if( ((0 <= r1_t) && (r1_t < n1_enrolled)) && (u <= s) )  //in stage one and nothing sure
    {
        int patLeft = n1 - enrolled;
        int n2 = n - n1;
        int x1, y1, x2, y2;

        for(x1 = 0; x1 <= r1_t; x1++) // propperbility to procceed to stage two
        {
            sum += exp(logFact(patLeft) - logFact(patLeft-x1) -logFact(x1)) * pow(pc1, x1) * pow((1- pc1), (patLeft - x1));
        }

        int sum1_boader = r_t < (s_u < patLeft ? s_u : patLeft) ? r_t : (s_u < patLeft ? s_u : patLeft) ;   // min(r_t, s_u, patLeft)
        int sum2_boader = s_u < patLeft ? s_u : patLeft;
        int sum3_boader = 0;
        int sum4_boader = 0;

        for(x1 = r1_t +1; x1 <= sum1_boader; x1++)
        {
            for(y1 = x1; y1 <= sum2_boader; y1++)
            {
                sum3_boader = (r_t - x1) < (s_u - y1) ? (r_t - x1) : (s_u - y1);
                sum4_boader = n2 < (s_u - y1) ? n2 : (s_u - y1);

                for(x2 = 0; x2 <= sum3_boader; x2++)
                {
                    for(y2 = x2; y2 <= sum4_boader; y2 ++)
                    {
                        sum += binomKoeff(patLeft, y1) * binomKoeff(y1, x1) * binomKoeff(n2, y2) * binomKoeff(y2, x2) * pow(pc1,(x1 + x2)) *
                                pow((pt1 - pc1), (y1 + y2 - x1 - x2)) * pow((1 - pt1), (n - enrolled - y1 - y2));
                    }
                }
            }
        }

        cp = 1 - sum;
        return cp;
    }
    else if( ((0 <= r1_t) && (r1_t < n1_enrolled)) && (u > s)) // if study procceeds to second stage  H0 can be rejected
    {
        int patLeft = n1 - enrolled;
        for(int x = 0; x <= r1_t; x++)
        {
            sum += exp(logFact(patLeft) - logFact(patLeft-x) -logFact(x)) * pow(pc1, x) * pow((1- pc1), (patLeft - x));
        }
        cp = 1 - sum;
        return cp;
    }
    else if( ((enrolled <= n1) && ((r1< t) &&(t <= r)) && (u <= s)) || ((enrolled >= n1) && (t<=r) && ((0 <= s_u) &&(s_u <= n_enrolled)))) // study in second stage or it is sure that the study will procceed to the second stage
    {
        int patLeft = n - enrolled;
        int sum1_boader = r_t < s_u ? r_t : s_u;
        for(int x = 0; x <= sum1_boader; x ++)
        {
            for(int y = x; y <= s_u; y ++)
            {
                sum += exp(logFact(patLeft) - logFact(patLeft-y) -logFact(y)) * exp(logFact(y) - logFact(y - x) -logFact(x)) *
                        pow(pc1, x) * pow((pt1 - pc1), (y-x)) * pow((1 - pt1),(patLeft - y));
            }
        }
        cp = 1 - sum;
        return cp;
    }
    else if( (enrolled >= n1) && ((0 <= r_t) && (r_t < n_enrolled) && (n_enrolled <= s_u))) // in second stage and only way to reject H0 are more responses in enpoint one (pc1)
    {
        int patLeft = n - enrolled;

        for(int x = 0; x <= patLeft; x++)
        {
            sum += exp(logFact(patLeft) - logFact(patLeft-x) -logFact(x)) * pow(pc1,x) * pow((1- pc1), (patLeft - x));
        }
        cp = 1 - sum;
        return cp;
    }
    else if( ((enrolled < n1) && ((t > r) || ((r1 < t) && (t <= r) && (u > s)))) || (((n1 <= enrolled) && (enrolled <= n)) && ((t > r) || (u > s))) ) //H0 can be rejected
    {
        cp = 1;
        return cp;
    }

    cp = 0;
    return cp;
}

void Sub1Design::calculateSC(int resID, double cut, int reps, bool all)
{
    if(allSub1Results->size() > 0)
    {
        std::vector<ResultSub1*>::iterator it = allSub1Results->begin();
        ResultSub1 *res;
        bool stop = false;

        while((it != allSub1Results->end()) && !stop )
        {
            res = *it;

            if(res->getID() == resID)
                stop = true;
            it++;
        }

        if(stop) // design with the ID "resID" was found
        {
            ResultSub1::Curtailment_SubD1 cur;
            if(all)
            {
                double cut_int;
                for(cut_int = cut; cut_int <= 1; cut_int = cut_int + 0.05)
                {
                    cur = this->calcSCintern(res->getR1(), res->getR(), res->getS(), res->getN1(), res->getN(), pc1, pt1, pc0, pt0, cut_int, reps);
                    res->addCurtailmentResult(cur);                    
                }
            }
            else
            {
                cur = this->calcSCintern(res->getR1(), res->getR(), res->getS(), res->getN1(), res->getN(), pc1, pt1, pc0, pt0, cut, reps);
                res->addCurtailmentResult(cur);
            }
        }

    }
}

ResultSub1::Curtailment_SubD1 Sub1Design::calcSCintern(int r1, int r, int s, int n1, int n, long double pc1, long double pt1, long double pc0, long double pt0, long double cut, int rep)
{
    ResultSub1::Curtailment_SubD1 cur;
    std::vector<ResultSub1::StoppingRule_SubD1> *stoppingRules = new std::vector<ResultSub1::StoppingRule_SubD1>();
    std::map<ResultSub1::StoppingRule_key, ResultSub1::StoppingRule_SubD1> *stoppingMap = new std::map<ResultSub1::StoppingRule_key, ResultSub1::StoppingRule_SubD1>();
    std::map<ResultSub1::StoppingRule_key, ResultSub1::StoppingRule_SubD1>::iterator it2;
    int enrolled, x, y;
    long double cp = 0;

    ResultSub1::StoppingRule_SubD1 tmp_sr;


    //calculate Stopping rules for the given cut.
    for(x = 0; x <= r; x ++)
    {
        for(y = x; y <= s; y ++)
        {
            for(enrolled = y; enrolled <= n; enrolled++)
            {
                cp = this->get_conditionalPower(x, y, enrolled, r1, n1, r, s, n, pc1, pt1);

                if(cp <= cut)
                {
                    ResultSub1::StoppingRule_SubD1 sr(x, y, enrolled, cp);
                    ResultSub1::StoppingRule_key key(enrolled, y);

                    it2 = stoppingMap->find(key);
                    if(it2 == stoppingMap->end())
                    {
                        stoppingMap->insert(std::pair<ResultSub1::StoppingRule_key, ResultSub1::StoppingRule_SubD1>(key, sr));
                    }
                    else
                    {
                        tmp_sr = it2->second;
                        if(tmp_sr.t_int < sr.t_int)
                        {
                            stoppingMap->erase(it2);
                            stoppingMap->insert(std::pair<ResultSub1::StoppingRule_key, ResultSub1::StoppingRule_SubD1>(key, sr));
                        }
                    }

                    enrolled = n;
                }
            }
        }
    }

    for(it2 = stoppingMap->begin(); it2 != stoppingMap->end(); it2++)
    {
        stoppingRules->push_back(it2->second);
    }

    // Start simulations to estimate the effect of (non-)stochastic curtailment
    float *numOfStopsForRules = new float[stoppingRules->size()];    
    Rcpp::NumericVector randomNumbers(n);
    int replicationCount, i, i1, count_ep1, count_ep2;
    int type1 = 0, type2 = 0;
    float pet_sc = 0;
    bool stop = false;

    std::vector<ResultSub1::StoppingRule_SubD1>::iterator it;
    ResultSub1::StoppingRule_SubD1 rule;
    int tmp = stoppingRules->size();

    for(i = 0; i < tmp; i++)
        numOfStopsForRules[i] = 0;

    Rcpp::RNGScope scope;

    for(replicationCount = 0; replicationCount < rep; replicationCount++)
    {

        randomNumbers = runif(n);

        stop = false;
        for (it = stoppingRules->begin(); (it != stoppingRules->end()) && !stop; it++)
        {
            rule = *it;

            count_ep1 = 0;
            count_ep2 = 0;

            for(i1 =0 ; i1 < rule.enrolled_int; i1++)       
            {
                if (randomNumbers[i1]<=pc1)  
                {
                    count_ep1++;
                    count_ep2++;
                }
                else if(randomNumbers[i1] < pt1)
                {
                    count_ep2++;
                }
            }

            if ((count_ep1 <= rule.t_int) && (count_ep2 <= rule.u_int))
            {
                    type2++;     // Type II error     
                    stop = true;
            }
        }

        stop = false;
        i = 0;
        for(it = stoppingRules->begin(); (it != stoppingRules->end()) && !stop; it++)
        {
            rule = *it;

            count_ep1 = 0;
            count_ep2 = 0;

            for(i1 = 0; i1 < rule.enrolled_int; i1++)
            {
                if(randomNumbers[i1] < pc0)
                {
                    count_ep1++;
                    count_ep2++;
                }
                else if(randomNumbers[i1] < pt0)
                {
                    count_ep2++;
                }
            }

            if( (count_ep1 <= rule.t_int) && (count_ep2 <= rule.u_int))
            {
                type1++;        // Type I error        
                if(rule.enrolled_int <= n1)
                    pet_sc++;               

                numOfStopsForRules[i]++; 
                stop = true; 
            }
            i++;
        }

    }

    //Analyse the simulation results
    float freps = (float)rep;
    cur.cut = cut;
    cur.type1_errorRate = (freps - (float)type1) / freps;
    cur.type2_errorRate = (float)type2 / freps;
    cur.pet_sc = pet_sc / freps;
    cur.en_sc = 0;
    cur.stoppingRulesNSC = stoppingRules;

    i1 = 0;
    for(it = stoppingRules->begin(); it != stoppingRules->end(); it++)
    {
        rule = *it;

        cur.en_sc += rule.enrolled_int * numOfStopsForRules[i1] / freps;
        i1++;
    }

    return cur;

}

int Sub1Design::getSolutionCount()
{
  return allSub1Results->size();
}

int Sub1Design::getMinimaxPos(){return miniMaxPos;}
int Sub1Design::getOptimalPos(){return optimalPos;}


double Sub1Design::calculateIntersection(double slope1, double yIntercept1, double slope2, double yIntercept2)
{
    double slopeDiff, yInterceptDiff;
    slopeDiff = slope1 - slope2;
    yInterceptDiff = yIntercept2 - yIntercept1;
    if( slopeDiff != 0)
        return yInterceptDiff / slopeDiff;
    return -1000;
}

void Sub1Design::setAdmissible(std::vector<ResultSub1 *> *results)
{
  if(results->size() > 0)
  {
    int i = results->size()-1, admissiblePos = results->size()-1, currentAdmissiblePos = 0;
    double tmpX, oldX =0, x=1;
    while(admissiblePos != miniMaxPos)
    {
        for(i=results->size() -1; i >= miniMaxPos; i--)
        {
            if(!(((ResultSub1*)(results->at(i)))->getAdmissible()))
            {
                tmpX = calculateIntersection(results->at(admissiblePos)->getN() - results->at(admissiblePos)->getEnP0(),results->at(admissiblePos)->getEnP0(),
                                             results->at(i)->getN() - results->at(i)->getEnP0(),results->at(i)->getEnP0());
                if(tmpX > oldX && tmpX < x)
                {
                    x = tmpX;
                    currentAdmissiblePos = i;
                }
            }
        }
        //Enter identified admissible design.
        if(admissiblePos == optimalPos)
        {
            results->at(admissiblePos)->setAdmissible(oldX,x,"Optimal");
        }
        else
        {
            results->at(admissiblePos)->setAdmissible(oldX,x,"Admissible ");
        }
        oldX = x;
        x=1;
        admissiblePos = currentAdmissiblePos;
    }

    //Last found admissible design should allways be the minimax design.
    results->at(miniMaxPos)->setAdmissible(oldX,1,"MiniMax");
  }
}


SEXP Sub1Design::getResultsForR()
{
    Rcpp::List df = Rcpp::List::create();
    for(unsigned int i = 0; i < allSub1Results->size(); i++)
    {
        Rcpp::List tmp = Rcpp::List(((ResultSub1*)allSub1Results->at(i))->get_R_Representation());
        df.push_back(tmp);
    }
    
    return df;
}



