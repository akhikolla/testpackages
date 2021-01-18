#include "simon.h"
#include "simonModule.h"
using namespace Rcpp;

SimonDesign::SimonDesign()
{
  miniMaxPos = 0;
  optimalPos = 0;
  maxn = 0;
  allResults = new std::vector<Result*>();
}

SimonDesign::~SimonDesign()
{
    delete allResults;
}

void SimonDesign::setAlpha(double a){ Design::setAlpha(a);}
double SimonDesign::getAlpha(){return alpha;}

void SimonDesign::setBeta(double b){ Design::setBeta(b);}
double SimonDesign::getBeta(){return beta;}

void SimonDesign::setP0(double p0){ Design::setP0(p0);}
double SimonDesign::getP0(){return p0;}

void SimonDesign::setP1(double p1){ Design::setP1(p1);}
double SimonDesign::getP1(){return p1;}

void SimonDesign::setMaxN(int maxn){this->maxn = maxn;}
int SimonDesign::getMaxN(){return maxn;}

int SimonDesign::getMinimaxPos(){return this->miniMaxPos;}

int SimonDesign::getOptimalPos(){return this->optimalPos;}

int SimonDesign::getSolutionCount(){return allResults->size();}

long double SimonDesign::logFact(int n)
{
    return Design::logFact(n);
}

long double SimonDesign::bin(int n, int r, long double p)
{
    return Design::bin(n, r, p);
}

long double SimonDesign::binsum(int n, int r, long double p)
{
    return Design::binsum(n, r, p);
}

long double SimonDesign::calcAlpha(int n1, int r1, int n, int r, long double p0)
{
    //Plausibility checks
    if(r1 > n1 || r1 > r || n1 > n || r >n || p0 > 1)
        return 0;
    long double sum = 0, alpha;
    int x = r1 +1;
    int n2 = n - n1;
    //Calculate alpha
    while((x<=r)  && (x<=n) && ((r-x+1)<=n2) && (x<=n1))
    {
        sum += bin(n1,x,p0)*binsum(n2,r-x,p0);
        x++;
    }
    alpha = 1- (binsum(n1,r1,p0) +sum);

    return alpha;
}

long double SimonDesign::calcBeta(int n1, int r1, int n, int r, long double p1)
{
    //Plausibility checks
    if(r1 > n1 || r1 > r || n1 > n || r >n || p1 > 1)
        return 0;
    long double sum = 0, beta;
    int x = r1 +1;
    int n2 = n - n1;
    //Calculate beta
    while((x<=r)  && (x<=n) && ((r-x+1)<=n2) && (x<=n1))
    {
        sum += bin(n1,x,p1)*binsum(n2,r-x,p1);
        x++;
    }
    beta = binsum(n1,r1,p1) +sum;

    return beta;
}


long double SimonDesign::aproximateMaxN()
{
    long double p1c;
    int minn1;
    //get minimal n1
    //if(beta == (double)(1.0L -p1))
    if(std::abs(beta - (1.0L -p1)) < 0.0001)
    {
        minn1=2;
        p1c = p1 - 0.025;
        //p1c = p1;
    }
    else
    {
        p1c = p1;
        minn1= ceil(log(beta)/log(1-p1c));
    }

    return this->aproximateMaxNInternal(alpha, beta, p0, p1c, minn1);
}

long double SimonDesign::aproximateMaxNInternal(long double alpha, long double beta, long double p0, long double p1, int n1)
{
    int maxn=n1 +1, r;
    bool stop = false;
    long double n2,a,b;
    //estimate maxn
    while(!stop)
    {
        n2 = maxn - n1;
        r = 1;
        while(r<=n2 && !stop)
        {
            a= this->calcAlpha(n1,0,maxn,r,p0);
            b= this->calcBeta(n1,0,maxn,r,p1);
            if(maxn < 50){
            }
            
            if( (a<alpha) && (b<beta))
            {
                maxn += 9;
                stop = true;
            }else
              r++;
        }
        maxn++;
    }
    this->maxn = maxn;
    return maxn;
}

void SimonDesign::calculateStudySolutions()
{
    //clear all
    allResults->clear();
    miniMaxPos = -1;
    optimalPos = -1;
    
    long double p1c;
    double Rp0, PETp0, ENp0, Rp1, PETp1,ENp1;
    int n1,r1,n,r, minn1, minn;

    //get minimal n1
    if(beta == 1 -p1)
    {
        minn1=2;
        p1c = p1 - 0.025;
    }
    else
    {
        p1c = p1;
        minn1= ceil(log(beta)/log(1-p1c));
    }
    minn = minn1 +1;

    if(maxn < minn)
        maxn=minn +1;


    bool minMaxFound = false;
    int miniN=0;
    double currentmin = 100000000;
    n = minn;
    //start search algorithm (following the instruktions given in 
    //Kunz C (2011) To-Stage Designs for Phase II Trials with One or Two Endpoints. 
    //PH.D. thesis, Medical Faculty of Heidelberg, Rupechts-Karls_universitaet.)
    while(n<=maxn)
    {
        for(n1=minn1; n1<(n < currentmin ? n : currentmin); n1++)
                for(r1=0; r1<n1; r1++)
                    if(binsum(n1,r1,p1)<=beta)
                        for(r=r1+1; r<=(n-n1+r1); r++)
                            if(binsum(n,r,p1)<=beta)
                            {
                                Rp0 = this->calcAlpha(n1,r1,n,r,p0);
                                if(Rp0 <= alpha)
                                {
                                    Rp1 = this->calcBeta(n1,r1,n,r,p1);
                                    if(Rp1 <= beta)
                                    {
                                        PETp0 =binsum(n1,r1,p0);
                                        ENp0 =PETp0*n1 +(1-PETp0)*n;
                                        PETp1 = binsum(n1,r1,p1);
                                        ENp1 = PETp1*n1+(1-PETp1)*n;
                                        if(ENp0<=currentmin)
                                        {
                                            allResults->push_back(new Result(n, r, n1, r1, Rp0, Rp1, PETp1, PETp0, ENp1, ENp0, allResults->size(), p0, p1));
                                            currentmin = ENp0;

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
                                        }
                                    }
                                }
                            }
                            else
                                r = n-n1+r1 +1;
        n++;
    }
    //Identify admissible designs
    this->setAdmissible(allResults);
}

// This function is used in "setAdmissible" to identify admissible designs.
double SimonDesign::calculateIntersection(double slope1, double yIntercept1, double slope2, double yIntercept2)
{
    double slopeDiff, yInterceptDiff;
    slopeDiff = slope1 - slope2;
    yInterceptDiff = yIntercept2 - yIntercept1;
    if( slopeDiff != 0)
        return yInterceptDiff / slopeDiff;
    return -1000;
}

void SimonDesign::setAdmissible(std::vector<Result *> *results)
{
  if(results->size() == 1)
  {
    results->at(miniMaxPos)->setAdmissible(0,1,"MiniMax");
  }
  else if(results->size() > 0)
  {
    int i = results->size()-1, admissiblePos = results->size()-1, currentAdmissiblePos = 0;
    double tmpX, oldX =0, x=1;
    while(admissiblePos != miniMaxPos)
    {
        for(i=results->size() -1; i >= miniMaxPos; i--)
        {
            if(!(((Result*)(results->at(i)))->getAdmissible()))
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
    //last found admissible design should allways be the minimax design
    results->at(miniMaxPos)->setAdmissible(oldX,1,"MiniMax");
  }
}

double SimonDesign::getConditionalPower(int rk, int k, int r1, int n1, int r, int n, double p1)
{
    if( ((0<(n1-k)) && ((n1-k)<(r1-rk+1))) || ((r-rk+1)>(n-k)))
        return 0;
    else if( 1<=(r1-rk+1) && (r1-rk+1)<=(n1-k))
    {
        double sum = 0;
        int x,min;
        min = (r-rk) <= (n1-k) ? r-rk : n1-k;
        for(x=r1-rk+1; x<=min;x++)
            sum += bin(n1-k,x,p1) * binsum(n-n1,r-x-rk,p1);
        return 1 - (binsum(n1-k,r1-rk,p1) + sum);
    }
    else if(r1<rk && rk <=r)
    {
        return 1 - binsum(n-k,r-rk,p1);
    }
    else if(rk>r)
        return 1;
    return 0;
}


Result::Curtailment SimonDesign::calcSCIntern(int resID, double cut, int reps)
{
    Result *result = allResults->at(resID);
    int ri, k;
    double cp;
    std::vector<float*> *ri_k_cp = new std::vector<float*>();
    //Identify the stopping rules resulting through the chosen "cut"
    for(int i=0; i <= result->getR(); i++)
    {
        float *entry = new float[3];
        entry[0] = 0;
        entry[1] = 0;
        entry[2] = 0;
        ri_k_cp->push_back(entry);
    }

    for(ri = 0; ri <= result->getR(); ri++)
    {
        for(k = ri; k <= result->getN(); k++)
        {
            if(k == result->getN1() && ri <= result->getR1())
               cp = 0;
           else
                cp = this->getConditionalPower(ri, k, result->getR1(), result->getN1(), result->getR(), result->getN(), p1);

           if(cp <= cut)
           {
               ((float*)ri_k_cp->at(ri))[0] = ri;
               ((float*)ri_k_cp->at(ri))[1] = k;
               ((float*)ri_k_cp->at(ri))[2] = 0;
               k = result->getN();
           }
        }
    }


    // Start simulations to estimate the effect of (non-)stochastic curtailment.
    int type1 = 0, type2 = 0;
    double pet_sc = 0, e1, e2;
    int w, i1, i2, j1, j2,count1, count2;

    Rcpp::RNGScope scope;
    Rcpp::NumericVector randomNumbers(result->getN());

    for(w = 0; w < reps; w++)
    {
        randomNumbers = runif(result->getN());

       for (i1=0; i1<=result->getR(); i1++)
       {
           count2 = 0;
           e1 = ((float*)ri_k_cp->at(i1))[1]; 
           for(j1=0 ; j1< e1; j1++)       
           {
               if (randomNumbers[j1]<=p1)  
                   count2 += 1;
           }
           if (count2<=((float*)ri_k_cp->at(i1))[0])
           {
                   type2 += 1;          
                   i1 = result->getR() + 1;
           }
       }

       for (i2=0; i2<=result->getR(); i2++)
       {
           count1 = 0;
           e2 = ((float*)ri_k_cp->at(i2))[1];
           for(j2=0 ; j2< e2; j2++)
           {
                   if (randomNumbers[j2]<=p0)
                       count1 += 1;
           }
           if (count1 <= ((float*)ri_k_cp->at(i2))[0])
           {
               type1 += 1;
               if (((float*)ri_k_cp->at(i2))[1] <= result->getN1())
                   pet_sc += 1;
               ((float*)ri_k_cp->at(i2))[2] += 1;
               i2 = result->getR() + 1; 
           }
        }
    }
    //Analyse the simulation results.
    float freps = (float)reps;

    Result::Curtailment curResult;
    curResult.cut = cut;
    curResult.type1_errorRate = (freps - (float)type1) / freps;
    curResult.type2_errorRate = (float)type2 / freps;
    curResult.pet_sc = (float)pet_sc / freps;
    curResult.en_sc = 0;
    float sd = 0;
    float temp = 0;
    int i;
    
    for(i = 0; i <= result->getR(); i++)
    {
        ((float*)ri_k_cp->at(i))[2] = ((float*)ri_k_cp->at(i))[2] / freps;
        sd += freps * ((float*)ri_k_cp->at(i))[2] * ((float*)ri_k_cp->at(i))[1] * ((float*)ri_k_cp->at(i))[1];
        curResult.en_sc += ((float*)ri_k_cp->at(i))[2] * ((float*)ri_k_cp->at(i))[1];
    }
    curResult.en_sc += (freps - type1) / freps * result->getN();
    sd = sqrt(((sd + freps * (freps - type1) / freps * result->getN() * result->getN()) - freps *  curResult.en_sc * curResult.en_sc)/((freps - 1) * freps));

    temp = 1.96 * sd;
    curResult.en_lower = curResult.en_sc - temp;
    curResult.en_upper = curResult.en_sc + temp;

    temp = 1.96 * sqrt((curResult.pet_sc) * (1 - curResult.pet_sc) / freps);
    curResult.pet_lower = curResult.pet_sc - temp;
    curResult.pet_upper = curResult.pet_sc + temp;

    temp = 1.96 * sqrt((freps - type1) / freps * (1 - (freps - type1) / freps) / freps);
    curResult.alpha_lower = (freps - type1) / freps - temp;
    curResult.alpha_upper = (freps - type1) / freps + temp;

    temp = sqrt((type2 / freps) * (1 - (type2 / freps)) / freps);
    curResult.beta_lower = type2 / freps - temp;
    curResult.beta_upper = type2 / freps + temp;

    curResult.stoppingRulesNSC = ri_k_cp;
    return curResult;
}

void SimonDesign::calculateSC(int resID, double cut, int reps, bool all)
{
    Result *result = allResults->at(resID);
    if(all)
    {
        int i = 0;
        while(cut < 1)
        {
            if(result->getCurtailmentResults()->find((int)(cut*100 + 0.5)) == result->getCurtailmentResults()->end())
                result->addCurtailmentResult(calcSCIntern(resID,cut,reps));
            cut = cut + 0.05;
            i++;
        }
    }
    else
        result->addCurtailmentResult( calcSCIntern(resID, cut, reps));
}


SEXP SimonDesign::getResultsForR()
{
    Rcpp::List df = Rcpp::List::create();
    for(unsigned int i = 0; i < allResults->size(); i++)
    {
        Rcpp::List tmp = Rcpp::List(((Result*)allResults->at(i))->get_R_Representation());
        df.push_back(tmp);
    }
    
    return df;
}

SEXP SimonDesign::getCurResultForR(int id)
{
    Result *res = (Result*)(allResults->at(id));
    DataFrame tmp;
    if(res->getCurtailmentResults()->size() > 0)
    {
        std::map<int, Result::Curtailment>::iterator it = res->getCurtailmentResults()->begin();
        Result::Curtailment cur = it->second;
        tmp = DataFrame::create( _["Cut"] = NumericVector::create(cur.cut),
                                           _["En_SC"] = NumericVector::create(cur.en_sc),
                                           _["Pet_SC"] = NumericVector::create(cur.pet_sc),
                                           _["Type_1_Errorrate"] = NumericVector::create(cur.type1_errorRate),
                                           _["Type_2_Errorrate"] = NumericVector::create(cur.type2_errorRate));
    }
    return tmp;
}
