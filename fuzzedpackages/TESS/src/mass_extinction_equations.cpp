#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector equations_pSurvival_rateshift_CPP(NumericVector lambda,
                                                NumericVector mu,
                                                NumericVector rateChangeTimes,
                                                NumericVector massExtinctionSurvivalProbabilities,
                                                double samplingProbability,
                                                NumericVector t_low,
                                                double t_high,
                                                double T,
                                                bool log){

    NumericVector prev_time(clone(t_low));
    NumericVector den(t_low.length(),1.0);
    NumericVector accumulatedRateTime(t_low.length(),0.0);

    for(int j = 0; j < rateChangeTimes.length(); j++){

        double rate = mu[j] - lambda[j];

        for(int i = 0; i < t_low.length(); i++){
            if(t_low[i] < rateChangeTimes[j] && t_high >= rateChangeTimes[j]){
                den[i] += exp(-rate*prev_time[i])*mu[j]/rate*exp(accumulatedRateTime[i])*(exp(rate*rateChangeTimes[j])-exp(rate*prev_time[i]));
                accumulatedRateTime[i] += rate*(rateChangeTimes[j]-prev_time[i])-std::log(massExtinctionSurvivalProbabilities[j]);
                prev_time[i] = rateChangeTimes[j];
                den[i] -= (massExtinctionSurvivalProbabilities[j]-1)*exp(accumulatedRateTime[i]);
            }
        }

    }

    int index = 0;
    int interval = 0;
    for(int i = 0; i < rateChangeTimes.length(); i++){
        if(t_high > rateChangeTimes[i]){
            interval ++;
        }
    }
    if(rateChangeTimes.length() > 0){
        if( interval < lambda.length()-1 ){
            index = interval;
        } else {
            index = lambda.length()-1;
        }
    }
    double rate = mu[index]-lambda[index];

    for(int i = 0; i < prev_time.length(); i++){
        den[i] += exp(-rate*prev_time[i])*exp(accumulatedRateTime[i])*mu[index]/rate*( exp(rate*t_high) - exp(rate*prev_time[i]) );
        accumulatedRateTime[i] += rate*(t_high-prev_time[i])-std::log(samplingProbability);
    }

    for(int i = 0; i < t_low.length(); i++){
        if( t_low[i] < T && t_high >= T ){
            den[i] -= (samplingProbability-1)*exp(accumulatedRateTime[i]);
        }
    }

    NumericVector res(den.length(),0.0);
    for(int i = 0; i < den.length(); i++){
        if(log){
            res[i] = std::log(1.0/den[i]);
        } else {
            res[i] = 1.0/den[i];
        }
    }

    return res;


}

//[[Rcpp::export]]
NumericVector equations_p1_rateshift_CPP(NumericVector lambda,
                                         NumericVector mu,
                                         NumericVector rateChangeTimes,
                                         NumericVector massExtinctionSurvivalProbabilities,
                                         double samplingProbability,
                                         NumericVector t,
                                         double T,
                                         bool log){

    NumericVector a = equations_pSurvival_rateshift_CPP(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,
                                                            samplingProbability,t,T,T,log);

    NumericVector prev_time(clone(t));
    NumericVector rate(t.length(),0.0);

    for(int j = 0; j < rateChangeTimes.length(); j++){
        if( T >= rateChangeTimes[j]){
            for(int i = 0; i < t.length(); i++){
                if( t[i] < rateChangeTimes[j]){
                    rate[i] += (mu[j]-lambda[j])*(rateChangeTimes[j]-prev_time[i])-std::log(massExtinctionSurvivalProbabilities[j]);
                    prev_time[i] = rateChangeTimes[j];
                }
            }
        }
    }

    for(int i = 0; i < rate.length(); i++){
        rate[i] -= std::log(samplingProbability);
        if( T > prev_time[i] ){
            rate[i] += (mu[mu.length()-1] - lambda[lambda.length()-1])*(T - prev_time[i]);
        }
    }

    NumericVector p(t.length(),0.0);
    for(int i = 0; i < t.length(); i++){
        if(log){
            p[i] = 2*a[i]+rate[i];
        } else {
            p[i] = std::log(2*a[i]+rate[i]);
        }
    }

    return p;

}
