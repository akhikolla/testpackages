// prep.h

#ifndef PREP_H_
#define PREP_H_

#include <iostream>
#include <vector>
#include <tuple>

#include <cstddef>
// Note: Rcpp ignores assert()
// #include <cassert>
#include <algorithm>
#include <memory> // unique_ptr

#include <Rcpp.h>

#include "grid.h"
#include "model.h"

namespace prep {
////////////////////////////////////////////////////////////////
// Compute C from clock_phases
////////////////////////////////////////////////////////////////
/*
using Cs = grid::Grid<double, double>;
std::unique_ptr<Cs> compCs_(Rcpp::IntegerVector const& times_of_day,
                            std::vector<double> const& clock_phases) {
  std::size_t const nsamples = times_of_day.size();
  std::unique_ptr<Cs> cs { new Cs(nsamples, clock_phases) };
  auto c0 = cs->begin();
  for (auto ph = clock_phases.begin(); ph != clock_phases.end(); ++ph) {
    auto c = c0.begin();
    for (auto d = times_of_day.begin(); d != times_of_day.end(); ++d)
      *c++ = model::clockC(*d, *ph);
    ++c0;
  }
  return cs;
}

inline std::unique_ptr<Cs> makeCs(Rcpp::IntegerVector const& times_of_day,
                                  Rcpp::NumericVector const& clock_phases) {
  return compCs_(times_of_day, Rcpp::as<std::vector<double>>(clock_phases));
}

inline std::unique_ptr<Cs> makeC(Rcpp::IntegerVector const& times_of_day,
                                 double clock_phases) {
  return compCs_(times_of_day, std::vector<double> { clock_phases });
}
*/
////////////////////////////////////////////////////////////////
// Compute E(e) from env.e.* and gate.e.*
////////////////////////////////////////////////////////////////

// Computes envF(*w, amp, th) for w in [weather_begin, weather_end]
// for all combinations of amp and th.
// Grid: (fs, ampl, th)
std::unique_ptr<grid::Grid<double, double,double>>
compFs_(Rcpp::NumericVector::const_iterator const& weather_begin,
        Rcpp::NumericVector::const_iterator const& weather_end,
        std::vector<double>                 const& amplitude,
        std::vector<double>                 const& threshold) {
  if (weather_begin >= weather_end)
    throw Rcpp::exception("Inconsistent args. (weather_begin >= weather_end)");
  std::size_t const blocksize = weather_end - weather_begin;
  std::unique_ptr<grid::Grid<double, double,double>> fs
  { new grid::Grid<double, double,double>(blocksize, amplitude, threshold) };
  auto f0 = fs->begin();
  for (auto& amp : amplitude)
    for (auto& th : threshold) {
      auto f = f0.begin();
      for (auto w = weather_begin; w != weather_end; ++w)
        *f++ = model::envF(*w, amp, th);
      ++f0;
    }
  return fs;
}

// Computes gateG(t, phase:0, amp, th) for t in [0:periodG-1]
// for all combinations of amp and th.
// Note:
// - pre-integrating gs(t, amp, th) over [t:t+timeStep-1] (for each t)
//   was net negative: complicates code, modest speedup of init, slows down optim
std::unique_ptr<grid::Grid<double, double,double>>
compGs_(std::vector<double> const& amplitude,
        std::vector<double> const& threshold) {
  std::size_t const blocksize = model::periodG;
  std::unique_ptr<grid::Grid<double, double,double>> gs
  { new grid::Grid<double, double,double>(blocksize, amplitude, threshold) };

  auto g0 = gs->begin();
  double const phase = 0;
  for (auto& amp : amplitude)
    for (auto& th : threshold) {
      auto g = g0.begin();
      for (int t = 0; t < model::periodG; ++t)
        *g++ = model::gateG(t, phase, amp, th);
      ++g0;
    }
  return gs;
}

void normalise(std::vector<double>::iterator const& begin,
               std::vector<double>::iterator const& end) {
  // assert(begin < end);
  double sum = 0.0; double min = *begin; double max = *begin;
  for (auto u = begin; u != end; ++u) {
    sum += *u;
    if (*u < min) min = *u;
    if (*u > max) max = *u;
  }
  double a = (min != max) ? 1.0/(max-min) : 1.0;
  double mean = sum/(end-begin);
  for (auto u = begin; u != end; ++u) {
    *u = a * (*u - mean);
  }
}

void rescale(std::vector<double>::iterator const& begin,
             std::vector<double>::iterator const& end,
             double a) {
  // assert(begin < end);
  for (auto u = begin; u != end; ++u) *u *= a;
}

inline int mod(int a, int b) {
  if (b < 0) { a = -a; b = -b; }
  int r = a % b;
  return (r > 0) ? r : b + r;
}

// Optimizaion memo:
// loop size for a temperature grid:
//   4321(period.max+1)*461(samples)*24*2*3 = 286,845,264 points
// All with -O0
// - a call to gateG() takes ~100ns
// - naive loop: 36.34s
//   - dominated by gateG() (4321*461*(24*2*3) * 100/1e9 = 28.68s)
// - with cache:  8.12s
//   - hits: 286,836,624, misses: 8,640 =1440*2*3
//   - calls to gateG() is now irrelevant: 8,640*100/1e9 = 0.001
// - further exploit the fact that t is contiguous: 5.59s
//   - 20ns/tightest loop

// Computes the values of E(e) for all combinations of
// (gate.amp, gate.th, gate.ph, env.amp, env.th, period, samples),
// and stores the result in the Grid data structure defined in grid.h
// - Depends on the chosen weather factor e through fs and gs.
// (val*samples, gate.amp,gate.th,phase,env.amp,env.th,period)
using Es = grid::Grid<double, double,double,int,double,double,int>;
std::unique_ptr<Es> compEs_(bool show_progress,
                            std::vector<int>                  const& times_pickup,
                            Rcpp::IntegerVector               const& times_of_day,
                            grid::Grid<double, double,double> const& fs,
                            grid::Grid<double, double,double> const& gs,
                            std::vector<int>                  const& phase,
                            std::vector<int>                  const& period,
                            int dataStep, int timeStep) {
  if (times_pickup.size() != times_of_day.size())
    throw Rcpp::exception("mismatch of sizes (ngenes) for times_pickup and times_of_day.");

  std::size_t nsamples = times_pickup.size();
  std::vector<double> const& env_amp  = std::get<0>(fs.coords);
  std::vector<double> const& env_th   = std::get<1>(fs.coords);
  std::vector<double> const& gate_amp = std::get<0>(gs.coords);
  std::vector<double> const& gate_th  = std::get<1>(gs.coords);
  if (show_progress) Rcpp::Rcout << "- nsamples(blocksize): " << nsamples << '\n';
  // this ordering of coords comes from performance considerations
  std::unique_ptr<Es> es
  { new Es(nsamples, std::tie(gate_amp, gate_th, phase, env_amp, env_th, period)) };

  auto e0 = es->begin();
  for (auto g0 = gs.cbegin(); g0 != gs.cend(); ++g0) {            // 2*3
    for (auto ph = phase.cbegin(); ph != phase.cend(); ++ph) {    // 24
      for (auto f0 = fs.cbegin(); f0 != fs.cend(); ++f0) {        // 2*1
        for (auto p = period.cbegin(); p != period.cend(); ++p) { // 7
          // fill block
          auto e = e0.begin();
          auto s = times_pickup.cbegin();
          auto d = times_of_day.begin(); // IntegerVector has no cbegin()
          for ( ; e != e0.end(); ++e, ++s, ++d) {   // for 461 (samples)
            // Units:
            // - *s, *d, *ph, *p: in the unit of min
            // - *f: dataStep [f0 + (*s - *p)/dataStep <= f < f0 + *s/dataStep]
            // - *g: is available for all 1440min (not integrated over timeStep)
            // - caller (makeEs_()) is responsible for ensuring timeStep%dataStep==0
            //   so timeStep/dataStep is an integer
            auto f = f0.cbegin() + (*s - *p)/dataStep;
            auto g = g0.cbegin() + mod(*d - *ph - *p, gs.blocksize);
            int fstep = timeStep/dataStep; // pointer advance per timeStep
            // sum over period (innermost tick: in the unit of min)
            double sum = 0.0;
            for (int i = 0; i < *p; i += timeStep, f += fstep) {
              // Here are two options:
              // // Integrate g over [t:t+timeStep-1]
              // // (As noted in the comment for compGs_(), precomputing this
              // // sum over g is actually not a good idea.)
              // for (int j = 0; j < timeStep; ++j, ++g) {
              //   if (g == g0.cend()) g = g0.cbegin();
              //   sum += *f * (*g);
              // }
              // Or can just use the g value at t.
              // (makes init phase very fast, but not so relevant in optim phase)
              if (g >= g0.cend()) g = g0.cbegin() + (g - g0.cend());
              sum += *f * (*g) * timeStep / (double)*p;
              g += timeStep;
            }
            *e = sum;
          }
          // normalise(e0.begin(), e); // normalise for samples
          // rescale(e0.begin(), e, 0.01);
          ++e0;
        }}}}
  return es;
}

// Note:
// - dataStep: time step of the weather_e_data (min), must be a divisor of timeStep
// - timeStep: time step of the model (min)
// - period_e is measured in min, no matter what dataStep/timeStep are.
// - magnitudes of Es should be mostly independent of dataStep/timeStep.
std::unique_ptr<Es> makeEs_(bool show_progress,
                            Rcpp::IntegerVector const& times_pickup,
                            Rcpp::IntegerVector const& times_of_day,
                            Rcpp::NumericVector const& weather_e_data,
                            std::vector<int>    const& period_e,
                            std::vector<double> const& env_e_amplitude,
                            std::vector<double> const& env_e_threshold,
                            std::vector<int>    const& gate_e_phase,
                            std::vector<double> const& gate_e_amplitude,
                            std::vector<double> const& gate_e_threshold,
                            int dataStep, int timeStep) {
  // In the unit of dataStep:
  //
  //  weather_e_data.begin()                                 weather_e_data.end()
  //  v                                                            v
  //  |<---------- w_end_offset -------------------->|             |
  //  |<- w_beg_offset ->|<- period_max ->|          |             |
  //  |                  |<======= to compFs =======>|             |
  //  ^                                   ^          ^             ^
  //  0                             pickup_min  pickup_max   weather_e_data.size()
  int const pickup_min = *std::min_element(times_pickup.begin(), times_pickup.end());
  int const pickup_max = *std::max_element(times_pickup.begin(), times_pickup.end());
  int const period_max = *std::max_element(period_e.cbegin(), period_e.cend());

  // ensure that enough weather data has been passed
  if (pickup_min < period_max || pickup_max < pickup_min ||
      weather_e_data.size()*dataStep < pickup_max)
    throw Rcpp::exception("Inconsistent args. (weather data too short?)");
  // timeStep must be an integral multiple of dataStep
  if (timeStep % dataStep != 0)
    throw Rcpp::exception("timeStep must be an integral multiple of dataStep.");

  if (show_progress) {
    Rcpp::Rcout << "# computing Fs..\n";
    Rcpp::Rcout << "# - weather_e_size: " << weather_e_data.size() << '\n';
    Rcpp::Rcout << "# - dataStep: " << dataStep << '\n';
    Rcpp::Rcout << "# - timeStep: " << timeStep << '\n';
    Rcpp::Rcout << "# - period_max: " << period_max << '\n';
    Rcpp::Rcout << "# - weather_begin: " << pickup_min - period_max
          << " weather_end: " << pickup_max << '\n';
  }
  int const weather_begin_offset = (pickup_min - period_max) / dataStep;
  int const weather_end_offset = pickup_max / dataStep;
  auto fs = compFs_(weather_e_data.begin() + weather_begin_offset,
                    weather_e_data.begin() + weather_end_offset,
                    env_e_amplitude, env_e_threshold);
  if (show_progress) Rcpp::Rcout << "# computing Gs..\n";
  auto gs = compGs_(gate_e_amplitude, gate_e_threshold);

  if (show_progress) Rcpp::Rcout << "# computing Es..\n";
  // convert times_pickup to offsets in fs
  auto times_pickup_ = Rcpp::as<std::vector<int>>(times_pickup);
  for (auto& t : times_pickup_) t = ((int)t/dataStep - weather_begin_offset) * dataStep;
  auto es = compEs_(show_progress,
                    times_pickup_, times_of_day, *fs, *gs, gate_e_phase, period_e,
                    dataStep, timeStep);
  return es;
}

std::unique_ptr<Es> makeEs(Rcpp::IntegerVector const& times_pickup,
                           Rcpp::IntegerVector const& times_of_day,
                           Rcpp::NumericVector const& weather_e_data,
                           Rcpp::IntegerVector const& period_e,
                           Rcpp::NumericVector const& env_e_amplitude,
                           Rcpp::NumericVector const& env_e_threshold,
                           Rcpp::IntegerVector const& gate_e_phase,
                           Rcpp::NumericVector const& gate_e_amplitude,
                           Rcpp::NumericVector const& gate_e_threshold,
                           int dataStep, int timeStep) {
  auto period_e_         = Rcpp::as<std::vector<int>>(period_e);
  auto env_e_amplitude_  = Rcpp::as<std::vector<double>>(env_e_amplitude);
  auto env_e_threshold_  = Rcpp::as<std::vector<double>>(env_e_threshold);
  auto gate_e_phase_     = Rcpp::as<std::vector<int>>(gate_e_phase);
  auto gate_e_amplitude_ = Rcpp::as<std::vector<double>>(gate_e_amplitude);
  auto gate_e_threshold_ = Rcpp::as<std::vector<double>>(gate_e_threshold);

  return makeEs_(false, // show_progress
                 times_pickup, times_of_day, weather_e_data,
                 period_e_, env_e_amplitude_, env_e_threshold_,
                 gate_e_phase_, gate_e_amplitude_, gate_e_threshold_,
                 dataStep, timeStep);
}

std::unique_ptr<Es> makeE(Rcpp::IntegerVector const& times_pickup,
                          Rcpp::IntegerVector const& times_of_day,
                          Rcpp::NumericVector const& weather_e_data,
                          int    period_e,
                          double env_e_amplitude,
                          double env_e_threshold,
                          int    gate_e_phase,
                          double gate_e_amplitude,
                          double gate_e_threshold,
                          int dataStep, int timeStep) {
  std::vector<int>    const period_e_         { period_e };
  std::vector<double> const env_e_amplitude_  { env_e_amplitude };
  std::vector<double> const env_e_threshold_  { env_e_threshold };
  std::vector<int>    const gate_e_phase_     { gate_e_phase };
  std::vector<double> const gate_e_amplitude_ { gate_e_amplitude };
  std::vector<double> const gate_e_threshold_ { gate_e_threshold };

  return makeEs_(false, // show_progress
                 times_pickup, times_of_day, weather_e_data,
                 period_e_, env_e_amplitude_, env_e_threshold_,
                 gate_e_phase_, gate_e_amplitude_, gate_e_threshold_,
                 dataStep, timeStep);
}

////////////////////////////////////////////////////////////////
} // namespace prep
#endif // PREP_H_
