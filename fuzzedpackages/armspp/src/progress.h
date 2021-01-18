#ifndef SRC_PROGRESS_HPP_
#define SRC_PROGRESS_HPP_

#include <chrono>
#include <cstdint>
#include <sstream>
#include <Rcpp.h>

namespace armspp {

class ProgressBar {
 public:
    explicit ProgressBar(uint64_t nSteps)
        : nSteps_(nSteps),
          currentStep_(0),
          lastCheckStep_(0),
          lastCheckStepTime_(0),
          lastOutputLength_(0) {
        lastCheckTime_ = std::chrono::system_clock::now();
    }

    uint64_t operator+=(uint64_t increment) {
        currentStep_ += increment;
        output();
        return currentStep_;
    }

    uint64_t operator++() {
        currentStep_ += 1;
        output();
        return currentStep_;
    }

 private:
    uint64_t nSteps_;
    uint64_t currentStep_;

    uint64_t lastCheckStep_;
    std::chrono::time_point<std::chrono::system_clock> lastCheckTime_;
    double lastCheckStepTime_;
    unsigned int lastOutputLength_;

    void output() {
        unsigned int percent = 100 * currentStep_ / nSteps_;

        if (currentStep_ == 1 || currentStep_ % 100 == 0) {
            std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsedSeconds = now - lastCheckTime_;
            lastCheckStepTime_ = 1000 * elapsedSeconds.count() / (currentStep_ - lastCheckStep_);

            lastCheckStep_ = currentStep_;
            lastCheckTime_ = now;
        }

        // Write out a string of spaces to overwrite the last output, because the
        // internal R terminal doesn't support ANSI codes
        Rcpp::Rcout << "\r";
        for (unsigned int i = 0; i < lastOutputLength_; ++i) {
            Rcpp::Rcout << " ";
        }

        std::ostringstream sstream;
        sstream << std::fixed << std::setprecision(2)
            << "\r"
            << currentStep_ << "/" << nSteps_ << " (" << percent << "%)"
            << " " << lastCheckStepTime_ << "ms/iteration"
            << " (" << ((nSteps_ - currentStep_) * lastCheckStepTime_) / 1000 << "s remaining)";
        std::string output = sstream.str();
        lastOutputLength_ = output.length();
        Rcpp::Rcout << output;

        if (currentStep_ == nSteps_) {
            Rcpp::Rcout << "\n";
        }

        Rcpp::Rcout.flush();
    }
};

}  // namespace armspp

#endif  // SRC_PROGRESS_HPP_
