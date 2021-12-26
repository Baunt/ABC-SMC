//
// Created by balint.galgoczi on 2021.09.27..
//

#ifndef ABC_SMC_ALGORITHM_PROBABILITY_DISTRIBUTION_H
#define ABC_SMC_ALGORITHM_PROBABILITY_DISTRIBUTION_H


#include <string>
#include <vector>
#include "../include/third-party-library/Eigen/Core"
#include "../include/third-party-library/pcg-cpp/pcg_random.hpp"

class ProbabilityDistribution {
public:
    double expected_value;
    double uncertainty;

    ProbabilityDistribution(double expectedValue, double uncertainty){
        this->expected_value = expectedValue;
        this->uncertainty = uncertainty;
    }
};

class NormalDistribution : public ProbabilityDistribution{
private:
    double p_tau;
public:
    NormalDistribution(double expectedValue, double uncertainty)
            : ProbabilityDistribution(expectedValue, uncertainty) {
        this->p_tau = 1 / (uncertainty * uncertainty);
    }

    double LogP(double value);
    Eigen::ArrayX<double> Sample(int draws, pcg32 & rng);
};


#endif //ABC_SMC_ALGORITHM_PROBABILITY_DISTRIBUTION_H
