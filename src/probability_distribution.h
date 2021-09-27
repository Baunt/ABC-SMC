//
// Created by balint.galgoczi on 2021.09.27..
//

#ifndef ABC_SMC_ALGORITHM_PROBABILITY_DISTRIBUTION_H
#define ABC_SMC_ALGORITHM_PROBABILITY_DISTRIBUTION_H


#include <string>
#include <vector>

class ProbabilityDistribution {
private:
    std::string p_name;
    std::string p_distributionType;
    double p_expected_value;
    double p_uncertainty;

public:
    ProbabilityDistribution(std::string name, std::string distributionType, double expectedValue, double uncertainty){
        this->p_distributionType = distributionType;
        this->p_name = name;
        this->p_expected_value = expectedValue;
        this->p_uncertainty = uncertainty;
    }
};

class NormalDistribution : public ProbabilityDistribution{
private:
    double p_sigma;
    double p_tau;
    double p_mean;
public:
    NormalDistribution(std::string name, std::string distributionType, double expectedValue, double uncertainty)
            : ProbabilityDistribution(name, distributionType, expectedValue, uncertainty) {
        this->p_mean = expectedValue;
        this->p_sigma =uncertainty;

        this->p_tau = 1 / (p_sigma * p_sigma);
    }

    double LogP(double value);
    std::vector<double> Sample(int draws);
};


#endif //ABC_SMC_ALGORITHM_PROBABILITY_DISTRIBUTION_H
