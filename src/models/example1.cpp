#include "example1.h"
#include <cmath>

Eigen::ArrayXX<double> Example1::GenerateInitialPopulation(int nsamples, int nparams, pcg32 & rng) {
    // int nparams = 2;
    Eigen::ArrayXX<double> priors(nsamples, 2);

    NormalDistribution x0 = NormalDistribution(0.3, 0.5);
    NormalDistribution y0 = NormalDistribution(-0.4, 0.5);
    InitialGuess.push_back(x0);
    InitialGuess.push_back(y0);

    for (int i = 0; i < InitialGuess.size(); ++i) {
        Eigen::ArrayX<double> tmp = InitialGuess[i].Sample(nsamples, rng);
        priors.col(i) = tmp;
    }

    return priors;
}

double Example1::ErrorCalculation(Eigen::ArrayX<double> parameters) {
    double x = parameters(0);
    double y = parameters(1);
    double tmp1 = x * sin(20.0 * y) + y * sin(20.0 * x);
    double tmp2 = x * cos(10.0 * y) - y * sin(10.0 * x);
    return cosh(sin(10.0 * x)*x) * tmp1 * tmp1 + cosh(cos(20.0 * y) * y) * tmp2*tmp2 + y*y;
}

double Example1::PriorLikelihood(Eigen::ArrayX<double> posterior) {
    double prior_likelihoods = 0.0;
    for (int i = 0; i < posterior.size(); ++i) {
        prior_likelihoods += InitialGuess[i].LogP(posterior(i));
    }
    return prior_likelihoods;
}

int Example1::NumberOfParameters() {
    return 2;
}


