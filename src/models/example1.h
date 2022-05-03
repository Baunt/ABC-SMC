#include "simulator.h"
#include "../probability_distribution.h"

class Example1 : public Simulator {
public:
    std::vector<NormalDistribution> InitialGuess;
    Eigen::ArrayXX<double> GenerateInitialPopulation(int nsamples, int nparams, pcg32 & rng);
    double ErrorCalculation(Eigen::ArrayX<double> diffSpectrum);
    double PriorLikelihood(Eigen::ArrayX<double> posterior);
    int NumberOfParameters();
};
