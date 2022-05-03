#ifndef SIMULATOR
#define SIMULATOR

#pragma once

#include "simulatorparameters.h"
#include <memory>
#include "../../include/third-party-library/Eigen/Core"
#include "../../include/third-party-library/pcg-cpp/pcg_random.hpp"


class Simulator {
public:
    virtual Eigen::ArrayXX<double> GenerateInitialPopulation(int nsamples, int nparams, pcg32 & rng) = 0;
    virtual double ErrorCalculation(Eigen::ArrayX<double> parameters) = 0;    
    virtual double PriorLikelihood(Eigen::ArrayX<double> posterior) = 0;
    virtual int NumberOfParameters() = 0;
    virtual ~Simulator() = default;

    static std::unique_ptr<Simulator> CreateSimulator(SimulatorParameters & params);    
};

#endif //SIMULATOR