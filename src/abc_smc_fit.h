//
// Created by balint.galgoczi on 2021.11.20..
//
#ifndef ABC_SMC_ALGORITHM_ABC_SMC_FIT_H
#define ABC_SMC_ALGORITHM_ABC_SMC_FIT_H

#include "../include/third-party-library/pcg-cpp/pcg_random.hpp"
#include "models/peak_model.h"
#include "models/simulator.h"
#include "job_data.cpp"

class AbcSmcFit {
private:
    int draws, nparams, n_steps;
    double beta, epsilon;
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> posteriors;
    Eigen::Array<double, Eigen::Dynamic, 1> prior_likelihoods;
    Eigen::Array<double, Eigen::Dynamic, 1> likelihoods;
    Eigen::Array<double, Eigen::Dynamic, 1> scalings;
    Eigen::Array<double, Eigen::Dynamic, 1> acc_per_chain;
    Eigen::Array<double, Eigen::Dynamic, 1> tempered_logp;
    Eigen::Array<double, Eigen::Dynamic, 1> importance_weights;
    std::unique_ptr<Simulator> model;

    void initializeArrays(double);
    void calculateInitialLikelihoods();
    void runMcMcChains(Eigen::LLT<Eigen::MatrixXd> &, pcg32 &);
    void resample(pcg32 &);
    double determineImportanceWeights(double);    
    void tuning(bool, int &, double &, double);
public:
    void Fit(JobData job);    
    AbcSmcFit(std::unique_ptr<Simulator> simulator) {
        this->model = move(simulator);
    };
};

#endif //ABC_SMC_ALGORITHM_ABC_SMC_FIT_H
