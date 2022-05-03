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
    double test;
    std::unique_ptr<Simulator> model;
    void runMcMcChains(int, int, double, double, int,
                   Eigen::ArrayXX<double> &, Eigen::ArrayX<double> &,
                   Eigen::ArrayX<double> &, Eigen::ArrayX<double> &,
                   Eigen::ArrayX<double> &, Eigen::LLT<Eigen::MatrixXd> &,
                   Eigen::ArrayX<double> &, pcg32 &);
public:
    void Fit(JobData job);    
    AbcSmcFit(std::unique_ptr<Simulator> simulator) {
        this->model = move(simulator);
    };
};

#endif //ABC_SMC_ALGORITHM_ABC_SMC_FIT_H
