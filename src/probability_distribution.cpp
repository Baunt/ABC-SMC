//
// Created by balint.galgoczi on 2021.09.27..
//
#define _USE_MATH_DEFINES

#include "probability_distribution.h"
#include "util.h"

double NormalDistribution::LogP(double value) {
    double piInverse = 1.0 / M_PI;

    return (
            -NormalDistribution::p_tau * (value - NormalDistribution::p_mean)*(value - NormalDistribution::p_mean) +
            log((NormalDistribution::p_tau * piInverse) * 0.5)
            )
            * 0.5;
}

Eigen::ArrayX<double> NormalDistribution::Sample(int draws, pcg32 &rng) {
    return getDistribution(NormalDistribution::p_mean, NormalDistribution::p_sigma, draws, rng);
}