//
// Created by balint.galgoczi on 2021.08.10..
//

#include "util.h"

std::vector<double> getDistribution(double x_mu, double x_sigma, size_t numberOfValues){
    // shortened compared to your example:
    std::mt19937 gen((std::random_device())());
    // create temporary (anonymous)     ^^
    // instance and call it immediately    ^^
    // afterwards

    std::normal_distribution<double> nd(x_mu, x_sigma);
    std::vector<double> result;
    result.reserve(numberOfValues);
    while(numberOfValues-- > 0)
    {
        // shorter than above: using result of previous
        // function (functor!) call directly as argument to next one
        result.push_back(nd(gen));
    }

    // finally something familiar from python:
    return result;
}