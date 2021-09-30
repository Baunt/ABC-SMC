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

std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &m)
{
    std::vector<std::vector<double>> result(m[0].size(), std::vector<double>(m.size()));

    for (std::vector<double>::size_type i(0); i < m[0].size(); ++i)
        for (std::vector<double>::size_type j(0); j < m.size(); ++j)
            result[i][j] = m[j][i];

    return result;
}

double arithmetic_mean(const std::vector<double> &vector){
    double sum = 0;
    for (int i = 0; i < vector.capacity(); i++)
    {
        sum += vector[i];
    }

    return sum / vector.capacity();
}