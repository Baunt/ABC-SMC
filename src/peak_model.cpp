//
// Created by balint.galgoczi on 2021.07.28..
//
#define _USE_MATH_DEFINES

#include "peak_model.h"
#include <cmath>

std::vector<double> PeakModel::Gaussian() {

    double sigma = abs(p_fwhm) / 2.35482;
    double c0 = 1 / pow(sigma*(2* M_PI), 0.5);
    double c1 = 1 / (0.5 / pow(sigma, 2));
    std::vector<double> intensity = p_intensity * c0 * exp(-c1 * pow(p_x - p_x0, 2));

    return intensity;
}

std::vector<double> PeakModel::Lorenzt() {

    double gamma = abs(p_fwhm) / 2;
    std::vector<double> intensity = p_intensity / (gamma * M_PI * (1 + pow((p_x - p_x0) / gamma, 2)));

    return intensity;
}