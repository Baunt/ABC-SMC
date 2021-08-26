//
// Created by balint.galgoczi on 2021.07.28..
//
#define _USE_MATH_DEFINES

#include "peak_model.h"
#include <cmath>

std::vector<double> subtractValueFromVector(double value, std::vector<double> vector){
    for (int i = 0; i < vector.capacity(); ++i) {
        vector[i] -= value;
    }

    return vector;
}

std::vector<double> addValueFromVector(double value, std::vector<double> vector){
    for (int i = 0; i < vector.capacity(); ++i) {
        vector[i] += value;
    }

    return vector;
}


std::vector<double> squaredVector(std::vector<double> vector){
    for (int i = 0; i < vector.capacity(); ++i) {
        vector[i] = pow(vector[i], 2);
    }

    return vector;
}

std::vector<double> multiplyValueFromVector(double value, std::vector<double> vector){
    for (int i = 0; i < vector.capacity(); ++i) {
        vector[i] *= value;
    }

    return vector;
}

std::vector<double> exponentVector(std::vector<double> vector){
    for (int i = 0; i < vector.capacity(); ++i) {
        vector[i] = exp(vector[i]);
    }

    return vector;
}

std::vector<double> divideValueFromVector(double value, std::vector<double> vector){
    for (int i = 0; i < vector.capacity(); ++i) {
        vector[i] /= value;
    }

    return vector;
}

std::vector<double> PeakModel::Gaussian() {

    double sigma = abs(p_fwhm) / 2.35482;
    double c0 = 1 / pow(sigma*(2* M_PI), 0.5);
    double c1 = 1 / (0.5 / pow(sigma, 2));

    std::vector<double> spectrum = subtractValueFromVector(p_x0, p_x);
    spectrum = squaredVector(spectrum);
    spectrum = subtractValueFromVector(-c1, spectrum);
    spectrum = exponentVector(spectrum);
    spectrum = multiplyValueFromVector(p_intensity * c0, spectrum);

    return spectrum;
}

std::vector<double> PeakModel::Lorenzt() {

    double gamma = abs(p_fwhm) / 2;

    std::vector<double> spectrum = subtractValueFromVector(p_x0, p_x);
    spectrum = divideValueFromVector(gamma, spectrum);
    spectrum = squaredVector(spectrum);
    spectrum = addValueFromVector(1.0, spectrum);
    spectrum = multiplyValueFromVector(gamma * M_PI, spectrum);
    spectrum = divideValueFromVector(p_intensity, spectrum);

    return spectrum;
}