//
// Created by balint.galgoczi on 2021.07.28..
//
#define _USE_MATH_DEFINES

#include "peak_model.h"
#include "../include/third-party-library/Eigen/Core"
#include <cmath>
#include <iostream>

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

std::vector<double> divideVectorFromValue(double value, std::vector<double> vector){
    std::vector<double> dividedVector(vector.capacity(), value);
    for (int i = 0; i < vector.capacity(); ++i) {
        dividedVector[i] /= vector[i];
    }

    return dividedVector;
}


//std::vector<double> PeakModel::GaussianNonOptimized() {
//
//    double sigma = abs(p_fwhm) / 2.35482;
//    double c0 = 1 / (sigma * pow(2 * M_PI, 0.5));
//    double c1 = 0.5 / pow(sigma, 2);
//
//    std::vector<double> spectrum = subtractValueFromVector(p_x0, p_x);
//    spectrum = squaredVector(spectrum);
//    spectrum = multiplyValueFromVector(-c1, spectrum);
//    spectrum = exponentVector(spectrum);
//    spectrum = multiplyValueFromVector(p_intensity * c0, spectrum);
//
//    return spectrum;
//}

Eigen::ArrayX<double> PeakModel::Gaussian()
{
    Eigen::ArrayX<double> spectrum(p_npix);
    double sigma = std::abs(p_fwhm) / 2.35482;
    double c0 = p_intensity / (sigma * 2.5066283);
    double c1 = 0.5 / (sigma * sigma);

    for (int i = 0; i < p_npix; ++i)
    {
        spectrum(i) = c0 * std::exp(-c1 * (p_x[i] - p_x0) * (p_x[i] - p_x0));
    }

    return spectrum;
}

//std::vector<double> PeakModel::LorenztNonOptimized() {
//
//    double gamma = abs(p_fwhm) / 2;
//
//    std::vector<double> spectrum = subtractValueFromVector(p_x0, p_x);
//    spectrum = divideValueFromVector(gamma, spectrum);
//    spectrum = squaredVector(spectrum);
//    spectrum = addValueFromVector(1.0, spectrum);
//    spectrum = multiplyValueFromVector(gamma * M_PI, spectrum);
//    spectrum = divideVectorFromValue(p_intensity, spectrum);
//
//    return spectrum;
//}

Eigen::ArrayX<double> PeakModel::Lorenzt() {
    Eigen::ArrayX<double> spectrum(p_npix);
    double gamma = std::abs(p_fwhm) / 2;
    double inverseGamma = 1.0 / gamma;
    for (int i = 0; i < p_npix; ++i) {
        spectrum(i) = p_intensity / (gamma * M_PI * (1.0 + ((p_x[i] - p_x0) * inverseGamma) * ((p_x[i] - p_x0) * inverseGamma)));
    }

    return spectrum;
}