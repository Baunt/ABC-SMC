//
// Created by balint.galgoczi on 2021.07.28..
//
#ifndef ABC_SMC_ALGORITHM_PEAK_MODEL_H
#define ABC_SMC_ALGORITHM_PEAK_MODEL_H

#include <vector>
#include <iostream>
#include "../../include/third-party-library/Eigen/Core"
// #include "spectrum_model.h"

enum PeakType{
    Gauss,
    Lorentz
};

class PeakModel {
public:
    PeakType peak;
    double x;
    double xUncertainty;
    double fwhm;
    double fwhmUncertainty;
    double intensity;
    double intensityUncertainty;

    PeakModel(double x, double xUncertainty, double fwhm, double fwhmUncertainty, double intensity, double intensityUncertainty, PeakType type)
    {
        this->x = x;
        this->xUncertainty = xUncertainty;
        this->fwhm = fwhm;
        this->fwhmUncertainty = fwhmUncertainty;
        this->intensity = intensity;
        this->intensityUncertainty = intensityUncertainty;
        this->peak = type;
    }
};


#endif //ABC_SMC_ALGORITHM_PEAK_MODEL_H
