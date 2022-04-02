//
// Created by Bálint on 4/2/2022.
//

#include <vector>
#include "../include/third-party-library/Eigen/Core"
#include "spectrum_model.h"

struct Spectrum{
    Eigen::ArrayX<double> spectrum;
    PeakType peak;
    double x;
    double xUncertainty;
    double fwhm;
    double fwhmUncertainty;
    double intensity;
    double intensityUncertainty;
};