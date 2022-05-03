//
// Created by Bálint on 4/2/2022.
//
#ifndef SPECTRUM
#define SPECTRUM

#include <vector>
#include "../../include/third-party-library/Eigen/Core"
#include "peak_model.h"

class PeakModel;

struct Spectrum{
    Eigen::ArrayX<double> spectrum;
    Eigen::ArrayX<double> spectrumWavelength;
    std::vector<PeakModel> peakModel;
};

#endif // SPECTRUM