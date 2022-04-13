//
// Created by Bálint on 4/2/2022.
//

#include <vector>
#include "../include/third-party-library/Eigen/Core"
#include "peak_model.h"

struct Spectrum{
    Eigen::ArrayX<double> spectrum;
    Eigen::ArrayX<double> spectrumWavelength;
    std::vector<PeakModel> peakModel;
};