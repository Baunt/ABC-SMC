//
// Created by balint.galgoczi on 2021.09.14..
//

#ifndef MATPLOTLIB_CPP_SPECTRUM_MODEL_H
#define MATPLOTLIB_CPP_SPECTRUM_MODEL_H

#include <vector>

struct spectrum_model_parameters{
    std::vector<double> x;
    double x0;
    double fwhm;
    double intensity;
    double npix;
}SpectrumModelParameters;


class SpectrumModel{
public:
    std::vector<double> Gaussian(spectrum_model_parameters modelParameters);
    std::vector<double> Lorentz(spectrum_model_parameters modelParameters);
};

#endif //MATPLOTLIB_CPP_SPECTRUM_MODEL_H
