//
// Created by balint.galgoczi on 2021.09.14..
//
#define _USE_MATH_DEFINES

#include "spectrum_model.h"

std::vector<double> SpectrumModel::Gaussian(spectrum_model_parameters modelParameters)
    {
        std::vector<double> spectrum(modelParameters.npix, 0.0);
        double sigma = abs(modelParameters.fwhm) / 2.35482;
        double c0 = modelParameters.intensity / sigma * pow((2* M_PI), 0.5) ;
        double c1 = 0.5 / (sigma * sigma);

        for (int i = 0; i < modelParameters.npix; ++i)
        {
            spectrum[i] = c0 * exp(-c1 * (modelParameters.x[i] - modelParameters.x0) * (modelParameters.x[i] - modelParameters.x0));
        }
        return spectrum;
    }

std::vector<double> SpectrumModel::Lorentz(spectrum_model_parameters modelParameters) {
    std::vector<double> spectrum(modelParameters.npix, 0.0);

    double gamma = abs(modelParameters.fwhm) / 2;
    double inverseGamma = 1 / gamma;

    for (int i = 0; i < modelParameters.npix; ++i) {
        spectrum[i] = modelParameters.intensity / (gamma * M_PI * (1.0 + ((modelParameters.x[i] - modelParameters.x0) * inverseGamma) * ((modelParameters.x[i] - modelParameters.x0) * inverseGamma)));
    }

    return spectrum;
}