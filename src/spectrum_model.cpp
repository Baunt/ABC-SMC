//
// Created by balint.galgoczi on 2021.09.14..
//
#define _USE_MATH_DEFINES

#include "spectrum_model.h"

std::vector<double> SpectrumModel::Gaussian()
    {
        std::vector<double> spectrum(p_npix, 0.0);
        double sigma = abs(p_fwhm) / 2.35482;
        double c0 = p_intensity / sigma * pow((2* M_PI), 0.5) ;
        double c1 = 0.5 / (sigma * sigma);

        for (int i = 0; i < p_npix; ++i)
        {
            spectrum[i] = c0 * exp(-c1 * (p_x[i] - p_x0) * (p_x[i] - p_x0));
        }
        return spectrum;
    }

std::vector<double> SpectrumModel::Lorentz() {
    std::vector<double> spectrum(p_npix, 0.0);

    double gamma = abs(p_fwhm) / 2;
    double inverseGamma = 1 / gamma;

    for (int i = 0; i < p_npix; ++i) {
        spectrum[i] = p_intensity / (gamma * M_PI * (1.0 + ((p_x[i] - p_x0) * inverseGamma) * ((p_x[i] - p_x0) * inverseGamma)));
    }

    return spectrum;
}