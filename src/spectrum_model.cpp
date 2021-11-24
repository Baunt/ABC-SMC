//
// Created by balint.galgoczi on 2021.11.21..
//
#define _USE_MATH_DEFINES
#include "spectrum_model.h"
#include "pcg_random_generator.h"
#include "util.h"

Eigen::ArrayX<double> lorenzt(const Eigen::ArrayX<double>& x, double x0, double fwhm, double intensity, int npix) {
    Eigen::ArrayX<double> spectrum(npix);
    double gamma = std::abs(fwhm) / 2;
    double inverseGamma = 1.0 / gamma;

    spectrum = intensity / (gamma * M_PI * (1.0 + ((x - x0) * inverseGamma) * ((x - x0) * inverseGamma)));

    return spectrum;
}

Eigen::ArrayX<double> gaussian(const Eigen::ArrayX<double>& x, double x0, double fwhm, double intensity, int npix)
{
    Eigen::ArrayX<double> spectrum(npix);
    double sigma = std::abs(fwhm) / 2.35482;
    double c0 = intensity / (sigma * 2.5066283);
    double c1 = 0.5 / (sigma * sigma);

    spectrum = c0 * exp(-c1 * (x - x0) * (x - x0));

    return spectrum;
}
Eigen::ArrayX<double> SpectrumModel::calculate(std::map<std::string, Eigen::ArrayX<double>> parameters, bool withNoise) {
    std::map<std::string, Eigen::ArrayX<double>>::iterator itr;
    Eigen::ArrayX<double> noise = getDistribution(0.0, 1.0, npix);
    Eigen::ArrayX<double> internalSpectrum(npix);
    internalSpectrum.setZero();
    
    for (itr = parameters.begin(); itr != parameters.end(); ++itr) {
        if (itr->first == "Gaussian"){
            internalSpectrum += gaussian(x, itr->second[0], itr->second[1], itr->second[2], npix);
        }
        else if(itr->first == "Lorentzian"){
            internalSpectrum += lorenzt(x, itr->second[0], itr->second[1], itr->second[2], npix);
        }
        else
            internalSpectrum.setZero();
    }

    if (withNoise){
        internalSpectrum += noise;
    }

    return internalSpectrum;
}

Eigen::ArrayXX<double> SpectrumModel::GenerateInitialPopulation(int nsamples, int nparams) {

    Eigen::ArrayXX<double> posteriors(nsamples , nparams);

    NormalDistribution x0 = NormalDistribution("x0", "normal", 0.63, 0.15);
    NormalDistribution fwhm0 = NormalDistribution("fwhm0", "normal", 0.09, 0.04);
    NormalDistribution int0 = NormalDistribution("int0", "normal", 2.65, 0.5);
    NormalDistribution x1 = NormalDistribution("x1", "normal", 0.45, 0.30);
    NormalDistribution fwhm1 = NormalDistribution("fwhm1", "normal", 0.25, 0.07);
    NormalDistribution int1 = NormalDistribution("int1", "normal", 0.35, 0.15);

    InitialGuess.push_back(x0);
    InitialGuess.push_back(fwhm0);
    InitialGuess.push_back(int0);
    InitialGuess.push_back(x1);
    InitialGuess.push_back(fwhm1);
    InitialGuess.push_back(int1);

    for (int i = 0; i < InitialGuess.size(); ++i) {
        Eigen::ArrayX<double> tmp = InitialGuess[i].Sample(nsamples);
        posteriors.col(i) = tmp;
    }

    return posteriors;
}


