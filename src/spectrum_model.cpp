//
// Created by balint.galgoczi on 2021.11.21..
//
#define _USE_MATH_DEFINES
#include "spectrum_model.h"
#include <utility>

Eigen::ArrayX<double> lorentz(const Eigen::ArrayX<double>& x, double x0, double fwhm, double intensity, int npix) {
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

Eigen::ArrayX<double> SpectrumModel::Calculate(Eigen::ArrayX<double> parameters) {
    Eigen::ArrayX<double> internalSpectrum(npix);
    internalSpectrum.setZero();

    for (int i = 0; i < peaks.size(); ++i) {
        PeakType actualPeak = peaks[i];
        switch (actualPeak) {
            case Gauss:
            {
                Eigen::ArrayX<double> gaussParametersSegment = parameters.segment(i * 3, 3);
                internalSpectrum += gaussian(energy, gaussParametersSegment[0], gaussParametersSegment[1], gaussParametersSegment[2], npix);
                break;
            }
            case Lorentz:
            {
                Eigen::ArrayX<double> lorentzParametersSegment = parameters.segment(i * 3, 3);
                internalSpectrum += lorentz(energy, lorentzParametersSegment[0], lorentzParametersSegment[1],
                                            lorentzParametersSegment[2], npix);
                break;
            }
            default:
                internalSpectrum.setZero();
                break;
        }
    }

    return internalSpectrum;
}

Eigen::ArrayXX<double> SpectrumModel::GenerateInitialPopulation(int nsamples, int nparams, pcg32 & rng) {

    Eigen::ArrayXX<double> priors(nsamples , nparams);

    NormalDistribution x0 = NormalDistribution(0.63, 0.15);
    NormalDistribution fwhm0 = NormalDistribution(0.09, 0.04);
    NormalDistribution int0 = NormalDistribution(2.65, 0.5);
    NormalDistribution x1 = NormalDistribution(0.45, 0.30);
    NormalDistribution fwhm1 = NormalDistribution(0.25, 0.07);
    NormalDistribution int1 = NormalDistribution(0.35, 0.15);

    InitialGuess.push_back(x0);
    InitialGuess.push_back(fwhm0);
    InitialGuess.push_back(int0);
    InitialGuess.push_back(x1);
    InitialGuess.push_back(fwhm1);
    InitialGuess.push_back(int1);

    for (int i = 0; i < InitialGuess.size(); ++i) {
        Eigen::ArrayX<double> tmp = InitialGuess[i].Sample(nsamples, rng);
        priors.col(i) = tmp;
    }

    return priors;
}

void SpectrumModel::SetPeakList(std::vector<PeakType>& peakList) {
    this->peaks = peakList;
}

double SpectrumModel::ErrorCalculation(Eigen::ArrayX<double> parameters) {
    return ((intensity - SpectrumModel::Calculate(parameters)).cwiseAbs()).mean();
}

double SpectrumModel::PriorLikelihood(Eigen::ArrayX<double> posterior) {

    double prior_likelihoods = 0.0;
    for (int i = 0; i < posterior.size(); ++i) {
    prior_likelihoods += InitialGuess[i].LogP(posterior(i));
    }

    return prior_likelihoods;
}

