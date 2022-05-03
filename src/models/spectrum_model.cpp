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

Eigen::ArrayX<double> gaussian(const Eigen::ArrayX<double>& x, double x0, double fwhm, double intensity, int npix) {
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
        PeakType actualPeak = peaks[i].peak;
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
    
    for (auto peakModel:peaks) {
        NormalDistribution xNormalDistribution = NormalDistribution(peakModel.x, peakModel.xUncertainty);
        NormalDistribution fwhmNormalDistribution = NormalDistribution(peakModel.fwhm, peakModel.fwhmUncertainty);
        NormalDistribution intNormalDistribution = NormalDistribution(peakModel.intensity, peakModel.intensityUncertainty);
        InitialGuess.push_back(xNormalDistribution);
        InitialGuess.push_back(fwhmNormalDistribution);
        InitialGuess.push_back(intNormalDistribution);
    }
    
    for (int i = 0; i < InitialGuess.size(); ++i) {
        Eigen::ArrayX<double> tmp = InitialGuess[i].Sample(nsamples, rng);
        priors.col(i) = tmp;
    }

    return priors;
}

void SpectrumModel::SetPeakList(std::vector<PeakModel>& peakList) {
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

int SpectrumModel::NumberOfParameters() {
    return peaks.size() * 3;
}