//
// Created by balint.galgoczi on 2021.11.21..
//

#ifndef MATPLOTLIB_CPP_SPECTRUM_MODEL_H
#define MATPLOTLIB_CPP_SPECTRUM_MODEL_H


#include "probability_distribution.h"
#include "peak_model.h"
#include <map>
#include <list>

enum PeakType{
    Gauss,
    Lorentz
};

class SpectrumModel {
private:
    std::vector<PeakType> peaks;
    int npix;

public:
    std::vector<NormalDistribution> InitialGuess;
    Eigen::ArrayXX<double> GenerateInitialPopulation(int nsamples, int nparams);
    Eigen::ArrayX<double> Calculate(Eigen::ArrayX<double> parameters, bool withNoise);
    double ErrorCalculation(Eigen::ArrayX<double> diffSpectrum);
    void SetPeakList(std::vector<PeakType>& peaks);
    Eigen::ArrayX<double> energy;
    Eigen::ArrayX<double> intensity;


    SpectrumModel(const Eigen::ArrayX<double>& energy, const Eigen::ArrayX<double>& intensity){
        this->npix = energy.size();
        this->energy = energy;
        this->intensity = intensity;
    }
};


#endif //MATPLOTLIB_CPP_SPECTRUM_MODEL_H
