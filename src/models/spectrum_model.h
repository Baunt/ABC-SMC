//
// Created by balint.galgoczi on 2021.11.21..
//

#ifndef MATPLOTLIB_CPP_SPECTRUM_MODEL_H
#define MATPLOTLIB_CPP_SPECTRUM_MODEL_H

#include "simulator.h"
#include "../probability_distribution.h"
#include <map>
#include <list>


class SpectrumModel : public Simulator {
private:
    std::vector<PeakModel> peaks;
    int npix;

public:
    std::vector<NormalDistribution> InitialGuess;
    Eigen::ArrayXX<double> GenerateInitialPopulation(int nsamples, int nparams, pcg32 & rng);
    Eigen::ArrayX<double> Calculate(Eigen::ArrayX<double> parameters);
    double ErrorCalculation(Eigen::ArrayX<double> parameters);
    int NumberOfParameters();
    void SetPeakList(std::vector<PeakModel>& peaks);
    double PriorLikelihood(Eigen::ArrayX<double> posterior);
    Eigen::ArrayX<double> energy;
    Eigen::ArrayX<double> intensity;

    SpectrumModel(SpectrumModelParameters & params){
        this->energy = params.Measurement.spectrumWavelength;
        this->intensity = params.Measurement.spectrum;
        this->npix = params.Measurement.spectrumWavelength.size();
        this->SetPeakList(params.Measurement.peakModel);
    }
};


#endif //MATPLOTLIB_CPP_SPECTRUM_MODEL_H
