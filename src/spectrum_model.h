//
// Created by balint.galgoczi on 2021.11.21..
//

#ifndef MATPLOTLIB_CPP_SPECTRUM_MODEL_H
#define MATPLOTLIB_CPP_SPECTRUM_MODEL_H


#include "probability_distribution.h"
#include "peak_model.h"
#include <map>
#include <list>


class SpectrumModel {
private:
public:
    int npix;
    std::vector<NormalDistribution> InitialGuess;
    Eigen::ArrayXX<double> GenerateInitialPopulation(int nsamples, int nparams);
    Eigen::ArrayX<double> calculate(std::map<std::string, Eigen::ArrayX<double>> parameters, bool withNoise);
    Eigen::ArrayX<double> x;

    // argumentumban megadott spektrummal initelni SpectrumModel(Eigen::ArrayX<double> energy, Eigen::ArrayX<double> intensity)
    SpectrumModel(int npix){
        this->npix = npix;
        this->x = Eigen::VectorXd::LinSpaced(npix, 0, 1);
    }
};


#endif //MATPLOTLIB_CPP_SPECTRUM_MODEL_H
