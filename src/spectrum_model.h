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
    std::list<PeakModel> peakModels;
public:
    int npix;
    std::vector<NormalDistribution> NormalDistributions;
    SpectrumModel * calculate(std::map<std::string, Eigen::ArrayX<double>> parameters, bool withNoise);
    Eigen::ArrayX<double> spectrum;
    Eigen::ArrayX<double> x;

    SpectrumModel(int npix){
        this->npix = npix;
        this->x = Eigen::VectorXd::LinSpaced(npix, 0, 1);
    }
};


#endif //MATPLOTLIB_CPP_SPECTRUM_MODEL_H
