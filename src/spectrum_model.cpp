//
// Created by balint.galgoczi on 2021.11.21..
//

#include "spectrum_model.h"
#include "pcg_random_generator.h"
#include "util.h"

SpectrumModel *SpectrumModel::calculate(std::map<std::string, Eigen::ArrayX<double>> parameters, bool withNoise) {
    std::map<std::string, Eigen::ArrayX<double>>::iterator itr;
    Eigen::ArrayX<double> noise = getDistribution(0.0, 1.0, npix);
    Eigen::ArrayX<double> spectrum(npix);

    for (itr = parameters.begin(); itr != parameters.end(); ++itr) {
        if (itr->first == "Gaussian"){
            PeakModel gauss = PeakModel(x, itr->second[0], itr->second[1], itr->second[2], npix, "Gaussian");
            peakModels.push_back(gauss);

            spectrum += gauss.Gaussian();
        }
        else if(itr->first == "Lorentzian"){
            PeakModel lorentz = PeakModel(x, itr->second[0], itr->second[1], itr->second[2], npix, "Lorentzian");
            peakModels.push_back(lorentz);

            spectrum += lorentz.Lorenzt();
        }
        else
            spectrum.setZero();
    }

    if (withNoise){
        spectrum += noise;
    }

    NormalDistribution x0 = NormalDistribution("x0", "normal", 0.63, 0.15);
    NormalDistribution fwhm0 = NormalDistribution("fwhm0", "normal", 0.09, 0.04);
    NormalDistribution int0 = NormalDistribution("int0", "normal", 2.65, 0.5);
    NormalDistribution x1 = NormalDistribution("x1", "normal", 0.45, 0.30);
    NormalDistribution fwhm1 = NormalDistribution("fwhm1", "normal", 0.25, 0.07);
    NormalDistribution int1 = NormalDistribution("int1", "normal", 0.35, 0.15);

    normalDistribution.push_back(x0);
    normalDistribution.push_back(fwhm0);
    normalDistribution.push_back(int0);
    normalDistribution.push_back(x1);
    normalDistribution.push_back(fwhm1);
    normalDistribution.push_back(int1);


    return this;
}
