//
// Created by balint.galgoczi on 2021.07.28..
//

#ifndef ABC_SMC_ALGORITHM_PEAK_MODEL_H
#define ABC_SMC_ALGORITHM_PEAK_MODEL_H

#include <vector>

class PeakModel {
private:
    std::vector<double> p_x;
    double p_x0;
    double p_fwhm;
    double p_intensity;
public:
    std::vector<double> Gaussian();

    std::vector<double> Lorenzt();

    PeakModel(std::vector<double> x, double x0, double fwhm, double intensity)
    {
        p_x = x;
        p_x0 = x0;
        p_fwhm = fwhm;
        p_intensity = intensity;
    }
};


#endif //ABC_SMC_ALGORITHM_PEAK_MODEL_H
