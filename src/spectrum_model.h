//
// Created by balint.galgoczi on 2021.09.14..
//

#ifndef MATPLOTLIB_CPP_SPECTRUM_MODEL_H
#define MATPLOTLIB_CPP_SPECTRUM_MODEL_H


#include <vector>

class SpectrumModel{
private:
    std::vector<double> p_x;
    double p_x0;
    double p_fwhm;
    double p_intensity;
    double p_npix;
public:
    std::vector<double> Gaussian();
    std::vector<double> Lorentz();
};

#endif //MATPLOTLIB_CPP_SPECTRUM_MODEL_H
