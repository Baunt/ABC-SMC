//
// Created by balint.galgoczi on 2021.11.20..
//

#ifndef ABC_SMC_ALGORITHM_ABC_SMC_FIT_H
#define ABC_SMC_ALGORITHM_ABC_SMC_FIT_H


#include "../include/third-party-library/pcg-cpp/pcg_random.hpp"
#include "spectrum_model.h"

class AbcSmcFit {
public:
    void Fit(SpectrumModel spectrumModel, const Eigen::ArrayX<double>& spectrum);
};


#endif //ABC_SMC_ALGORITHM_ABC_SMC_FIT_H
