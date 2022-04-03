//
// Created by balint.galgoczi on 2021.11.20..
//

#ifndef ABC_SMC_ALGORITHM_ABC_SMC_FIT_H
#define ABC_SMC_ALGORITHM_ABC_SMC_FIT_H

#include "../include/third-party-library/pcg-cpp/pcg_random.hpp"
#include "spectrum_model.h"

class AbcSmcFit {
public:
    void Fit(SpectrumModel spectrumModel, int nparams, int draws, double epsilon, double threshold, double acc_rate, int n_steps, double p_acc_rate, bool tune_steps, double factor, int rngSeed);
};


#endif //ABC_SMC_ALGORITHM_ABC_SMC_FIT_H
