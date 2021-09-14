//
// Created by balint.galgoczi on 2021.09.14..
//

#ifndef MATPLOTLIB_CPP_CONCRETE_SPECTRUM_MODEL_H
#define MATPLOTLIB_CPP_CONCRETE_SPECTRUM_MODEL_H

#include "spectrum_model.h"

class ConcreteSpectrumModel{
private:
    SpectrumModel* spectrumModel{};
public:
    ConcreteSpectrumModel();
    ~ConcreteSpectrumModel();
    void Reset();
    SpectrumModel* GetSpectrumModel();
    void GaussianModel();
    void LorentzianModel();
};

#endif //MATPLOTLIB_CPP_CONCRETE_SPECTRUM_MODEL_H
