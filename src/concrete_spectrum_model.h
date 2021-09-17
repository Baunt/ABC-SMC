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
    void GaussianModel(spectrum_model_parameters modelParameters);
    void LorentzianModel(spectrum_model_parameters modelParameters);
};

#endif //MATPLOTLIB_CPP_CONCRETE_SPECTRUM_MODEL_H
