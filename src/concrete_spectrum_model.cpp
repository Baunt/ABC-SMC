//
// Created by balint.galgoczi on 2021.09.14..
//

#include "concrete_spectrum_model.h"

ConcreteSpectrumModel::ConcreteSpectrumModel(){
    this->Reset();
}

ConcreteSpectrumModel::~ConcreteSpectrumModel(){
    delete spectrumModel;
}

void ConcreteSpectrumModel::Reset(){
    this->spectrumModel = new SpectrumModel();
}

void ConcreteSpectrumModel::GaussianModel(){
    this->spectrumModel->Gaussian();
}

void ConcreteSpectrumModel::LorentzianModel(){
    this->spectrumModel->Lorentz();
}

SpectrumModel* ConcreteSpectrumModel::GetSpectrumModel(){
    SpectrumModel* result = this->spectrumModel;
    this->Reset();
    return result;
}

