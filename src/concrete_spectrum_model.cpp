//
// Created by balint.galgoczi on 2021.09.14..
//

#include "concrete_spectrum_model.h"

#include <utility>

ConcreteSpectrumModel::ConcreteSpectrumModel(){
    this->Reset();
}

ConcreteSpectrumModel::~ConcreteSpectrumModel(){
    delete spectrumModel;
}

void ConcreteSpectrumModel::Reset(){
    this->spectrumModel = new SpectrumModel();
}

void ConcreteSpectrumModel::GaussianModel(spectrum_model_parameters modelParameters){
    this->spectrumModel->Gaussian(std::move(modelParameters));
}

void ConcreteSpectrumModel::LorentzianModel(spectrum_model_parameters modelParameters){
    this->spectrumModel->Lorentz(std::move(modelParameters));
}

SpectrumModel* ConcreteSpectrumModel::GetSpectrumModel(){
    SpectrumModel* result = this->spectrumModel;
    this->Reset();
    return result;
}

