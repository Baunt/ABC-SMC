#include "simulator.h"
#include "simulatorparameters.h"
#include "spectrum_model.h"
#include "../json_handler.h"

SimulatorType resolveSimulatorType(std::string input) {
    if( input == "SpectrumModel" ) return SimulatorType::SPECTRUMMODEL;
    if( input == "Example1" ) return SimulatorType::EXAMPLE1;
    if( input == "Example2" ) return SimulatorType::EXAMPLE2;
    return NOTIMPLEMENTED;
}

std::unique_ptr<SimulatorParameters> SimulatorParameters::ParameterFactory(std::string modelType, std::string inputFilePath) {
    switch (resolveSimulatorType(modelType)){
        case SimulatorType::SPECTRUMMODEL: {
            auto parameters = std::make_unique<SpectrumModelParameters>();
            parameters->Measurement = JsonHandler().LoadMeasuredSpectrum(inputFilePath);
            return parameters;
        }
        case SimulatorType::EXAMPLE1:
            return std::make_unique<Example1Parameters>();
        case SimulatorType::EXAMPLE2:
            return std::make_unique<Example1Parameters>();
        case SimulatorType::NOTIMPLEMENTED:
            return NULL;
    }    
    return NULL;
}
