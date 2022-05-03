#ifndef SIMULATORPARAMETERS
#define SIMULATORPARAMETERS

// #include "simulator.h"
// #include "spectrum_model.h"
#include "spectrum.cpp"
#include <memory>

enum SimulatorType {
    SPECTRUMMODEL = 1,
    EXAMPLE1 = 2,
    EXAMPLE2 = 3,
    NOTIMPLEMENTED = 99
};

struct SimulatorParameters {
    virtual SimulatorType simulatorType() const = 0;
    virtual ~SimulatorParameters() = default;
    static std::unique_ptr<SimulatorParameters> ParameterFactory(std::string modelType, std::string inputFilePath);
};

struct SpectrumModelParameters : public SimulatorParameters {
    SimulatorType simulatorType() const {return SimulatorType::SPECTRUMMODEL; };
    Spectrum Measurement;
};

struct Example1Parameters: public SimulatorParameters {
    SimulatorType simulatorType() const {return SimulatorType::EXAMPLE1; };

};

struct Example2Parameters: public SimulatorParameters {
    SimulatorType simulatorType() const {return SimulatorType::EXAMPLE2; };

};


#endif // SIMULATORPARAMETERS