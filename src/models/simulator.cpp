#include "simulator.h"
#include "spectrum_model.h"
#include "example1.h"

std::unique_ptr<Simulator> Simulator::CreateSimulator(SimulatorParameters & params) {
    switch(params.simulatorType()){
        case SimulatorType::SPECTRUMMODEL:{            
            return std::make_unique<SpectrumModel>(static_cast<SpectrumModelParameters&>(params)); 
        }
        case SimulatorType::EXAMPLE1:            
            return std::make_unique<Example1> ();
        case SimulatorType::EXAMPLE2:
            return std::make_unique<Example1> ();
        default:
            return std::make_unique<Example1> ();
    }
}