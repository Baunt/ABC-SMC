#include "util.h"
// #include "../include/third-party-library/matplotlib-cpp/matplotlibcpp.h"
#include "models/spectrum_model.h"
#include "json_handler.h"
#include "abc_smc_fit.h"
#include "models/simulatorparameters.h"

int main(int argc, char** argv) {

    Config config = JsonHandler().LoadConfigFile();
    JobData jobData = JsonHandler().LoadJobFile(config.jobFilePath);

    std::unique_ptr<SimulatorParameters> params = SimulatorParameters::ParameterFactory(jobData.modelType, jobData.inputFilePath);    
    std::unique_ptr<Simulator> simulator = Simulator::CreateSimulator(* params);
    
    AbcSmcFit optimizer(move(simulator)); 
    optimizer.Fit(jobData);
    return 0;
}
