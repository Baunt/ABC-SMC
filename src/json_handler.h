//
// Created by Bálint on 2/26/2022.
//

#ifndef ABC_SMC_ALGORITHM_CONFIG_FILE_HANDLER_H
#define ABC_SMC_ALGORITHM_CONFIG_FILE_HANDLER_H

#include "config.cpp"
#include "job_data.cpp"
#include "models/spectrum.cpp"

class JsonHandler {
public:
    Config LoadConfigFile();
    JobData LoadJobFile(std::string string);
    Spectrum LoadMeasuredSpectrum(std::string measuredSpectrumPath);
};


#endif //ABC_SMC_ALGORITHM_CONFIG_FILE_HANDLER_H
