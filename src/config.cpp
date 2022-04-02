//
// Created by Bálint on 2/26/2022.
//

#include <string>

struct Config{
    std::string jobFilePath;
    std::string resultFilePath;
    std::string measuredSpectrumPath;
    int loggingLevel;
    bool simulated;
};