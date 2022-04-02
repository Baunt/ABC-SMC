//
// Created by Bálint on 2/26/2022.
//

#include "json_handler.h"
#include <iostream>
#include <fstream>
#include "../include/third-party-library/json/json.hpp"

using json = nlohmann::json;

Config JsonHandler::LoadConfigFile() {
    std::string path = "C:\\Balint\\Diplomamunka\\ABC_SMC\\configuration.json";
    std::fstream is(path);
    if (!is)
    {
        std::cout << "Cannot open " << path << std::endl;
    } else{
        json configurationFile;
        is >> configurationFile;

        Config Config;
        Config.jobFilePath = configurationFile["jobFilePath"];
        Config.loggingLevel = configurationFile["loggingLevel"];
        Config.measuredSpectrumPath = configurationFile["measuredSpectrumPath"];
        Config.resultFilePath = configurationFile["resultFilePath"];
        Config.simulated = configurationFile["simulated"];
        return Config;
    }
    return Config();
}

JobData JsonHandler::LoadJobFile(std::string jobFilePath) {
    std::fstream is(jobFilePath);
    if (!is)
    {
        std::cout << "Cannot open " << jobFilePath << std::endl;
    } else{

        json jobFile;
        is >> jobFile;

        JobData JobData;
        JobData.nparams = jobFile["nparams"];
        JobData.draws = jobFile["draws"];
        JobData.epsilon = jobFile["epsilon"];
        JobData.threshold = jobFile["threshold"];
        JobData.acc_rate = jobFile["acc_rate"];
        JobData.n_steps = jobFile["n_steps"];
        JobData.p_acc_rate = jobFile["p_acc_rate"];
        JobData.tune_steps = jobFile["tune_steps"];
        double factorNum = jobFile["factor"];
        double factor = factorNum / JobData.nparams;
        JobData.factor = factor;
        JobData.rngSeed = jobFile["rngSeed"];
        JobData.npix = jobFile["npix"];
        return JobData;
    }
    return JobData();
}

Spectrum JsonHandler::LoadMeasuredSpectrum(std::string measuredSpectrumPath) {
    std::fstream is(measuredSpectrumPath);
    if (!is)
    {
        std::cout << "Cannot open " << measuredSpectrumPath << std::endl;
    } else{

        json measuredFile;
        is >> measuredFile;

        Spectrum Data;
        json s = measuredFile["spectrum"];
        auto spectrumValues = s.get<std::vector<double>>();
        Eigen::ArrayX<double> spectrumArray(256);
        for (int i = 0; i < spectrumValues.size(); ++i) {
            spectrumArray[i] = spectrumValues[i];
        }
        Data.spectrum = spectrumArray;

        Data.peak = measuredFile[""]

        return Data;
    }
    return Spectrum();
}
