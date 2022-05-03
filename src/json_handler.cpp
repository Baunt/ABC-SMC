//
// Created by Bálint on 2/26/2022.
//

#include "json_handler.h"
#include <iostream>
#include <fstream>
#include "../include/third-party-library/json/json.hpp"

using json = nlohmann::json;

Config JsonHandler::LoadConfigFile() {
    std::string path = "./configuration.json";
    std::fstream is(path);
    if (!is) {
        std::cout << "Cannot open " << path << std::endl;
    } else {
        json configurationFile;
        is >> configurationFile;

        Config Config;
        Config.jobFilePath = configurationFile["jobFilePath"];
        Config.loggingLevel = configurationFile["loggingLevel"];
        Config.resultFilePath = configurationFile["resultFilePath"];
        Config.simulated = configurationFile["simulated"];
        return Config;
    }
    return Config();
}

JobData JsonHandler::LoadJobFile(std::string jobFilePath) {
    std::fstream is(jobFilePath);
    if (!is) {
        std::cout << "Cannot open " << jobFilePath << std::endl;
    } else {

        json jobFile;
        is >> jobFile;

        JobData JobData;
        JobData.draws = jobFile["draws"];
        JobData.epsilon = jobFile["epsilon"];
        JobData.threshold = jobFile["threshold"];
        JobData.n_steps = jobFile["n_steps"];
        JobData.p_acc_rate = jobFile["p_acc_rate"];
        JobData.tune_steps = jobFile["tune_steps"];
        JobData.rngSeed = jobFile["rngSeed"];
        JobData.modelType = jobFile["modelType"];
        JobData.inputFilePath = jobFile["measuredSpectrumPath"];

        return JobData;
    }
    return JobData();
}

Spectrum JsonHandler::LoadMeasuredSpectrum(std::string measuredSpectrumPath) {
    std::fstream is(measuredSpectrumPath);
    if (!is) {
        std::cout << "Cannot open " << measuredSpectrumPath << std::endl;
    } else {

        json measuredFile;
        is >> measuredFile;

        Spectrum Data;
        json s = measuredFile["spectrum"];
        json sw = measuredFile["spectrumWavelength"];

        auto spectrumValues = s.get<std::vector<double>>();
        auto spectrumWavelengthValues = sw.get<std::vector<double>>();

        if (spectrumValues.size() != spectrumWavelengthValues.size())
        {
            std::cout << "ERROR: Spectrum intensity and spectrum energy vector size must be equal! " << std::endl;
        }

        Eigen::ArrayX<double> spectrumArray(spectrumValues.size());
        for (int i = 0; i < spectrumValues.size(); ++i) {
            spectrumArray[i] = spectrumValues[i];
        }
        Data.spectrum = spectrumArray;

        Eigen::ArrayX<double> spectrumWavelengthArray(spectrumValues.size());
        for (int i = 0; i < spectrumWavelengthValues.size(); ++i) {
            spectrumWavelengthArray[i] = spectrumWavelengthValues[i];
        }
        Data.spectrumWavelength = spectrumWavelengthArray;

        std::vector<PeakModel> peakModelVector;
        auto peaksGuess = measuredFile["peaksGuess"];
        for (auto &peaks: peaksGuess.items()) {
            std::vector<double> x;
            std::vector<double> fwhm;
            std::vector<double> intensity;
            PeakType peakFromJson;

            for (auto &peak: peaks.value().items()) {
                auto peakKey = peak.key();
                if (peakKey == "x") {
                    x.clear();
                    x.push_back(peak.value()[0]);
                    x.push_back(peak.value()[1]);
                } else if (peakKey == "fwhm") {
                    fwhm.clear();
                    fwhm.push_back(peak.value()[0]);
                    fwhm.push_back(peak.value()[1]);
                } else if (peakKey == "int") {
                    intensity.clear();
                    intensity.push_back(peak.value()[0]);
                    intensity.push_back(peak.value()[1]);
                } else if (peakKey == "peak") {
                    if (peak.value() == "Gauss") {
                        PeakType gauss = Gauss;
                        peakFromJson = gauss;
                    } else if (peak.value() == "Lorentz") {
                        PeakType lorentz = Lorentz;
                        peakFromJson = lorentz;
                    } else {
                        //TODO error handling
                    }
                } else {
                    //TODO error handling
                }
            }
            PeakModel pm = PeakModel(x[0], x[1], fwhm[0], fwhm[1], intensity[0], intensity[1], peakFromJson);
            peakModelVector.push_back(pm);
        }

        Data.peakModel = peakModelVector;
        return Data;
    }
    return Spectrum();
}
