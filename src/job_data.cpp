//
// Created by Bálint on 2/26/2022.
//
#ifndef JOBDATA
#define JOBDATA

#include <string>

struct JobData{
    int draws;
    double epsilon;
    double threshold;
    int n_steps;
    double p_acc_rate;
    bool tune_steps;
    int rngSeed;
    std::string modelType;
    std::string inputFilePath;
};

#endif // JOBDATA