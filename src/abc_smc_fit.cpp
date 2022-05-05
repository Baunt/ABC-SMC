//
// Created by balint.galgoczi on 2021.11.20..
//
#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <utility>
#include "abc_smc_fit.h"
#include "util.h"
#include "models/peak_model.h"
#include "models/simulator.h"
#include "../include/third-party-library/Eigen/src/Cholesky/LLT.h"


void AbcSmcFit::initializeArrays(double init_scalings){
    acc_per_chain.resize(draws, Eigen::NoChange);
    scalings.resize(draws, Eigen::NoChange);
    prior_likelihoods.resize(draws, Eigen::NoChange);
    tempered_logp.resize(draws, Eigen::NoChange);
    likelihoods.resize(draws, Eigen::NoChange);
    importance_weights.resize(draws, Eigen::NoChange);
    posteriors.resize(draws, Eigen::NoChange);
    
    acc_per_chain = 1.0;
    scalings = init_scalings;
}
    
double AbcSmcFit::determineImportanceWeights(double threshold){
    int rN = int(round(likelihoods.size() * threshold));        
    double newBeta = 0.0;
    double lowBeta = beta;
    double oldBeta = beta;
    double upBeta = 2.0;

    Eigen::ArrayX<double> ll_diffs = likelihoods - likelihoods.maxCoeff();

    while ((upBeta - lowBeta) > 1e-6) {
        double sum_of_square_weights = 0;
        newBeta = (lowBeta + upBeta) / 2.0;

        importance_weights = exp((newBeta - oldBeta) * ll_diffs);
        importance_weights /= importance_weights.sum();
        sum_of_square_weights = importance_weights.square().sum();

        int ESS = int(round(1.0 / sum_of_square_weights));
        if (ESS == rN) {
            break;
        } else if (ESS < rN) {
            upBeta = newBeta;
        } else {
            lowBeta = newBeta;
        }
    }

    if (newBeta >= 1) {
        double sum_of_weights_un = 0;
        newBeta = 1;
        importance_weights = exp((newBeta - oldBeta) * ll_diffs);
        importance_weights /= importance_weights.sum();
    }

    beta = newBeta;
    return beta;
}

void AbcSmcFit::calculateInitialLikelihoods(){
    for (int i = 0; i < draws; ++i) {
        Eigen::ArrayX<double> parameters = posteriors.row(i);
        prior_likelihoods(i) = model->PriorLikelihood(parameters);
        double mean_abs_error = model->ErrorCalculation(parameters);
        likelihoods(i) = (-(mean_abs_error * mean_abs_error) / (epsilon * epsilon) +
                          log(1.0 / (2.0 * M_PI * epsilon * epsilon))) / 2.0;
    }
}

void AbcSmcFit::runMcMcChains(Eigen::LLT<Eigen::MatrixXd> &llt, pcg32 &rng) {
    std::normal_distribution<double> std_dist(0.0, 1.0); // normal distribution with mu=0, sigma=1
    std::uniform_real_distribution<double> uni_dist(0.0, 1.0); // uniform ditribution between: 0 <= r < 1
    auto rand_fn = [&]() { return std_dist(rng); };// lambda that generates random numbers

    for (int draw = 0; draw < draws; ++draw) {
        Eigen::MatrixXd random_matrix = Eigen::MatrixXd::NullaryExpr(nparams, n_steps, rand_fn);// matrix expression
        Eigen::MatrixXd deltas = llt.matrixLLT() * random_matrix;
        deltas *= scalings[draw];

        int accepted = 0;

        double old_likelihood = likelihoods[draw];
        double old_prior = prior_likelihoods[draw];
        double old_tempered_logp = tempered_logp[draw];
        Eigen::ArrayX<double> q_old = posteriors.row(draw);
        Eigen::ArrayX<double> simulatedSpectrum;

        for (int i = 0; i < n_steps; ++i) {
            Eigen::ArrayX<double> delta = deltas.col(i);
            Eigen::ArrayX<double> q_new = q_old + delta;
            double priorLikelihood = model->PriorLikelihood(q_new);
            double mean_abs_error  = model->ErrorCalculation(q_new);
            double likelihood = (-(mean_abs_error * mean_abs_error) / (epsilon * epsilon) +
                                 log(1.0 / (2.0 * M_PI * epsilon * epsilon))) / 2.0;
            double new_tempered_logp = priorLikelihood + likelihood * beta;
            double delta_tempered_logp = new_tempered_logp - tempered_logp[draw];

            double r = log(uni_dist(rng));
            if (r < delta_tempered_logp) {
                q_old = q_new;
                accepted += 1;
                old_prior = priorLikelihood;
                old_likelihood = likelihood;
                old_tempered_logp = new_tempered_logp;
            }
        }
        prior_likelihoods(draw) = old_prior;
        likelihoods(draw) = old_likelihood;
        tempered_logp(draw) = old_tempered_logp;
        acc_per_chain(draw) = double(accepted) / double(n_steps);
        posteriors.row(draw) = q_old;

        for (int i = 0; i < nparams; ++i) {
            posteriors(draw, i) = q_old[i];
        }
    }
}

void AbcSmcFit::resample(pcg32 & rng){
        Eigen::ArrayX<int> resamplingIndexes = randomWeightedIndices(draws, importance_weights, rng);

        posteriors = resampling_rows(posteriors, resamplingIndexes);
        prior_likelihoods = resampling(prior_likelihoods, resamplingIndexes);
        likelihoods = resampling(likelihoods, resamplingIndexes);
        scalings = resampling(scalings, resamplingIndexes);
        acc_per_chain = resampling(acc_per_chain, resamplingIndexes);

        tempered_logp = prior_likelihoods + likelihoods * beta;
}

void AbcSmcFit::tuning(bool tune_steps, int & proposed, double & acc_rate, double p_acc_rate){
    double ave_scalings = exp(log(scalings.mean() + acc_per_chain.mean() - 0.234));

    scalings = 0.5 * (ave_scalings + exp(scalings.log() + (acc_per_chain - 0.234)));

    if (tune_steps) {
        acc_rate = std::max(1.0 / proposed, acc_rate);
        int t = (int) round((log(1.0 - p_acc_rate) / log(1.0 - acc_rate)));
        n_steps = std::min(n_steps, std::max(2, t));
    }
    proposed = draws * n_steps;
}

void AbcSmcFit::Fit(JobData job) {    
    // AbcSmcFit members
    nparams = model->NumberOfParameters();
    draws = job.draws;
    n_steps = job.n_steps;
    epsilon = job.epsilon;
    beta = 0.0;

    // locals
    pcg32 rng(job.rngSeed);
    double threshold = job.threshold;
    double acc_rate = 1.0;    
    double p_acc_rate = 0.99;
    bool tune_steps = job.tune_steps;
    int proposed = draws * n_steps;
    int stage = 0;


    initializeArrays(std::min(1.0, (2.38 * 2.38) / nparams));
    posteriors = model->GenerateInitialPopulation(draws, nparams, rng); 
    calculateInitialLikelihoods();


    std::cout << "Starting population statistics: \n    ";
    populationStatistics(posteriors);
    std::cout << "\n";


    while (beta < 1) {
        std::cout << "Stage: " << std::setw(3) << stage << "    ";
        std::cout << "Beta: " << std::fixed << std::setprecision(6) << beta << "    ";
        std::cout << "Steps: " << std::setw(2) << n_steps << "    ";
        std::cout << "Acc. rate: " << std::fixed << std::setprecision(6) << acc_rate << std::endl;

        beta = determineImportanceWeights(threshold);
        resample(rng);

        // calculate the covariance matrix, and it's Cholesky decomposition
        Eigen::MatrixXd centered = posteriors.rowwise() - posteriors.colwise().mean();
        Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(posteriors.rows() - 1);
        Eigen::LLT<Eigen::MatrixXd> llt(cov);

        if (stage > 0) {
            tuning(tune_steps, proposed, acc_rate, p_acc_rate);
        }

        runMcMcChains(llt, rng);

        acc_rate = acc_per_chain.mean();
        stage += 1;
    }

    std::cout << "\nFinal population statistics: \n    ";
    populationStatistics(posteriors);

}