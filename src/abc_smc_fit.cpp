//
// Created by balint.galgoczi on 2021.11.20..
//
#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include "abc_smc_fit.h"
#include "util.h"
#include "../include/third-party-library/Eigen/src/Cholesky/LLT.h"

void defineImportanceWeights(double beta, const Eigen::ArrayX<double>& ll_diffs, Eigen::ArrayX<double>& weights, double& newBeta, int rN){
    double lowBeta = beta;
    double oldBeta = beta;
    double upBeta = 2.0;

    while ((upBeta - lowBeta) > 1e-6) {
        double sum_of_square_weights = 0;
        newBeta = (lowBeta + upBeta) / 2.0;

        weights = exp((newBeta - oldBeta) * ll_diffs);
        weights /= weights.sum();
        sum_of_square_weights = weights.square().sum();

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
        weights = exp((newBeta - oldBeta) * ll_diffs);
        weights /= weights.sum();
    }
}

void runMcMcChains(int draws, int n_steps, double epsilon, double beta, int nparams,
                   Eigen::ArrayXX<double> &posteriors,
                   Eigen::ArrayX<double> &likelihoods,
                   Eigen::ArrayX<double> &tempered_logp,
                   Eigen::ArrayX<double> &acc_per_chain,
                   Eigen::ArrayX<double> &scalings,
                   Eigen::LLT<Eigen::MatrixXd> &llt,
                   Eigen::ArrayX<double> &prior_likelihoods,
                   SpectrumModel spectrumModel,
                   pcg32 &rng) {

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
            simulatedSpectrum = spectrumModel.Calculate(q_new);
            double priorLikelihood = spectrumModel.PriorLikelihood(q_new);
            double mean_abs_error = spectrumModel.ErrorCalculation(simulatedSpectrum);
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

void AbcSmcFit::Fit(SpectrumModel spectrumModel, int nparams, int draws, double epsilon, double threshold, double acc_rate, int n_steps, double p_acc_rate, bool tune_steps, double factor, int rngSeed) {
//    int nparams = 6;
//    int draws = 1000;
//    double epsilon = 0.01;
//    double threshold = 0.5;
//    double acc_rate = 1.0;
//    int n_steps = 25;
//    double p_acc_rate = 0.99;
//    bool tune_steps = true;
//    double factor = (2.38 * 2.38) / nparams;
    int proposed = draws * n_steps;

    int stage = 0;
    double beta = 0.0;

    pcg32 rng(rngSeed);
    Eigen::ArrayX<double> acc_per_chain(draws);
    Eigen::ArrayX<double> scalings(draws);
    Eigen::ArrayX<double> prior_likelihoods(draws);

    Eigen::ArrayXX<double> posteriors = spectrumModel.GenerateInitialPopulation(draws, nparams, rng);
    Eigen::ArrayX<double> likelihoods(draws);

    std::cout << "Starting population statistics: \n    ";

    if (factor < 1) {
        scalings = factor;
    } else {
        scalings = 1;
    }

    populationStatistics(posteriors);
    std::cout << "\n";

    for (int i = 0; i < draws; ++i) {
        Eigen::ArrayX<double> parameters = posteriors.row(i);
        prior_likelihoods(i) = spectrumModel.PriorLikelihood(parameters);
        double mean_abs_error = spectrumModel.ErrorCalculation(parameters);
        likelihoods(i) = (-(mean_abs_error * mean_abs_error) / (epsilon * epsilon) +
                          log(1.0 / (2.0 * M_PI * epsilon * epsilon))) / 2.0;
    }


    while (beta < 1) {
        std::cout << "Stage: " << std::setw(3) << stage << "    ";
        std::cout << "Beta: " << std::fixed << std::setprecision(6) << beta << "    ";
        std::cout << "Steps: " << std::setw(2) << n_steps << "    ";
        std::cout << "Acc. rate: " << std::fixed << std::setprecision(6) << acc_rate << std::endl;

#pragma region importance_weights
        int rN = int(round(likelihoods.size() * threshold));
        Eigen::ArrayX<double> ll_diffs = likelihoods - likelihoods.maxCoeff();
        Eigen::ArrayX<double> weights(ll_diffs.size());
        double newBeta = 0.0;

        defineImportanceWeights(beta, ll_diffs, weights, newBeta, rN);

        beta = newBeta;
#pragma endregion importance_weights

#pragma region resampling
        Eigen::ArrayX<int> resamplingIndexes = randomWeightedIndices(draws, weights, rng);

        posteriors = resampling_rows(posteriors, resamplingIndexes);
        prior_likelihoods = resampling(prior_likelihoods, resamplingIndexes);
        likelihoods = resampling(likelihoods, resamplingIndexes);
        scalings = resampling(scalings, resamplingIndexes);
        acc_per_chain = resampling(acc_per_chain, resamplingIndexes);

        Eigen::ArrayX<double> tempered_logp(draws);
        tempered_logp = prior_likelihoods + likelihoods * beta;

        Eigen::MatrixXd centered = posteriors.rowwise() - posteriors.colwise().mean();
        Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(posteriors.rows() - 1);
        Eigen::LLT<Eigen::MatrixXd> llt(cov);

#pragma endregion resampling

#pragma region tuning
        if (stage > 0) {
            double ave_scalings = exp(log(scalings.mean() + acc_per_chain.mean() - 0.234));

            scalings = 0.5 * (ave_scalings + exp(scalings.log() + (acc_per_chain - 0.234)));

            if (tune_steps) {
                acc_rate = std::max(1.0 / proposed, acc_rate);
                int t = (int) round((log(1.0 - p_acc_rate) / log(1.0 - acc_rate)));
                n_steps = std::min(n_steps, std::max(2, t));
            }
            proposed = draws * n_steps;
        }

#pragma endregion tuning

#pragma region MCMC

        runMcMcChains(draws,
                      n_steps,
                      epsilon,
                      beta,
                      nparams,
                      posteriors,
                      likelihoods,
                      tempered_logp,
                      acc_per_chain,
                      scalings,
                      llt,
                      prior_likelihoods,
                      spectrumModel,
                      rng);

        acc_rate = acc_per_chain.mean();

#pragma endregion MCMC

        stage += 1;
    }

    std::cout << "\nFinal population statistics: \n    ";
    populationStatistics(posteriors);

}