//
// Created by balint.galgoczi on 2021.11.20..
//
#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include "abc_smc_fit.h"
#include "util.h"
#include "peak_model.h"
#include "../include/third-party-library/Eigen/src/Cholesky/LLT.h"
#include "pcg_random_generator.h"

bool DEBUG;

void runMcMcChains(int draws, int n_steps, double epsilon, double beta, int nparams,
                   Eigen::ArrayXX<double>& posteriors,
                   Eigen::ArrayX<double>& likelihoods,
                   Eigen::ArrayX<double>& tempered_logp,
                   Eigen::ArrayX<double>& acc_per_chain,
                   Eigen::ArrayX<double>& scalings,
                   Eigen::LLT<Eigen::MatrixXd>& llt,
                   Eigen::ArrayX<double>& prior_likelihoods,
                   SpectrumModel simulatedSpectrumModel,
                   SpectrumModel spectrumModel){
    for (int draw = 0; draw < draws; ++draw) {
        if (DEBUG)
        {
            std::cout << "  Sample " << draw << "\n";
            std::cout << posteriors(draw, 0) << " " << posteriors(draw, 1) << " " << posteriors(draw, 2) << " " << posteriors(draw, 3) << " " << posteriors(draw, 4) << " " << posteriors(draw, 5) << "\n";
            std::cout << "    err        pl            ll          beta     logp_0    logp_1          d_logp         r\n";
        }

        //TODO melyik legyen?
        //Eigen::MatrixXd randomsForDeltas = Eigen::MatrixXd::Random(n_steps, 6);
        pcg32 rng = PcgRandomGenerator::GetInstance()->value();
        std::normal_distribution<double> std_dist(0.0, 1.0); // normal distribution with mu=0, sigma=1
        std::uniform_real_distribution<double> uni_dist(0.0, 1.0); // uniform ditribution between: 0 <= r < 1
        auto rand_fn = [&](){return std_dist(rng);};// lambda that generates random numbers
        Eigen::MatrixXd random_matrix = Eigen::MatrixXd::NullaryExpr(n_steps, 6, rand_fn);// matrix expression

        Eigen::MatrixXd randomsForDeltas = random_matrix;// saving the elements to a matrix
        Eigen::MatrixXd deltas = randomsForDeltas * llt.matrixLLT();
        deltas *= scalings[draw];

        int accepted = 0;

        double old_likelihood = likelihoods[draw];
        double old_prior = prior_likelihoods[draw];
        double old_tempered_logp = tempered_logp[draw];
        Eigen::ArrayX<double> q_old = posteriors.row(draw);

        for (int i = 0; i < n_steps; ++i) {
            Eigen::ArrayX<double> delta = deltas.row(i);
            Eigen::ArrayX<double> q_new = q_old + delta;
            simulatedSpectrumModel.spectrum = staticPeakModel(spectrumModel.x, q_new);
            Eigen::ArrayX<double> spectrumDiffs(spectrumModel.npix);

            spectrumDiffs = abs(spectrumModel.spectrum - simulatedSpectrumModel.spectrum);

            double pl = 0.0;
            //TODO
            for (int dist = 0; dist < spectrumModel.NormalDistributions.capacity(); ++dist) {
                pl += spectrumModel.NormalDistributions[dist].LogP(q_new[dist]);
            }
            double mean_abs_error = spectrumDiffs.mean();
            double ll = (-(mean_abs_error * mean_abs_error) / (epsilon * epsilon) + log(1.0 / (2.0 * M_PI * epsilon * epsilon))) / 2.0;
            double new_tempered_logp = pl + ll * beta;
            double delta_tempered_logp = new_tempered_logp - tempered_logp[draw];

            double r = log(uni_dist(rng));
            // std::cout << r << "  " << delta_tempered_logp << "\n";
            if (r < delta_tempered_logp){
                q_old = q_new;
                accepted += 1;
                old_prior = pl;
                old_likelihood = ll;
                old_tempered_logp = new_tempered_logp;
            }

            if (DEBUG){
                char sacc = '-';
                std::cout << std::setw(8) << std::fixed << std::setprecision(3) << mean_abs_error << " | ";
                std::cout << std::setw(8) << std::fixed << std::setprecision(3) << pl << " + ";
                std::cout << std::setw(12) << std::fixed << std::setprecision(2) << ll << " * ";
                std::cout << std::setw(8) << std::fixed << std::setprecision(6) << beta << " = ";
                std::cout << std::setw(8) << std::fixed << std::setprecision(3) << new_tempered_logp << " - ";
                std::cout << std::setw(8) << std::fixed << std::setprecision(3) << tempered_logp[draw];
                std::cout << std::setw(10) << std::fixed << std::setprecision(3) << " -> "  << delta_tempered_logp << " ?? ";
                std::cout << std::setw(8) << std::fixed << std::setprecision(3) << r << " ";
                if (r < delta_tempered_logp)
                    sacc = '+';
                std::cout << sacc << "\n";
            }

        }
        prior_likelihoods(draw) = old_prior;
        likelihoods(draw) = old_likelihood;
        tempered_logp(draw) = old_tempered_logp;
        acc_per_chain(draw) = double(accepted) / double(n_steps);
        posteriors.row(draw) = q_old;

        for (int i = 0; i < nparams; ++i)
        {
            posteriors(draw, i) = q_old[i];
        }
    }

}

void AbcSmcFit::Fit(SpectrumModel spectrumModel) {

    //abc-smc fitting
    int nparams = 6;
    int draws = 1000;
    double epsilon = 0.01;
    int stage = 0;
    double beta = 0.0;
    int marginal_likelihood = 1;
    double threshold=0.5;
    double acc_rate = 1.0;
    int n_steps = 25;
    double p_acc_rate=0.99;
    bool tune_steps= true;
    int max_steps = n_steps;
    int dimension = nparams;
    double factor = (2.38 * 2.38) / nparams;
    int proposed = draws * n_steps;

    Eigen::ArrayX<double> acc_per_chain(draws);
    Eigen::ArrayX<double> scalings(draws);
    Eigen::ArrayX<double> prior_likelihoods(draws);
    Eigen::ArrayXX<double> posteriors(draws , nparams);
    /*   x0 fwhm0 int0 x1 fwhm1 int2
       [ 0    0    0   0    0    0 ]
       [ 0    0    0   0    0    0 ]
       [ 0    0    0   0    0    0 ]
       [ 0    0    0   0    0    0 ]
       [ 0    0    0   0    0    0 ]

       Posteriors vector:
       columns = nparams
       row = drwas
     * */
    SpectrumModel simulatedSpectrumModel = SpectrumModel(spectrumModel.npix);
    std::vector<NormalDistribution> normalDistribution = spectrumModel.NormalDistributions;

    if (factor < 1){
        scalings = factor;
    } else{
        scalings = 1;
    }

    /*' ___________________________________________'
    ' |                                         |'
    ' | Kezdeti populacio meghatarozasa         |'
    ' |_________________________________________|'
    '--------------------------------------------------------------------------------------------------------'*/
    for (int i = 0; i < normalDistribution.size(); ++i) {
        Eigen::ArrayX<double> tmp = normalDistribution[i].Sample(draws);
        posteriors.col(i) = tmp;
    }

    Eigen::ArrayXX<double> init_population = posteriors;
    Eigen::ArrayX<double> likelihoods(draws);

    std::cout << "Starting population statistics: \n    ";

    Eigen::MatrixXd copyPosteriors = posteriors;
    populationStatistics(copyPosteriors);
    std::cout << "\n";

    for (int i = 0; i < draws; ++i) {
        //TODO
        for (int j = 0; j < posteriors.row(i).size(); ++j) {
            prior_likelihoods(i) += normalDistribution[j].LogP(init_population(i, j));
        }

        std::map<std::string, Eigen::ArrayX<double>> simSpectrum;
        Eigen::ArrayX<double> peak1FromInitPopulation(3);
        peak1FromInitPopulation << init_population(i, 0), init_population(i, 1), init_population(i, 2);

        Eigen::ArrayX<double> peak2FromInitPopulation(3);
        peak2FromInitPopulation << init_population(i, 3), init_population(i, 4), init_population(i, 5);

        simSpectrum.insert(std::pair<std::string, Eigen::ArrayX<double>>("Gaussian", peak1FromInitPopulation));
        simSpectrum.insert(std::pair<std::string, Eigen::ArrayX<double>>("Lorentzian", peak2FromInitPopulation));
        simulatedSpectrumModel.calculate(simSpectrum, false);

        Eigen::ArrayX<double> spectrumDiffs(spectrumModel.npix);
        spectrumDiffs = (spectrumModel.spectrum - simulatedSpectrumModel.spectrum).cwiseAbs();

        double mean_abs_error = spectrumDiffs.mean();

        likelihoods(i) = (-(mean_abs_error * mean_abs_error) / (epsilon * epsilon) + log(1.0 / (2.0 * M_PI * epsilon * epsilon))) / 2.0;
    }


    while (beta < 1){
        std::cout << "Stage: " << std::setw(3) << stage << "  ";
        std::cout << "Beta: " << std::fixed << std::setprecision(8) << beta << "  ";
        std::cout << "Steps: " << std::setw(2) << n_steps << "  ";
        std::cout << "Acceptance rate: " << std::fixed << std::setprecision(8) << acc_rate << std::endl;

        double lowBeta = beta;
        double oldBeta = beta;
        double upBeta = 2.0;
        double newBeta;

        int rN = int(round(int(likelihoods.size()) * threshold));

        Eigen::ArrayX<double> ll_diffs(likelihoods.size());
        double max = likelihoods.maxCoeff();
        ll_diffs = likelihoods - max;

//    ' ___________________________________________'
//    ' |                                         |'
//    ' | Fontossagi sulyok es beta meghatarozasa |'
//    ' |_________________________________________|'
#pragma region importance_weights
        std::cout << "    Stage " << stage << " - Importance weights\n";
        Eigen::ArrayX<double> weights_un(ll_diffs.size());
        Eigen::ArrayX<double> weights(ll_diffs.size());

        while ((upBeta - lowBeta) > 1e-6 ){
            double sum_of_weights_un = 0;
            double sum_of_square_weights = 0;
            newBeta = (lowBeta + upBeta) / 2.0;

            weights_un = exp(newBeta - oldBeta) * ll_diffs;
            sum_of_weights_un = weights_un.sum();

            weights = weights_un / sum_of_weights_un;
            sum_of_square_weights = sum_of_square_weights + weights.square().sum();

            int ESS = int(round(1.0 / sum_of_square_weights));
            if (ESS == rN){
                break;
            }
            else if(ESS < rN){
                upBeta = newBeta;
            }
            else{
                lowBeta = newBeta;
            }
        }

        if (newBeta >= 1){
            double sum_of_weights_un = 0;
            newBeta = 1;
            weights_un = exp(newBeta - oldBeta) * ll_diffs;
            sum_of_weights_un = weights_un.sum();
            weights = weights_un / sum_of_weights_un;
        }

        beta = newBeta;
        std::cout << "    Stage " << stage << " - new beta = " << beta << "\n";
#pragma endregion importance_weights
//    ' ______________________________________________'
//    ' |                                            |'
//    ' | Mintavetelezes a fontossagi sulyok alapjan |'
//    ' |____________________________________________|'

#pragma region resampling
        std::cout << "    Stage " << stage << " - Resampling\n";
        //TODO
        Eigen::ArrayX<int> resamplingIndexes = randomWeightedIndices(draws, weights);

        posteriors = resampling(posteriors, resamplingIndexes);
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
//    ' ______________________________________________'
//    ' |                                            |'
//    ' | Algoritmus hangolasa                       |'
//    ' |____________________________________________|'

#pragma region tuning
        std::cout << "    Stage " << stage << " - Tuning\n";
        if (stage > 0){
            double ave_scalings = exp(log(scalings.mean() + acc_per_chain.mean() - 0.234));

            std::cout << "      <scalings> = " << scalings.mean() << "\n";
            std::cout << "      <acc_per_chain> = " << acc_per_chain.mean() << "\n";
            std::cout << "      ave_scaling = " << ave_scalings << "\n";

            scalings = 0.5 * (ave_scalings + exp(scalings.log() + (acc_per_chain - 0.234)));
            std::cout << "      <scalings> = " << scalings.mean() << "\n";

            if (tune_steps){
                acc_rate = std::max(1.0 / proposed, acc_rate);
                int t = (int) round((log(1.0 - p_acc_rate) / log(1.0 - acc_rate)));
                n_steps = std::min(max_steps, std::max(2, t));
            }
            proposed = draws * n_steps;
            std::cout << "      acc_rate = " << acc_rate << "\n";
        }

#pragma endregion tuning


//    ' ______________________________________________'
//    ' |                                            |'
//    ' | A monte carlo lancok futtatasa             |'
//    ' |____________________________________________|'
#pragma region MCMC
        std::cout << "    Stage " << stage << " - MCMC chains\n";

        Eigen::MatrixXd copyPosteriorsInMC = posteriors;
        populationStatistics(copyPosteriorsInMC);

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
                      simulatedSpectrumModel,
                      spectrumModel);

        Eigen::MatrixXd copyPosteriorsInLoop = posteriors;
        populationStatistics(copyPosteriorsInLoop);
        acc_rate = acc_per_chain.mean();

        #pragma endregion MCMC

        stage += 1;
        if (stage > 20) {break;}
    }
}