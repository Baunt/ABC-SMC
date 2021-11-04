#define _USE_MATH_DEFINES
#include <iomanip>
#include <random>
#include "util.h"
#include "peak_model.h"
#include "../include/third-party-library/matplotlib-cpp/matplotlibcpp.h"
#include "../include/third-party-library/Eigen/Core"
#include "probability_distribution.h"
#include "../include/third-party-library/Eigen/src/Cholesky/LLT.h"
#include <algorithm>
#include <math.h> 

int main() {

    //peak 1
    double real_x0 = 0.67;
    double real_fwhm0 = 0.11;
    double real_intensity0 = 2.3;
    //peak 2
    double real_x1 = 0.51;
    double real_fwhm1 = 0.20;
    double real_intensity1 = 0.5;

    std::vector<double> real_y(256);
    Eigen::VectorXd real_y_eigen(256);
    std::vector<double> y_sim(256);

    //simulated spectrum
    int npix = 256;
    std::vector<double> x = linspace(0, 1, npix);

    std::vector<double> noise = getDistribution(0.0, 1.0, npix);

    PeakModel gaussianPeakModel = PeakModel(x, real_x0, real_fwhm0, real_intensity0, npix);
    std::vector<double> gaussianModel = gaussianPeakModel.Gaussian();

    PeakModel lorentzianPeakModel = PeakModel(x, real_x1, real_fwhm1, real_intensity1, npix);
    std::vector<double> lorentzianModel = gaussianPeakModel.Lorenzt();

    for (int i = 0; i < gaussianModel.capacity(); ++i) {
        real_y[i] = gaussianModel[i] + lorentzianModel[i];
        real_y_eigen(i) = gaussianModel[i] + lorentzianModel[i];
    }

    for (int i = 0; i < real_y.capacity(); ++i) {
        real_y[i] += noise[i];
        real_y_eigen(i) += noise[i];
    }

    //Initial tips
    NormalDistribution x0 = NormalDistribution("x0", "normal", 0.63, 0.15);
    NormalDistribution fwhm0 = NormalDistribution("fwhm0", "normal", 0.09, 0.04);
    NormalDistribution int0 = NormalDistribution("int0", "normal", 2.65, 0.5);
    NormalDistribution x1 = NormalDistribution("x1", "normal", 0.45, 0.30);
    NormalDistribution fwhm1 = NormalDistribution("fwhm1", "normal", 0.25, 0.07);
    NormalDistribution int1 = NormalDistribution("int1", "normal", 0.35, 0.15);

    std::vector<NormalDistribution> normalDistribution = {x0, fwhm0, int0, x1, fwhm1, int1};
    //abc-smc fitting
    int nparams = 6;
    int draws = 1000;
    double epsilon = 0.01;

    std::vector<double> priors(draws);
    std::vector<std::vector<double>> posteriors( draws , std::vector<double> (nparams, 0));

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

    std::vector<double>acc_per_chain(draws, 1);
    std::vector<double> scalings(draws);
    for (int i = 0; i < draws; ++i) {
        double factor = (2.38 * 2.38) / nparams;
        if (factor < 1){
            scalings[i] = factor;
        } else{
            scalings[i] = 1;
        }
    }

    int proposed = draws * n_steps;

    /*' ___________________________________________'
    ' |                                         |'
    ' | Kezdeti populacio meghatarozasa         |'
    ' |_________________________________________|'
    '--------------------------------------------------------------------------------------------------------'
    'A modell parameterekhez tartozo eloszlasokbol (ebben a kodban csak NormalDistribution opcio van) veszunk'
    ' "draws" darab mintat '*/

    std::vector<std::vector<double>> transposedPosterior = transpose(posteriors);

    // for (int i = 0; i < normalDistribution.capacity(); ++i) {
    //     for (int j = 0; j < normalDistribution[i].Sample(draws).capacity(); ++j) {
    //         transposedPosterior[i][j] = normalDistribution[i].Sample(draws)[j];
    //     }
    // }
    

    
    for (int i = 0; i < normalDistribution.capacity(); ++i) {
        auto tmp = normalDistribution[i].Sample(draws);
        for (int j = 0; j < draws; ++j) {
            transposedPosterior[i][j] = tmp[j];
        }
        std::cout << i << "\n";
        histogram(tmp, 20);
    }


    std::vector<std::vector<double>> init_population = transpose(transposedPosterior);
    posteriors = init_population;
    std::vector<double> likelihoods(draws, 0);

    //matplotlibcpp::plot(real_y);
//    'A kezdeti populacionak kiszamoljuk a logp es likelihood_logp ertekeit.'
    for (int i = 0; i < draws; ++i) {
        for (int j = 0; j < posteriors[i].capacity(); ++j) {
            priors[i] += normalDistribution[j].LogP(init_population[i][j]);
        }

        PeakModel gaussModel = PeakModel(x, init_population[i][0], init_population[i][1], init_population[i][2], npix);
        std::vector<double> gauss = gaussModel.Gaussian();
        PeakModel lorentzModel = PeakModel(x, init_population[i][3], init_population[i][4], init_population[i][5], npix);
        std::vector<double> lorentz = lorentzModel.Lorenzt();
        for (int i = 0; i < gauss.capacity(); ++i) {
            y_sim[i] = gauss[i] + lorentz[i];
        }

        //matplotlibcpp::plot(y_sim);

        std::vector<double> spectrumDiffs(256);
        for (int j = 0; j < real_y.capacity(); ++j) {
            spectrumDiffs[j] = abs(real_y[j] - y_sim[j]);
        }

        double mean_abs_error = arithmetic_mean(spectrumDiffs);

        likelihoods[i] = (-(mean_abs_error * mean_abs_error) / epsilon * epsilon + log(1 / (2 * M_PI * epsilon * epsilon))) / 2.0;
    }

    //    'Itt kezdodik az iteracio'
    while (beta < 1){
        std::cout << std::setprecision(5) << "Stage: " << stage << std::endl;
        std::cout << std::setprecision(5) << "Beta: " << beta << std::endl;
        std::cout << std::setprecision(5) << "Steps: " << n_steps << std::endl;
        std::cout << std::setprecision(5) << "Accelerate: " << acc_rate << std::endl;

        double lowBeta = beta;
        double oldBeta = beta;
        double upBeta = 2.0;
        double newBeta;

        double rN = likelihoods.capacity() * threshold;

        std::vector<double> ll_diffs(likelihoods.capacity());
        double max = *std::max_element(likelihoods.begin(), likelihoods.end());
        for (int i = 0; i < likelihoods.capacity(); ++i) {
            ll_diffs[i] = likelihoods[i] - max;
        }

//    ' ___________________________________________'
//    ' |                                         |'
//    ' | Fontossagi sulyok es beta meghatarozasa |'
//    ' |_________________________________________|'

        std::vector<double> weights_un(ll_diffs.capacity());
        std::vector<double> weights(ll_diffs.capacity());

        while ((upBeta - lowBeta) > 1e-6 ){
            double sum_of_weights_un = 0;
            double sum_of_weights = 0;
            newBeta = (lowBeta + upBeta) / 2.0;

            for (int i = 0; i < ll_diffs.capacity(); ++i) {
                weights_un[i] = exp((newBeta - oldBeta) * ll_diffs[i]);
                sum_of_weights_un += weights_un[i];
            }

            for (int i = 0; i < weights_un.capacity(); ++i) {
                weights[i] = weights_un[i] / sum_of_weights_un;
                sum_of_weights += weights[i] * weights[i];
            }

            double ESS = 1 / sum_of_weights;
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

            for (int i = 0; i < ll_diffs.capacity(); ++i) {
                weights_un[i] = exp((newBeta - oldBeta) * ll_diffs[i]);
                sum_of_weights_un += weights_un[i];
            }

            for (int i = 0; i < weights_un.capacity(); ++i) {
                weights[i] = weights_un[i] / sum_of_weights_un;
            }
        }

        //marginal_likelihood????

        beta = newBeta;

//    ' ______________________________________________'
//    ' |                                            |'
//    ' | Mintavetelezes a fontossagi sulyok alapjan |'
//    ' |____________________________________________|'

        
        // auto resamplingIndexes = randomWeightedIndices(draws, weights);

        // testing resampling


        /* drawing 10000 numbers from a set of 100 weighted indices */
        int nweights = 100;
        int ndraws = 10000;

        // creating random weights
        std::random_device randomDevice;
        std::mt19937 randomGen(randomDevice());
        std::normal_distribution<double> distribution(0.0, 1.0);
        std::vector<double> randomWeights(nweights, 0.0);
        double norm = 0.0;
        for (int i = 0; i < nweights; ++i)
        {
            double r = abs(distribution(randomGen)); // half normal distribution (most probabilities are small, there are a few larger ones)
            randomWeights[i] = r;
            norm += r;
        }
        
        // normalizing weights
        norm = 1.0 / norm;
        for (int i = 0; i < nweights; ++i)
        {
            randomWeights[i] *= norm;
        }


        // drawing random indexes (~ randomchoice of numpy)
        auto resamplingIndexes = randomWeightedIndices(ndraws, randomWeights);


        // finally, check the results
        std::vector<int> numberoftimes(nweights, 0); // this will store how many times each index was drawn

        for (int i = 0; i < ndraws; ++i)
        {
            numberoftimes[resamplingIndexes[i]] += 1;
        }
        std::cout << "Testing the resampling algorithm:\n Index    Probability    Number of draws\n";


        for (int i = 0; i < nweights; ++i)
        {
            std::cout << std::setw(4) << i << std::fixed << std::setw(16) << std::setprecision(8) << randomWeights[i] << std::setw(14) << numberoftimes[i] << "\n";
        }

        posteriors = resampling(posteriors, resamplingIndexes);
        priors = resampling(priors, resamplingIndexes);
        likelihoods = resampling(likelihoods, resamplingIndexes);
        scalings = resampling(scalings, resamplingIndexes);
        acc_per_chain = resampling(acc_per_chain, resamplingIndexes);

        std::vector<double>tempered_logp(draws);
        for (int i = 0; i < draws; ++i) {
            tempered_logp[i] = priors[i] + likelihoods[i] * beta;
        }


        Eigen::MatrixXd post(draws, nparams);
        for (int i = 0; i < draws; ++i) {
            for (int j = 0; j < nparams; ++j) {
                post(i, j) = posteriors[i][j];
            }
        }

        Eigen::MatrixXd centered = post.rowwise() - post.colwise().mean();
        Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(post.rows() - 1);
        Eigen::LLT<Eigen::MatrixXd> llt(cov);

//        Eigen::MatrixXd L = lltOfCovMatrix.matrixL();
//        std::cout << "The Cholesky factor L is" << std::endl << L << std::endl;
//        std::cout << "To check this, let us compute L * L.transpose()" << std::endl;
//        std::cout << L * L.transpose() << std::endl;
//        std::cout << "This should equal the matrix A" << std::endl << cov << std::endl;

//    ' ______________________________________________'
//    ' |                                            |'
//    ' | Algoritmus hangolasa                       |'
//    ' |____________________________________________|'

        if (stage > 0){
            auto ave_scalings = exp(log(arithmetic_mean(scalings) + arithmetic_mean(acc_per_chain) - 0.234));

            for (int i = 0; i < scalings.capacity(); ++i) {
                scalings[i] = 0.5 * (ave_scalings + exp(log(scalings[i])) + (acc_per_chain[i] - 0.234));
            }

            if (tune_steps){
                acc_rate = std::max(1.0 / proposed, acc_rate);
                int t = (int) (log(1-p_acc_rate) / log(1- acc_rate));
                n_steps = std::min(max_steps, std::max(2, t));
            }
            proposed = draws * n_steps;
        }

//    ' ______________________________________________'
//    ' |                                            |'
//    ' | A monte carlo lancok futtatasa             |'
//    ' |____________________________________________|'

        Eigen::MatrixXd new_posteriors(post);
        Eigen::VectorXd new_priors(priors.capacity());
        Eigen::VectorXd new_likelihoods(likelihoods.capacity());
        Eigen::VectorXd new_acc_per_chain(acc_per_chain.capacity());
        double new_acc_rate = 0.0;

        for (int draw = 0; draw < draws; ++draw) {

            Eigen::MatrixXd randomsForDeltas = Eigen::MatrixXd::Random(n_steps, 6);
            Eigen::MatrixXd deltas = randomsForDeltas * llt.matrixLLT();
            deltas *= scalings[draw];

            int accepted = 0;

            for (int i = 0; i < n_steps; ++i) {
                Eigen::VectorXd delta = deltas.row(n_steps - 1);
                std::cout << delta << std::endl;
                Eigen::VectorXd q_old(posteriors[draw].capacity());
                for (int j = 0; j < posteriors[draw].capacity(); ++j) {
                    q_old(j) = posteriors[draw][j];
                }

                Eigen::VectorXd q_new = q_old + delta;
                std::cout << q_new << std::endl;

                double pl = 0.0;
                for (int dist = 0; dist < normalDistribution.capacity(); ++dist) {
                    pl += normalDistribution[dist].LogP(q_new[dist]);
                }

                PeakModel gaussModel = PeakModel(x, q_new[0], q_new[1], q_new[2], npix);
                std::vector<double> gauss = gaussModel.Gaussian();
                PeakModel lorentzModel = PeakModel(x, q_new[3], q_new[4], q_new[5], npix);
                std::vector<double> lorentz = lorentzModel.Lorenzt();
                for (int peak = 0; peak < gauss.capacity(); ++peak) {
                    y_sim[i] = gauss[peak] + lorentz[peak];
                }
                std::vector<double> spectrumDiffs(256);
                for (int j = 0; j < real_y.capacity(); ++j) {
                    spectrumDiffs[j] = abs(real_y[j] - y_sim[j]);
                }

                double mean_abs_error = arithmetic_mean(spectrumDiffs);
                double ll = (-(mean_abs_error * mean_abs_error) / epsilon * epsilon + log(1 / (2 * M_PI * epsilon * epsilon))) / 2.0;
                double new_tempered_logp = pl + ll * beta;
                double sub_new_old_tempered_logp = new_tempered_logp - tempered_logp[draw];

                double old_prior = 0.0;
                double old_likelihood = 0.0;
                double old_tempered_logp = 0.0;

                if (getUniformRandomNumber() < sub_new_old_tempered_logp){
                    q_old = q_new;
                    accepted += 1;
                    old_prior = pl;
                    old_likelihood = ll;
                    old_tempered_logp = new_tempered_logp;
                }
                
            }
        }

//    results = [metrop_kernel(x, posterior[draw], tempered_logp[draw], priors[draw], likelihoods[draw], draw, *parameters) for draw in range(draws)]
//
//    posterior, acc_list, priors, likelihoods = zip(*results)
//    posterior = np.array(posterior)
//    priors = np.array(priors)
//    likelihoods = np.array(likelihoods)
//    acc_per_chain = np.array(acc_list)
//    acc_rate = np.mean(acc_list)

        stage += 1;
    }



    // matplotlibcpp::show();

    return 0;
}
