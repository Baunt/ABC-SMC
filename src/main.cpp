#define _USE_MATH_DEFINES
#include <iomanip>
#include <random>
#include "util.h"
#include "peak_model.h"
// #include "../include/third-party-library/matplotlib-cpp/matplotlibcpp.h"
#include "../include/third-party-library/Eigen/Core"
#include "probability_distribution.h"
#include "../include/third-party-library/Eigen/src/Cholesky/LLT.h"
#include <algorithm>
#include <cmath>
// #include <math.h> 

int main(int argc, char** argv) {
    bool DEBUG = false;
    if (argc > 1)
        DEBUG = true;
    int npix = 256;

    //peak 1
    double real_x0 = 0.67;
    double real_fwhm0 = 0.11;
    double real_intensity0 = 2.3;
    //peak 2
    double real_x1 = 0.51;
    double real_fwhm1 = 0.20;
    double real_intensity1 = 0.5;

    std::vector<double> real_y(npix);
    Eigen::VectorXd real_y_eigen(npix);
    std::vector<double> y_sim(npix);

    //simulated spectrum
    
    std::vector<double> x = linspace(0, 1, npix);

    std::vector<double> noise = getDistribution(0.0, 1.0, npix);

    PeakModel gaussianPeakModel = PeakModel(x, real_x0, real_fwhm0, real_intensity0, npix);
    std::vector<double> gaussianModel = gaussianPeakModel.Gaussian();
    PeakModel lorentzianPeakModel = PeakModel(x, real_x1, real_fwhm1, real_intensity1, npix);
    std::vector<double> lorentzianModel = lorentzianPeakModel.Lorenzt();

    for (int i = 0; i < gaussianModel.capacity(); ++i) {
        real_y[i] = gaussianModel[i] + lorentzianModel[i];
        real_y_eigen(i) = gaussianModel[i] + lorentzianModel[i];
    }

    if (DEBUG){
        for (int i = 0; i < gaussianModel.capacity(); ++i) {
            std::cout << x[i] << "  " << gaussianModel[i] << "  " << lorentzianModel[i] << "  " << noise[i] << "\n";
        }        
        int asdff  = 0;
        std::cin >> asdff; }

    double realerror = 0.0;
    for (int i = 0; i < real_y.capacity(); ++i) {
        real_y[i] += noise[i];
        real_y_eigen(i) += noise[i];
        realerror += std::abs(noise[i]);
    }
    realerror /= double(npix);
    std::cout << "Mean abs. error of model: " << realerror << "\n";
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

    std::vector<double> prior_likelihoods(draws);
    std::vector<std::vector<double>> posteriors(draws , std::vector<double> (nparams, 0));

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
    
    for (int i = 0; i < normalDistribution.capacity(); ++i) {
        auto tmp = normalDistribution[i].Sample(draws);
        for (int j = 0; j < draws; ++j) {
            posteriors[j][i] = tmp[j];
        }
        // histogram(tmp, 20);
    }

    std::vector<std::vector<double>> init_population = posteriors;
    posteriors = init_population;
    std::vector<double> likelihoods(draws, 0);

    std::cout << "Starting population statistics: \n    ";
    populationstatistics(posteriors);
    std::cout << "\n";

    //matplotlibcpp::plot(real_y);
//    'A kezdeti populacionak kiszamoljuk a logp es likelihood_logp ertekeit.'
    for (int i = 0; i < draws; ++i) {
        for (int j = 0; j < posteriors[i].capacity(); ++j) {
            prior_likelihoods[i] += normalDistribution[j].LogP(init_population[i][j]);
        }

        PeakModel gaussModel = PeakModel(x, init_population[i][0], init_population[i][1], init_population[i][2], npix);
        std::vector<double> gauss = gaussModel.Gaussian();
        PeakModel lorentzModel = PeakModel(x, init_population[i][3], init_population[i][4], init_population[i][5], npix);
        std::vector<double> lorentz = lorentzModel.Lorenzt();
        for (int i = 0; i < gauss.capacity(); ++i) {
            y_sim[i] = gauss[i] + lorentz[i];
        }

        //matplotlibcpp::plot(y_sim);

        std::vector<double> spectrumDiffs(npix);
        for (int j = 0; j < real_y.capacity(); ++j) {
            spectrumDiffs[j] = std::abs(real_y[j] - y_sim[j]);
        }

        double mean_abs_error = arithmetic_mean(spectrumDiffs);

        likelihoods[i] = (-(mean_abs_error * mean_abs_error) / (epsilon * epsilon) + log(1.0 / (2.0 * M_PI * epsilon * epsilon))) / 2.0;
    }

    //    'Itt kezdodik az iteracio'

    // this is a uniform random generator used for the random steps:
    std::random_device randomDevice;
    std::mt19937 randomGen(randomDevice());
    std::normal_distribution<double> std_dist(0.0, 1.0); // normal distribution with mu=0, sigma=1
    std::uniform_real_distribution<double> uni_dist(0.0, 1.0); // uniform ditribution between: 0 <= r < 1

    if (DEBUG){
        int asd  = 0;
        std::cin >> asd; }

    while (beta < 1){
        std::cout << "Stage: " << std::setw(3) << stage << "  ";
        std::cout << "Beta: " << std::fixed << std::setprecision(8) << beta << "  ";
        std::cout << "Steps: " << std::setw(2) << n_steps << "  ";
        std::cout << "Acceptance rate: " << std::fixed << std::setprecision(8) << acc_rate << std::endl;
        
        double lowBeta = beta;
        double oldBeta = beta;
        double upBeta = 2.0;
        double newBeta;

        int rN = int(round(likelihoods.capacity() * threshold));

        std::vector<double> ll_diffs(likelihoods.capacity());
        double max = *std::max_element(likelihoods.begin(), likelihoods.end());
        for (int i = 0; i < likelihoods.capacity(); ++i) {
            ll_diffs[i] = likelihoods[i] - max;
        }

//    ' ___________________________________________'
//    ' |                                         |'
//    ' | Fontossagi sulyok es beta meghatarozasa |'
//    ' |_________________________________________|'
        #pragma region importance_weights        
        std::cout << "    Stage " << stage << " - Importance weights\n";
        std::vector<double> weights_un(ll_diffs.capacity());
        std::vector<double> weights(ll_diffs.capacity());

        while ((upBeta - lowBeta) > 1e-6 ){
            double sum_of_weights_un = 0;
            double sum_of_square_weights = 0;
            newBeta = (lowBeta + upBeta) / 2.0;

            for (int i = 0; i < ll_diffs.capacity(); ++i) {
                weights_un[i] = exp((newBeta - oldBeta) * ll_diffs[i]);
                sum_of_weights_un += weights_un[i];
            }

            for (int i = 0; i < weights_un.capacity(); ++i) {
                weights[i] = weights_un[i] / sum_of_weights_un;
                sum_of_square_weights += weights[i] * weights[i];
            }

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
        std::cout << "    Stage " << stage << " - new beta = " << beta << "\n";        
        #pragma endregion importance_weights
//    ' ______________________________________________'
//    ' |                                            |'
//    ' | Mintavetelezes a fontossagi sulyok alapjan |'
//    ' |____________________________________________|'        

        #pragma region resampling
        std::cout << "    Stage " << stage << " - Resampling\n";
        auto resamplingIndexes = randomWeightedIndices(draws, weights);
        posteriors = resampling(posteriors, resamplingIndexes);
        prior_likelihoods = resampling(prior_likelihoods, resamplingIndexes);
        likelihoods = resampling(likelihoods, resamplingIndexes);
        scalings = resampling(scalings, resamplingIndexes);
        acc_per_chain = resampling(acc_per_chain, resamplingIndexes);

        std::vector<double>tempered_logp(draws);
        for (int i = 0; i < draws; ++i) {
            tempered_logp[i] = prior_likelihoods[i] + likelihoods[i] * beta;
        }


        Eigen::MatrixXd posteriors_matrix(draws, nparams);
        for (int i = 0; i < draws; ++i) {
            for (int j = 0; j < nparams; ++j) {
                posteriors_matrix(i, j) = posteriors[i][j];
            }
        }

        Eigen::MatrixXd centered = posteriors_matrix.rowwise() - posteriors_matrix.colwise().mean();
        Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(posteriors_matrix.rows() - 1);
        Eigen::LLT<Eigen::MatrixXd> llt(cov);

//        Eigen::MatrixXd L = lltOfCovMatrix.matrixL();
//        std::cout << "The Cholesky factor L is" << std::endl << L << std::endl;
//        std::cout << "To check this, let us compute L * L.transpose()" << std::endl;
//        std::cout << L * L.transpose() << std::endl;
//        std::cout << "This should equal the matrix A" << std::endl << cov << std::endl;

        #pragma endregion resampling             
//    ' ______________________________________________'
//    ' |                                            |'
//    ' | Algoritmus hangolasa                       |'
//    ' |____________________________________________|'

        #pragma region tuning
        std::cout << "    Stage " << stage << " - Tuning\n";   
        if (stage > 0){
            auto ave_scalings = std::exp(std::log(arithmetic_mean(scalings)) + (arithmetic_mean(acc_per_chain) - 0.234));


            std::cout << "      <scalings> = " << arithmetic_mean(scalings) << "\n";
            std::cout << "      <acc_per_chain> = " << arithmetic_mean(acc_per_chain) << "\n";
            std::cout << "      ave_scaling = " << arithmetic_mean(scalings) << "\n";


            for (int i = 0; i < scalings.capacity(); ++i) {
                scalings[i] = 0.5 * (ave_scalings + std::exp(std::log(scalings[i]) + (acc_per_chain[i] - 0.234)));
            }

            std::cout << "      <scalings> = " << arithmetic_mean(scalings) << "\n";


            if (tune_steps){
                acc_rate = std::max(1.0 / proposed, acc_rate);
                int t = (int) round((std::log(1.0 - p_acc_rate) / std::log(1.0 - acc_rate)));
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
        
        
        populationstatistics(posteriors_matrix);
        double new_acc_rate = 0.0;

        for (int draw = 0; draw < draws; ++draw) {
            if (DEBUG)
            {
                std::cout << "  Sample " << draw << "\n";
                std::cout << posteriors[draw][0] << " "<< posteriors[draw][1] << " "<< posteriors[draw][2] << " "<< posteriors[draw][3] << " "<< posteriors[draw][4] << " "<< posteriors[draw][5] <<"\n";
                std::cout << "    err        pl            ll          beta     logp_0    logp_1          d_logp         r\n";
            }


            //Eigen::MatrixXd randomsForDeltas = Eigen::MatrixXd::Random(n_steps, 6);
            auto rand_fn = [&](){return std_dist(randomGen);};                                 // lambda that generates random numbers
            Eigen::MatrixXd random_matrix = Eigen::MatrixXd::NullaryExpr(n_steps, 6, rand_fn); // matrix expression
            Eigen::MatrixXd randomsForDeltas = random_matrix;                                  // saving the elements to a matrix
            Eigen::MatrixXd deltas = randomsForDeltas * llt.matrixLLT();            
            deltas *= scalings[draw];


            int accepted = 0;

            double old_likelihood = likelihoods[draw];
            double old_prior = prior_likelihoods[draw];
            double old_tempered_logp = tempered_logp[draw];
            Eigen::VectorXd q_old = posteriors_matrix.row(draw);

            for (int i = 0; i < n_steps; ++i) {
                Eigen::VectorXd delta = deltas.row(i);
                Eigen::VectorXd q_new = q_old + delta;                
                // std::cout << delta[0] << " " << delta[1] << " "<< delta[2] << "" << delta[3] << " "<< delta[4] << " "<< delta[0] << "\n";


                // Eigen::VectorXd q_old(posteriors[draw].capacity());
                // for (int j = 0; j < posteriors[draw].capacity(); ++j) {
                //     q_old(j) = posteriors[draw][j];
                // }

                // std::cout << "\n\n" << q_new[0] << " " << q_new[1] << " " << q_new[2] << " " << q_new[3] << " " << q_new[4] << " " << q_new[5] << "\n";
                // PeakModel gaussModel = PeakModel(x, q_new[0], q_new[1], q_new[2], npix);
                // std::vector<double> gauss = gaussModel.Gaussian();
                // PeakModel lorentzModel = PeakModel(x, q_new[3], q_new[4], q_new[5], npix);
                // std::vector<double> lorentz = lorentzModel.Lorenzt();
                // for (int peak = 0; peak < gauss.capacity(); ++peak) {
                //     // std::cout << x[peak]<< "    " << gauss[peak] << "    " << lorentz[peak] << "\n";
                //     y_sim[i] = gauss[peak] + lorentz[peak];
                // }

                y_sim = staticpeakmodel(x, q_new);
                std::vector<double> spectrumDiffs(npix);
                for (int j = 0; j < real_y.capacity(); ++j) {
                    spectrumDiffs[j] = std::abs(real_y[j] - y_sim[j]);
                }
                double pl = 0.0;
                for (int dist = 0; dist < normalDistribution.capacity(); ++dist) {
                    pl += normalDistribution[dist].LogP(q_new[dist]);
                }
                double mean_abs_error = arithmetic_mean(spectrumDiffs);
                double ll = (-(mean_abs_error * mean_abs_error) / (epsilon * epsilon) + log(1.0 / (2.0 * M_PI * epsilon * epsilon))) / 2.0;
                double new_tempered_logp = pl + ll * beta;
                double delta_tempered_logp = new_tempered_logp - tempered_logp[draw];


                // double old_prior = 0.0;
                // double old_likelihood = 0.0;
                // double old_tempered_logp = 0.0;

                /*
                q_old, accept = metrop_select(new_tempered_logp - old_tempered_logp, q_new, q_old)
                if accept:
                    accepted += 1
                    old_prior = pl
                    old_likelihood = ll
                    old_tempered_logp = new_tempered_logp
                */
                
                double r = log(uni_dist(randomGen));
                // std::cout << r << "  " << delta_tempered_logp << "\n";
                if (r < delta_tempered_logp){
                    q_old = q_new;
                    accepted += 1;
                    old_prior = pl;
                    old_likelihood = ll;
                    old_tempered_logp = new_tempered_logp;
                }
                // new_likelihoods(draw) = old_likelihood;
                
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
            prior_likelihoods[draw] = old_prior;
            likelihoods[draw] = old_likelihood;
            tempered_logp[draw] = old_tempered_logp;
            acc_per_chain[draw] = double(accepted) / double(n_steps);
            posteriors_matrix.row(draw) = q_old;        
            for (int i = 0; i<nparams; ++i)
            {
                posteriors[draw][i] = q_old[i];
            }            
        }
        
        populationstatistics(posteriors_matrix);
        acc_rate = arithmetic_mean(acc_per_chain);

        #pragma endregion MCMC

        stage += 1;
        if (stage > 20) {break;}
    }



    // matplotlibcpp::show();

    return 0;
}
