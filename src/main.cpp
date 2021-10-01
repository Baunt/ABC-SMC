#define _USE_MATH_DEFINES
#include <iomanip>
#include "util.h"
#include "peak_model.h"
#include "../include/third-party-library/matplotlib-cpp/matplotlibcpp.h"
#include "probability_distribution.h"

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
    }

    for (int i = 0; i < real_y.capacity(); ++i) {
        real_y[i] += noise[i];
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

    for (int i = 0; i < normalDistribution.capacity(); ++i) {
        for (int j = 0; j < normalDistribution[i].Sample(draws).capacity(); ++j) {
            transposedPosterior[i][j] = normalDistribution[i].Sample(draws)[j];
        }
    }

    std::vector<std::vector<double>> init_population = transpose(transposedPosterior);
    std::vector<double> likelihoods(draws, 0);

    matplotlibcpp::plot(real_y);
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

        matplotlibcpp::plot(y_sim);

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
        double max = *max_element(likelihoods.begin(), likelihoods.end());
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



//       resampling_indexes = np.random.choice(np.arange(draws), size=draws, p=weights) # veletlen valasztas, implementalni kell...
//            posterior = posterior[resampling_indexes]
//    priors = priors[resampling_indexes]
//    likelihoods = likelihoods[resampling_indexes]
//    tempered_logp = priors + likelihoods * beta
//    acc_per_chain = acc_per_chain[resampling_indexes]
//    scalings = scalings[resampling_indexes]


    }



//#region update_proposal
//#  === smc.update_proposal() ==================================================================
//       cov = np.cov(posterior, bias=False, rowvar=0) # <- kovariancia matrix, ezt implementalni kell...
//            cov = np.atleast_2d(cov)
//    cov += 1e-6 * np.eye(cov.shape[0])
//    if np.isnan(cov).any() or np.isinf(cov).any():
//    raise ValueError('Sample covariances not valid! Likely "draws" is too small!')
//    proposal = MultivariateNormalProposal(cov)
//#endregion
//
//
//    ' ______________________________________________'
//    ' |                                            |'
//    ' | Algoritmus hangolasa                       |'
//    ' |____________________________________________|'
//#region tune
//#  === smc.tune() ============================================================================
//    if stage > 0:
//    ave_scaling = np.exp(np.log(scalings.mean()) + (acc_per_chain.mean() - 0.234))
//    scalings = 0.5 * (ave_scaling + np.exp(np.log(scalings) + (acc_per_chain - 0.234)))
//
//    if tune_steps:
//        acc_rate = max(1.0 / proposed, acc_rate)
//    n_steps = min(max_steps, max(2, int(np.log(1 - p_acc_rate) / np.log(1 - acc_rate))))
//    proposed = draws * n_steps
//#endregion
//
//
//    ' ______________________________________________'
//    ' |                                            |'
//    ' | A monte carlo lancok futtatasa             |'
//    ' |____________________________________________|'
//    ' A metrop_kernel() fuggveny itt 1db MCMC lancot futtat le'
//#region mutate
//#  === smc.mutate() ===========================================================================
//
//       parameters = (proposal, scalings, n_steps, beta)
//
//    results = [metrop_kernel(x, posterior[draw], tempered_logp[draw], priors[draw], likelihoods[draw], draw, *parameters) for draw in range(draws)]
//
//    posterior, acc_list, priors, likelihoods = zip(*results)
//
//
//    posterior = np.array(posterior)
//    priors = np.array(priors)
//    likelihoods = np.array(likelihoods)
//    acc_per_chain = np.array(acc_list)
//    acc_rate = np.mean(acc_list)
//# endregion
//
//    stage += 1
//

    matplotlibcpp::show();

    return 0;
}
