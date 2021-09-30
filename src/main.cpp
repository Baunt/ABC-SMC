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

    //simulated spectrum
    int npix = 256;
    std::vector<double> x = linspace(0, 1, npix);

    std::vector<double> noise = getDistribution(0.0, 1.0, npix);

    PeakModel gaussianPeakModel = PeakModel(x, real_x0, real_fwhm0, real_intensity0, npix);
    std::vector<double> gaussianModel = gaussianPeakModel.Gaussian();

    PeakModel lorentzianPeakModel = PeakModel(x, real_x1, real_fwhm1, real_intensity1, npix);
    std::vector<double> lorentzianModel = gaussianPeakModel.Lorenzt();

    std::vector<double> real_y(256);
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

    /* [ 0 0 0 0 0 0 ]
       [ 0 0 0 0 0 0 ]
       [ 0 0 0 0 0 0 ]
       [ 0 0 0 0 0 0 ]
       [ 0 0 0 0 0 0 ]

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

    return 0;
}
