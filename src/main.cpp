#include <random>
#include "util.h"
#include "peak_model.h"
// #include "../include/third-party-library/matplotlib-cpp/matplotlibcpp.h"
#include "probability_distribution.h"
#include "pcg_random_generator.h"
#include "spectrum_model.h"
#include "abc_smc_fit.h"
#include <algorithm>

PcgRandomGenerator* PcgRandomGenerator::singleton_= nullptr;
extern bool DEBUG;

int main(int argc, char** argv) {
    if (argc > 1)
        DEBUG = true;
    int npix = 256;

    std::map<std::string, Eigen::ArrayX<double>> realSpectrum;
    //peak 1
    Eigen::ArrayX<double> peak1(3);
    peak1 << 0.67, 0.11, 2.3;
    //peak 2
    Eigen::ArrayX<double> peak2(3);
    peak2 << 0.51, 0.20, 0.5;
    realSpectrum.insert(std::pair<std::string, Eigen::ArrayX<double>>("Gaussian", peak1));
    realSpectrum.insert(std::pair<std::string, Eigen::ArrayX<double>>("Lorentzian", peak2));

    Eigen::ArrayX<double> real_y(npix);
    Eigen::ArrayX<double> y_sim(npix);

    //simulated spectrum
    SpectrumModel simulatedSpectrumModel = SpectrumModel(npix);
    simulatedSpectrumModel.calculate(realSpectrum, true);

    AbcSmcFit().Fit(simulatedSpectrumModel);

    return 0;
}
