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

    Eigen::ArrayX<double> realSpectrum(6);
    realSpectrum << 0.67, 0.11, 2.3, 0.51, 0.20, 0.5;
    std::vector<PeakType> peaks {Gauss, Lorentz};
    Eigen::ArrayX<double> energy = Eigen::VectorXd::LinSpaced(npix, 0, 1);
    Eigen::ArrayX<double> intensity = getSimulatedSpectrum(realSpectrum, peaks, 256, true);
    SpectrumModel spectrumModel = SpectrumModel(energy, intensity);
    spectrumModel.SetPeakList(peaks);

    AbcSmcFit().Fit(spectrumModel);

    return 0;
}
