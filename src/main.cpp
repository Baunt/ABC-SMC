#include <random>
#include "util.h"
// #include "../include/third-party-library/matplotlib-cpp/matplotlibcpp.h"
#include "spectrum_model.h"
#include "abc_smc_fit.h"
#include "json_handler.h"
#include <algorithm>

int main(int argc, char** argv) {

    Config config = JsonHandler().LoadConfigFile();
    JobData jobData = JsonHandler().LoadJobFile(config.jobFilePath);
    Spectrum realSpectrum = JsonHandler().LoadMeasuredSpectrum(config.measuredSpectrumPath);

    if (config.simulated){
        int npix = 256;
        Eigen::ArrayX<double> simulatedSpectrumParamteres(6);
        simulatedSpectrumParamteres << 0.67, 0.11, 2.3, 0.51, 0.20, 0.5; //x0 fwhm0 int0 x1 fwhm1 int2
        std::vector<PeakType> peaks {Lorentz, Gauss};
        Eigen::ArrayX<double> energy = Eigen::VectorXd::LinSpaced(npix, 0, 1);
        Eigen::ArrayX<double> intensity = getSimulatedSpectrum(simulatedSpectrumParamteres, peaks, npix, true);
        SpectrumModel spectrumModel = SpectrumModel(energy, intensity);
        spectrumModel.SetPeakList(peaks);
        AbcSmcFit().Fit(spectrumModel);
    } else{

    }


    return 0;
}
