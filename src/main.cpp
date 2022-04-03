#include "util.h"
// #include "../include/third-party-library/matplotlib-cpp/matplotlibcpp.h"
#include "spectrum_model.h"
#include "json_handler.h"
#include "abc_smc_fit.h"

int main(int argc, char** argv) {

    Config config = JsonHandler().LoadConfigFile();
    JobData jobData = JsonHandler().LoadJobFile(config.jobFilePath);
    Spectrum realSpectrum = JsonHandler().LoadMeasuredSpectrum(config.measuredSpectrumPath);

    if (config.simulated){
        int npix = 256;
        Eigen::ArrayX<double> simulatedSpectrumParameters(6);
        simulatedSpectrumParameters << 0.67, 0.11, 2.3, 0.51, 0.20, 0.5; //x0 fwhm0 int0 x1 fwhm1 int2
        std::vector<PeakType> peaks {Lorentz, Gauss};
        Eigen::ArrayX<double> energy = Eigen::VectorXd::LinSpaced(npix, 0, 1);
        Eigen::ArrayX<double> intensity = getSimulatedSpectrum(simulatedSpectrumParameters, peaks, npix, true);
        SpectrumModel spectrumModel = SpectrumModel(energy, intensity);
        spectrumModel.SetPeakList(peaks);
        AbcSmcFit().Fit(spectrumModel, config.simulated, jobData.nparams, jobData.draws, jobData.epsilon, jobData.threshold, jobData.acc_rate, jobData.n_steps, jobData.p_acc_rate, jobData.tune_steps, jobData.factor, jobData.rngSeed, realSpectrum.peakModel);
    } else{
        Eigen::ArrayX<double> energy = Eigen::VectorXd::LinSpaced(realSpectrum.spectrum.size(), 0, 1);
        SpectrumModel spectrumModel = SpectrumModel(energy, realSpectrum.spectrum);
        std::vector<PeakType> peaks;
        for (auto peak:realSpectrum.peakModel) {
            peaks.push_back(peak.peak);
        }
        spectrumModel.SetPeakList(peaks);
        AbcSmcFit().Fit(spectrumModel, config.simulated, jobData.nparams, jobData.draws, jobData.epsilon, jobData.threshold, jobData.acc_rate, jobData.n_steps, jobData.p_acc_rate, jobData.tune_steps, jobData.factor, jobData.rngSeed, realSpectrum.peakModel);
    }
    return 0;
}
