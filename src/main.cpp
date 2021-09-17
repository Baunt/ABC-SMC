#include "util.h"
#include "peak_model.h"
#include "../include/third-party-library/matplotlib-cpp/matplotlibcpp.h"
#include "concrete_spectrum_model.h"


int main() {

    //peak 1
    double real_x0 = 0.67;
    double real_fwhm0 = 0.11;
    double real_intensity0 = 2.3;
    //peak 2
    double real_x1 = 0.51;
    double real_fwhm1 = 0.20;
    double real_intensity1 = 0.5;

    //measurement
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

//    matplotlibcpp::plot(real_y);
//    matplotlibcpp::show();
    auto* builder = new ConcreteSpectrumModel();

    SpectrumModelParameters.x = linspace(0,1,npix);
    SpectrumModelParameters.x0 = real_x0;
    SpectrumModelParameters.fwhm = real_fwhm0;
    SpectrumModelParameters.intensity = real_intensity0;
    builder->GaussianModel(SpectrumModelParameters);

    SpectrumModelParameters.x = linspace(0,1,npix);
    SpectrumModelParameters.x0 = real_x0;
    SpectrumModelParameters.fwhm = real_fwhm0;
    SpectrumModelParameters.intensity = real_intensity0;
    builder->LorentzianModel(SpectrumModelParameters);
    return 0;
}
