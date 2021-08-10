#include "util.h"
#include <random>

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
    //TODO what means x??
    std::vector<double> simulatedSpectrum = linspace(0,1, npix);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);
    int noise[256]={};

    for (int i=0; i<npix; ++i) {
        double number = distribution(generator);
        if ((number>=0.0)&&(number<npix)) ++noise[int(number)];
    }



    return 0;
}
