//
// Created by balint.galgoczi on 2021.11.21..
//

#ifndef MATPLOTLIB_CPP_SPECTRUM_MODEL_H
#define MATPLOTLIB_CPP_SPECTRUM_MODEL_H


#include "probability_distribution.h"
#include "peak_model.h"
#include <map>
#include <list>

//A Spectrumodel-t én a következőképpen képzelem el:
//
//- elég csak a szimulált spektrumnak létrehozni egy ilyen osztályt az smc osztályban
//Most a tesztelésnél a mérés szimulációját (amit mi igazi spektrumnak nevezünk) is így számoljuk ki, de igazából az
//nem lesz dolga az smc osztálynak. Magát a mérési eredmény szimulációját egyelőre csinálhatjuk a main-ben (itt létrehozthatunk
//        amúgy egy másik a SpectrumModell osztályt erre a célra, de ez csak arra kell hogy generáljunk egy mérést, utána kuka), később
//        pedig majd valami csv-ből kellene beolvasni a mért, vagy szimulált mérési spektrumokat.
//
//- amint azt írtam a kommentekben is, nem hiszem hogy van értelme eltárolni bármilyen belső állapotot ami nem állandó
//(pl. szimulált spektrum aminek most "spectrum" a neve)
//- amit viszont érdemes lenne tárolni azok a nem változó dolgok:
//        * a mért hullámossz tengely (ez most energy) és mért spektrum (ez most hiányzik).
//        * a modell: a csúcstípusok listája: ez lehetne akár egy vector<string> is,
//        de a Gaussian és Lorentz típusoknak érdemes lenne egy enumot bevezetni,
//        vagy számokkal kódolni őket
//        * a Calculate a csúcstípusok listáját iterálná, amikor meghívod a függvényt akkor ezeket már nem kellene megadni
//        * én a paramétereket egyetlen tömbbe tenném "ömlesztve" és amikor iteráljuk a csúcsok listáját,
//        akkor mindig léptetni hogy melyik elemeket vesszük ki : parameters.segment(i*3,3)
//        Ez nem olyan elegáns, de így csak a paramétereket tartalmazó 2d tömbből elég lenne egy oszlopot
//        odaadni a Calculate függvénynek
//        * a Calculate visszatérési értéke szerintem legyen a kiszámolt spektrum ArrayX<double>
//        * illetve be lehetne vezetni egy double MeanAbsError(ArrayX params) függvényt is, ami rögtön
//        kiszámolná a hibát, és ekkor ezt már nem SmcFit-ben kellene csinálni. Ez meghívhatná a
//        Calculate függvényt (amit én azért meghagynék publikusnak, de akkor nem kellene hívogatni
//        az smcfit osztályból)
//

enum PeakType{
    Gauss,
    Lorentz
};

class SpectrumModel {
private:
    std::vector<PeakType> peaks;
    int npix;

public:
    std::vector<NormalDistribution> InitialGuess;
    Eigen::ArrayXX<double> GenerateInitialPopulation(int nsamples, int nparams);
    Eigen::ArrayX<double> Calculate(Eigen::ArrayX<double> parameters, bool withNoise);
    double ErrorCalculation(Eigen::ArrayX<double> diffSpectrum);
    void SetPeakList(std::vector<PeakType>& peaks);
    Eigen::ArrayX<double> energy;
    Eigen::ArrayX<double> intensity;


    SpectrumModel(const Eigen::ArrayX<double>& energy, const Eigen::ArrayX<double>& intensity){
        this->npix = energy.size();
        this->energy = energy;
        this->intensity = intensity;
    }
};


#endif //MATPLOTLIB_CPP_SPECTRUM_MODEL_H
