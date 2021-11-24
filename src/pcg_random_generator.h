//
// Created by balint.galgoczi on 2021.11.21..
//

#ifndef ABC_SMC_ALGORITHM_PCG_RANDOM_GENERATOR_H
#define ABC_SMC_ALGORITHM_PCG_RANDOM_GENERATOR_H

#include <random>
#include "../include/third-party-library/pcg-cpp/pcg_random.hpp"


// manualisan rogziteni a seedet teszteleshez

class PcgRandomGenerator
{
protected:
    PcgRandomGenerator(const pcg32 value): value_(value)
    {
    }

    static PcgRandomGenerator* singleton_;

    pcg32 value_;

public:
    PcgRandomGenerator(PcgRandomGenerator &other) = delete;
    void operator=(const PcgRandomGenerator &) = delete;

    static PcgRandomGenerator *GetInstance(){
        pcg_extras::seed_seq_from<std::random_device> seed_source;
        pcg32 rng(seed_source);

        if(singleton_==nullptr){
            singleton_ = new PcgRandomGenerator(rng);
        }
        return singleton_;
    };

    pcg32 value() const{
        return value_;
    }
};

#endif //ABC_SMC_ALGORITHM_PCG_RANDOM_GENERATOR_H
