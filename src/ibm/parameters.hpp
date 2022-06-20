#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

// header file with basic ingredients for the simulation
// such as parameters and enumerators (enums)

enum Sex
{
    female = 0,
    male = 1
};

class Parameters
{
    public:
        // baseline mortality probability
        double mortality_baseline = 0.1;

        // number of males and females per patch
        int n[2] = {5,5};

        // dispersal probability 
        double d[2] = {0.1,0.5};

        double envt_change[2] = {0.1,0.5};

        // maximum number of time steps the simulation
        // should be running
        unsigned long max_time_steps = 10;

        // number of pods
        unsigned int npatches = 50;

        std::string base_name = "sim_network";

        double mu_pp = 0.01;
        double mu_pc = 0.01;
        double mu_pr = 0.01;
        double sdmu = 0.02;
};

#endif
