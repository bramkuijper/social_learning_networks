#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

// header file with basic ingredients for the simulation
// such as parameters and enumerators (enums)


// make a data type that translates 
// male and female to array indices
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
        int n[2] = {10,10};

        // sex-biased dispersal
        // dispersal probability 
        double d[2] = {0.1,0.5};

        double envt_change[2] = {0.1,0.5};

        // number of learning attempts to create proficiency in any trait
        int n_learning_attempts = 25;

        // if an offspring is immigrant, should it learn
        // remotely or locally

        // number of traits in repertoire
        int n_traits = 30;

        // variable to specify max search attempts
        // e.g., when changing an environment
        int max_search_attempts = 30;

        // environmental switch rates
        double switch_rate[2] = {0.0,0.0};

        // maximum number of time steps the simulation
        // should be running
        unsigned long max_time_steps = 10;

        // number of pods
        unsigned int npatches = 50;

        // file names:
        // the file name for the averages and other stats
        std::string base_name = "sim_network_stats";
        // the file name for the network data
        std::string base_name_matrix_file = "sim_network_data";

        // mutation probabilities for each of the traits
        // see individual.hpp for a description
        double mu_il = 0.01;
        double mu_pp = 0.01;
        double mu_pc = 0.01;
        double mu_pr = 0.01;

        // standard deviation of the mutational effect size distribution
        double sdmu = 0.02;
};

#endif
