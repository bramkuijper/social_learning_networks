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

        double envt_change[2] = {1,0.1};

        // number of learning attempts to create proficiency in any trait
        int n_learning_attempts = 100;

        // initial_network connection probability
        double p_network_init = 0.1;

        // if an offspring is immigrant, should it learn
        // remotely or locally

        // number of traits in repertoire
        int n_traits = 150;

        // variable to specify max search attempts
        // e.g., when changing an environment
        int max_search_attempts = 100;

        // environmental switch rates
        double switch_rate[2] = {0.9,1};

        // maximum number of time steps the simulation
        // should be running
        unsigned long max_time_steps = 1e08;

        // number of pods
        unsigned int n_patches = 150;

        // innovation success probability
        // see eq. (1) in Smolla & Akcay
        double gamma = 0.75;
        // see eq. (2) in Smolla & Akcay
        double sigma = 0.01;

        // file names:
        // the file name for the averages and other stats
        std::string base_name = "sim_network_stats";
        // the file name for the network data
        std::string base_name_matrix_file = "sim_network_data";

        // mutation probabilities for each of the traits
        // see individual.hpp for a description
        double mu_il = 1;
        double mu_pp = 1;
        double mu_pc = 1;
        double mu_pr = 1;

        // initial values of everything
        double il_init = 0.75;
        double pp_init = 0.75;
        double pc_init = 0.75;
        double pr_init = 0.75;

        // standard deviation of the mutational effect size distribution
        double sdmu = 0.1;

        // output interval of stats
        unsigned long data_output_interval = 1000;


};

#endif
