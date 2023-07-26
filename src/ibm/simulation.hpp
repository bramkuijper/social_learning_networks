#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "patch.hpp"
#include "parameters.hpp"

// Simulation class represents the heart of the code
// with functions governing the life cycle 
// and the spatial distribution
// of the organisms in question

class Simulation
{
    private:
        // the file to which summary data (simulation stats, parameters)
        // will be written
        std::ofstream data_file;
        
        // the file to which network data (simulation stats, parameters)
        // will be written
        std::ofstream data_file_network;

        // keep track of the time step of the simulation
        long unsigned time_step = 0;
        
        // random device which is used to generate
        // proper random seeds
        std::random_device rd;

        // store the random seed
        // we need to store this so that we can output the 
        // random seed, so that we could 'replay' the exact
        // same sequence of random numbers for debugging purposes etc
        unsigned int seed;
        
        // random number generator
        std::mt19937 rng_r;

        // uniform distribution to compare against probabilities
        std::uniform_real_distribution<double> uniform;

        // uniform distribution to sample a random patch
        std::uniform_int_distribution<int> patch_sampler;
        
        // parameter object
        // containing all the parameters for this run
        Parameters par;
        
        // for stats purposes collect
        // the mean survival probability and other things
        int n_patches_2 = 0;

        // the metapopulation consisting of pods, here called patches
        std::vector < Patch > metapop;

        // the total fitness sum over all patches and individuals
        double W_global_total = 0.0;

        // vector to save stats on the average repertoire size
        std::vector <double> latest_repertoire_sizes;

        // vector to save stats on the average numbers of neighbours with trait t
        std::vector <double> latest_pi_ts;
        
        // vector to save stats on the average numbers of neighbours with trait t
        std::vector <double> latest_pi_ts_var;

    public:
        // the simulation constructor - building a simulation object
        Simulation(Parameters const &params);

        // function that writes data headers to the output file
        void write_data_headers();
        // function that writes data to the output file
        void write_data();
        // function that writes parameters to the output file
        void write_parameters();
        // write all the networks to the output file
        void write_all_networks();

        // randomly select individual to die and then
        // seek replacement
        void death_birth();

        // change the environment of a local pod
        void environmental_change(bool const is_envt_2);

        void make_new_individual(
            int const local_patch_idx // index of patch where offspring will establish as breeder
            ,bool const offspring_is_male // offspring will be male or female 
            ,int const individual_idx // index in vector of breeders that this individual will take up
        );

        void learn();

        // function that actually runs the simulation
        void run();

        // generate the social network 
        // based on pp, pr, pc etc
        void generate_network(
                int const local_patch_idx
                ,int const patch_origin_idx
                ,int const offspring_idx
                ,Sex const offspring_sex
                ,int const mother_idx
                ,int const father_idx
                );

        // actually learn from others 
        // and increase repertoire
        void learn(
                int const local_patch_idx // patch where this individual eventually lives
                ,int const patch_of_origin_idx  // patch where individual is born (i.e., where its parents are)
                ,Sex const offspring_sex // whether individual is male or female
                ,int const individual_idx // its index in the vector of breeders
                ,int const mother_id // the idx of the mom
                ,int const father_id // the idx of the dad
                );
};

#endif
