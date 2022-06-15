#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "parameters.hpp"

// Simulation class represents the heart of the code
// with functions governing the life cycle 
// and the spatial distribution
// of the organisms in question

class Simulation
{
    private:
        // the file to which data (simulation stats, parameters)
        // will be written
        std::ofstream data_file;

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
        
        // parameter object
        // containing all the parameters for this run
        Parameters par;
        
        // for stats purposes collect
        // the mean survival probability and other things
        double mean_surv_prob[2] = {0.0,0.0};
        int nsurvivors[2] = {0,0};
        double mean_offspring[2] = {0.0,0.0};

    public:
        // the simulation constructor - building a simulation object
        Simulation(Parameters &params);

        // function that writes data headers to the output file
        void write_data_headers();
        // function that writes data to the output file
        void write_data();
        // function that writes parameters to the output file
        void write_parameters();

        // function that actually runs the simulation
        void run();
};

#endif
