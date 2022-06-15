#include "simulation.hpp"
#include "patch.hpp"

// the simulation constructor which 
// initializes everything
Simulation::Simulation(Parameters &par) :
    rd{} // initialize random device, see *.hpp file
    ,seed{rd()} // initialize seed
    ,rng_r{seed} // initialize the random number generator
    ,uniform{0.0,1.0} // initialize the normal distribution
    ,data_file{params.base_name.c_str()} // initialize the data file by giving it a name
    ,par{params} // initialize the parameter data member with the constructor argument
    ,metapop{par.npatches, Patch(par.npp[0],par.npp[1])} // initialize a meta population each with n1 individuals of species 1 and n2 individuals of species 2
{} // end IBM_Mutualism

// run the simulation
Simulation::run()
{
    // write the headers to the output file
    write_data_headers();

    // we loop over all the time steps and perform a whole life cycle
    // each and every time step
    for (time_step = 0; time_step < par.max_time_steps; ++time_step)
    {

    }
} // end Simulation::run()

void Simulation::write_data_headers()
{}

void Simulation::write_data()
{}
