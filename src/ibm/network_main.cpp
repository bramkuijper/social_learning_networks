#include "simulation.hpp"
#include "parameters.hpp"

// the main part of the code
int main(int argc, char **argv)
{
    Parameters params;

    Simulation sim(params);

    sim.run();
}
