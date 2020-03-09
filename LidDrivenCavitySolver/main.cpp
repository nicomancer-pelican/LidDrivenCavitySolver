#include <iostream>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    //inputs here for now
    double xlen = 1.0;
    double ylen = 1.0;
    int nx = 10;
    int ny = 10;
    double deltat = 0.1;
    double finalt = 1;
    double Re = 1000;
    
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();

    // Configure the solver here...
    // ...
    solver->Initialise(xlen, ylen, nx, ny, deltat, finalt, Re);

    // Run the solver
    solver->Integrate();

	return 0;
}