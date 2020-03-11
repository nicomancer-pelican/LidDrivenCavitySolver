#include <iostream>
#include <iomanip>
using namespace std;

#include "LidDrivenCavity.h"

//print functions for testing
void printMatrix(double* s, int Nx, int Ny){
    int k = 0;
    for(int j=0; j<Ny; j++){
        for(int i=0; i<Nx; i++){
            cout << setw(7) << setprecision(2) << fixed << *(s+k) << " ";
            k++;
        }
        cout << endl;
    }
    cout << endl;
}




int main(int argc, char **argv)
{
    //inputs here for now
    double xlen = 0.8;
    double ylen = 0.8;
    int nx = 5;
    int ny = 5;
    double deltat = 0.1;
    double finalt = 1;
    double Re = 100;
    
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();

    // Configure the solver here...
    solver->SetDomainSize(xlen, ylen);
    solver->SetGridSize(nx, ny);
    solver->SetTimeStep(deltat);
    solver->SetFinalTime(finalt);
    solver->SetReynoldsNumber(Re);
    
    solver->Initialise();
    
    // Run the solver
    solver->Integrate();
    
    solver->FirstPart();
    double* v = solver->getV();
    
    printMatrix(v,nx,ny);
    cout << endl;
    
    
    delete solver;
	return 0;
}