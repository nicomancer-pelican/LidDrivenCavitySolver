#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <getopt.h>
#include <map>
using namespace std;

#include "LidDrivenCavity.h"
#include "PoissonSolver.h"

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

//GET ARGUMENTS
//instructions to user
void usage(){
    cout << "required parameters: \n"
            "--Lx <number> length of the domain in the x direction.\n"
            "--Ly <number> length of the domain in the y direction.\n"
            "--Nx <number> number of grid points in the x direction.\n"
            "--Ny <number> number of grid points in the y direction.\n"
            "--Px <number> number of partitions in the x direction (parallel)"
            "--Py <number> number of partitions in the y direction (parallel)"
            "--dt <number> time step size.\n"
            "--T  <number> final time.\n"
            "--Re <number> reynolds number.";
    exit(1);
}
//function to get arguments
std::map<string, double> getArgs(int argc, char **argv){
    std::map<string, double> args;
    
    const char* const short_opts = "h";
    static struct option long_options[] = {
        {"Lx", required_argument, 0, 1},
        {"Ly", required_argument, 0, 2},
        {"Nx", required_argument, 0, 3},
        {"Ny", required_argument, 0, 4},
        {"Px", required_argument, 0, 5},
        {"Py", required_argument, 0, 6},
        {"dt", required_argument, 0, 7},
        {"T",  required_argument, 0, 8},
        {"Re", required_argument, 0, 9},
        {0, 0, 0, 0}
    };
    
    int arg;
    do{
        arg = getopt_long(argc, argv, short_opts, long_options, 0);
        
        switch(arg){
            case -1: break; //end of argument string
            case 1: if(isdigit(optarg[0]))
                        args["Lx"] = atof(optarg);
                    else usage();
                    break;
            case 2: if(isdigit(optarg[0]))
                        args["Ly"] = atof(optarg);
                    else usage();
                    break;
            case 3: if(isdigit(optarg[0]))
                        args["Nx"] = atof(optarg);
                    else usage();
                    break;
            case 4: if(isdigit(optarg[0]))
                        args["Ny"] = atof(optarg);
                    else usage();
                    break;
            case 5: if(isdigit(optarg[0]))
                        args["Px"] = atof(optarg);
                    else usage();
                    break;
            case 6: if(isdigit(optarg[0]))
                        args["Py"] = atof(optarg);
                    else usage();
                    break;
            case 7: if(isdigit(optarg[0]))
                        args["dt"] = atof(optarg);
                    else usage();
                    break;
            case 8: if(isdigit(optarg[0]))
                        args["T"] = atof(optarg);
                    else usage();
                    break;
            case 9: if(isdigit(optarg[0]))
                        args["Re"] = atof(optarg);
                    else usage();
                    break;
            default: usage(); break;    //inavalid argument
        }
    }
    while(arg != -1);
        return args;
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
    std::map<string, double> args = getArgs(argc, argv);
    
    solver->SetDomainSize(xlen, ylen);
    solver->SetGridSize(nx, ny);
    solver->SetTimeStep(deltat);
    solver->SetFinalTime(finalt);
    solver->SetReynoldsNumber(Re);
    
    solver->Initialise();
    
    // Run the solver
    solver->FirstPart();
    PoissonSolver* poisson = new PoissonSolver();
    poisson->test(*solver);
    poisson->SetA(*solver);
    poisson->SetY(*solver);
    //poisson->SetY(*solver);
    //poisson->newInteriorS(*solver);
    double* s = solver->getS();
    
    printMatrix(s,nx,ny);
    cout << endl;
    
    solver->Integrate();
    
    delete solver;
	return 0;
}