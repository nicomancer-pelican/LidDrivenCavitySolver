#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <getopt.h>
#include <map>
#include <mpi.h>
using namespace std;

#include "LidDrivenCavity.h"
#include "PoissonSolver.h"

//print functions for testing
void printMatrix(double* s, int Nx, int Ny){
    int k = 0;
    for(int j=0; j<Ny; j++){
        for(int i=0; i<Nx; i++){
            cout << setw(15) << setprecision(5) << fixed << *(s+k) << " ";
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
    //gather inputs
    std::map<string, double> args = getArgs(argc, argv);
    //test user inputs for validity - look into using assert
    //Nx, Ny > 2 as there must be at least one internal point
    //check Courant–Friedrichs–Lewy condition is satisfied (Cmax = 1 as this is explicit)
    
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();

    // Configure the solver here...
    solver->SetDomainSize(args["Lx"], args["Ly"]);
    solver->SetGridSize(args["Nx"], args["Ny"]);
    solver->SetTimeStep(args["dt"]);
    solver->SetFinalTime(args["T"]);
    solver->SetReynoldsNumber(args["Re"]);
    solver->SetPartitions(args["Px"], args["Py"]);
    
    solver->Initialise();
    
    
    //initialise MPI
    int retval = MPI_Init(&argc, &argv);
    if(retval != MPI_SUCCESS){
        cout << "an error occured initialising MPI" << endl;
    }
    
    //find rank
    int Rank, Size, retval_Rank, retval_Size;
    retval_Rank = MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    retval_Size = MPI_Comm_size(MPI_COMM_WORLD, &Size);
    if(retval_Rank == MPI_ERR_COMM || retval_Size == MPI_ERR_COMM){
        cout << "invalid communicator" << endl;
        return 1;
    }
    solver->SetRank(Rank);
    
    // Run the solver
    solver->Integrate();
    
    
    double* v = solver->getV();
    int dim1 = args["Nx"]/args["Px"] + 2;
    int dim2 = args["Ny"]/args["Py"] + 2;
    double* out;
    double* disp;
    int Px = args["Px"];
    int Py = args["Py"];
    
    if(Rank == 0){
        out = new double[dim1*dim2*Size];
    }
    
    MPI_Gather(v, dim1*dim2, MPI_DOUBLE, out, dim1*dim2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if(Rank == 0){
        int r = 0;
        int k = 0;
        int s = 0;
        while(k<Px && k<Py){
            for(int j=1; j<dim2-1; j++){
                s=0;
                r=0;
                while(r<Px){
                    for(int i=1; i<dim1-1; i++){
                        cout << setw(12) << setprecision(4) << *(out + dim1*dim2*r + (i + j*dim1));
                    }
                    r++;
                    s++;
                    
                }
                cout << endl;
            }
            //r=0;
            s++;
            while(r<Py){
                s=0;
                while(s<dim2-2){
                    s++;
                    for(int i=1; i<dim1-1; i++){
                        cout << setw(12) << setprecision(4) << *(out + dim1*dim2*r + (i + s*dim1));
                    }
                    cout << endl;
                    //s++;
                    //r++;
                }
                r++;
            }
            k++;
        }
    }
    
    
    //finalise MPI
    MPI_Finalize();
    
    delete solver;
	return 0;
}