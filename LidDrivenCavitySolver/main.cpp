#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <getopt.h>
#include <map>
#include <mpi.h>
#include <fstream>
using namespace std;

#include "LidDrivenCavity.h"
#include "PoissonSolver.h"


/* np3217 01333401
 * AERO96021 High Performance Computing 2020 Coursework
 * 
 * LID-DRIVEN CAVITY PROBLEM
 * 
 * Objectives: write a parallel numerical code for solving the voriticity-stream
 * funtion formulation of the incompressible Navier-Stokes equations in 2D using
 * the finite difference method.
 * 
 * Step 1: vorticity at time t is updated on the boundaries
 * Step 2: vorticity at time t is updated on the interior
 * Step 3: vorticity at time t+dt is computed on the interior
 * Step 4: stream function at time t+dt is computed (solve Poisson equation)
 * 
 * The algorithm was split such that step 4 is taken care of in the PoissonSolver
 * class. This is used in the LidDrivenCavity class to run the algorithm through
 * all four steps and loop through them.
 * 
 * main.cpp is run to gather the user inputs which are fed into the LidDrivenCavity
 * class and then into the PoissonClass to perform the iterative calculations. The
 * data for the final timestep, T, is then returned to main.cpp where it is then
 * reconstructed and written to their respective files.
*/

//FUNCTIONS

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

//instructions to user - accessed by -h on command line
void usage(){
    cout << "required parameters: \n"
            "--Lx <number> length of the domain in the x direction.\n"
            "--Ly <number> length of the domain in the y direction.\n"
            "--Nx <number> number of grid points in the x direction.\n"
            "--Ny <number> number of grid points in the y direction.\n"
            "--Px <number> number of partitions in the x direction (parallel).\n"
            "--Py <number> number of partitions in the y direction (parallel).\n"
            "--dt <number> time step size.\n"
            "--T  <number> final time.\n"
            "--Re <number> reynolds number.\n";
    exit(1);
}

//get arguments
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



//MAIN
int main(int argc, char **argv)
{
    //gather inputs
    std::map<string, double> args = getArgs(argc, argv);
    
    double Lx = args["Lx"];
    double Ly = args["Ly"];
    int Nx = args["Nx"];
    int Ny = args["Ny"];
    int Px = args["Px"];
    int Py = args["Py"];
    double dt = args["dt"];
    double T  = args["T "];
    double Re = args["Re"];
    
    //check inputs
    if(Nx%Px != 0 || Ny%Py != 0){
        cout << "number of processor divisions must be a factor of the number of elements" << endl;
        return 1;
    }
    if(Nx <= 2 || Ny <= 2){
        cout << "domain must have a size of 3x3 or greater" << endl;
        return 1;
    }
    if(Lx <= 0 || Ly <= 0){
        cout << "dimensions must be greater than 0" << endl;
        return 1;
    }
    /*if(T <= 0){
        cout << "final time must be greater than 0" << endl;
        return 1;
    }*/
    if(Re <=0){
        cout << "Reynolds number must be greater than 0" << endl;
        return 1;
    }
    if(dt*(Nx-1)/(Lx) > 1){
        cout << "Courant-Friedrichs-Lewy condition not satisfied - decrease the timestep" << endl;
        return 1;
    }
    
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();

    // Configure the solver
    solver->SetDomainSize(args["Lx"], args["Ly"]);
    solver->SetGridSize(args["Nx"], args["Ny"]);
    solver->SetTimeStep(args["dt"]);
    solver->SetFinalTime(args["T"]);
    solver->SetReynoldsNumber(args["Re"]);
    solver->SetPartitions(args["Px"], args["Py"]);
    
    //initialise solver
    solver->Initialise();
    
    
    //initialise MPI
    int retval = MPI_Init(&argc, &argv);
    if(retval != MPI_SUCCESS){
        cout << "an error occured initialising MPI" << endl;
    }
    
    //find rank and size
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
    
    //prep to get results
    double* v = solver->getV();
    double* s = solver->getS();
    int dim1 = args["Nx"]/args["Px"] + 2;
    int dim2 = args["Ny"]/args["Py"] + 2;
    double* outV;
    double* outS;
    
    //only allocate memory on root rank - unecessary to do so on the others
    if(Rank == 0){
        outV = new double[dim1*dim2*Size];
        outS = new double[dim1*dim2*Size];
    }
    
    //gather to root rank
    MPI_Gather(v, dim1*dim2, MPI_DOUBLE, outV, dim1*dim2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(s, dim1*dim2, MPI_DOUBLE, outS, dim1*dim2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if(Rank == 0){
        //create a file for vorticity data
        ofstream vOut("vorticity.txt", ios::out);
        vOut.precision(4);
        vOut << endl << "FINAL VORTICITY - in row major format (rows divided by a blank line)" << endl << endl;
        
        //create a file for streamfunction data
        ofstream sOut("streamfunction.txt", ios::out);
        sOut.precision(4);
        sOut << endl << endl << "FINAL STREAMFUNCTION - in row major format (rows divided by a blank line)" << endl << endl;
        
        //set counters
        int k = 0;
        int count = 0;
        int j = 0;
        int r = 0;
        
        //reconstruct global matrices
        while(count < Size){
            for(int a = 1; a<dim2-1; a++){
                j = 0;
                r = count;
                while(j<Px){
                    for(int i=1; i<dim1-1; i++){
                        vOut << setw(12) << *(outV + dim1*dim2*r + (i + a*dim1));
                        sOut << setw(12) << *(outS + dim1*dim2*r + (i + a*dim1));
                    }
                    j++;
                    r++;
                }
                //space between each row of data
                vOut << endl << endl;
                sOut << endl << endl;
                k++;
            }
            count += Px;
        }
        vOut.close();
        sOut.close();
    }
    
    
    //finalise MPI
    MPI_Finalize();
    
    delete solver;
	return 0;
}