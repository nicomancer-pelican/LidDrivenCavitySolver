#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

using namespace std;

class PoissonSolver
{
public:
    //CONSTRUCTORS
    PoissonSolver();
    ~PoissonSolver();
    
    //MEMBER FUNCTIONS    
    double* SetA(double Lx, double Ly, int Nx, int Ny);
    double* SetY(int Nx, int Ny, double* v);
    double* SetX(int Nx, int Ny);
    
    //step 4 member function
    double* Execute(double Lx, double Ly, int Nx, int Ny, double* v, double* s);    //interior stream function at time t+dt

private:
    double* a = nullptr;    //matrix of coefficients
    double* y = nullptr;    //input: vector of vorticities, output: vector of streamfunctions
    double* x = nullptr;    //output of conjugate gradient method - the updated streamfunctions
};

#endif //POISSON_SOLVER_H

/* Poisson Solver class used to solve a set of linear equations in the form y = Ax
 * y = vector of known interior vorticities at time t+dt
 * x = vector of unknown interior streamfunctions at time t+dt
 * A = symmetric banded matrix of constant coefficients*/