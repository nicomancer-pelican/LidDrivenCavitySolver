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
    double* SetGlobalA(double Lx, double Ly, int Nx, int Ny);
    double* SetLocalA(int Nx, int Ny, int Px, int Py, int startCol, int endCol, int startRow, int endRow, int rank);
    double* SetY(int Nx, int Ny, int Px, int Py, int startCol, int endCol, int startRow, int endRow, double* v);
    double* SetX(int Nx, int Ny, int Px, int Py, int startCol, int endCol, int startRow, int endRow);
    
    
    //step 4 member function
    double* Execute(double Lx, double Ly, int Nx, int Ny, int Px, int Py, int startCol, int endCol, int startRow, int endRow, double* v, double* s, int rank);

private:
    double* A = nullptr;    //matrix of coefficients - global
    double* a = nullptr;    //matrix of coefficients - local to processor
    double* y = nullptr;    //input: vector of vorticities, output: vector of streamfunctions
    double* x = nullptr;    //output of conjugate gradient method - the updated streamfunctions
};

#endif //POISSON_SOLVER_H

/* Poisson Solver class used to solve a set of linear equations in the form y = Ax
 * y = vector of known interior vorticities at time t+dt
 * x = vector of unknown interior streamfunctions at time t+dt
 * A = symmetric banded matrix of constant coefficients*/