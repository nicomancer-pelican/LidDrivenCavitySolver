#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

using namespace std;

class PoissonSolver
{
public:
    friend class LidDrivenCavity;   //to access things in LidDrivenCavity
    
    //CONSTRUCTORS
    PoissonSolver();
    ~PoissonSolver();
    
    //MEMBER FUNCTIONS
    void SetA(int Nx, int Ny, double* a);
    void SetY(int Nx, int Ny, double* y, double* v);
    void SetX(int Nx, int Ny, double* x);

private:
    double* a = nullptr;
    double* y = nullptr;
    double* x = nullptr;
};

#endif //POISSON_SOLVER_H

/* Poisson Solver class used to solve a set of linear equations in the form y = Ax
 * y = vector of known interior vorticities at time t+dt
 * x = vector of unknown interior streamfunctions at time t+dt
 * A = symmetric banded matrix of constant coefficients*/