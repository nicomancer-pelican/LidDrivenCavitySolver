#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

using namespace std;

class PoissonSolver
{
public:
    //LidDrivenCavity* temp = new LidDrivenCavity();
    //CONSTRUCTORS
    PoissonSolver();
    ~PoissonSolver();
    
    //MEMBER FUNCTIONS
    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetV(double* V);
    void SetS(double* S);
    
    void SetA();
    void SetY();
    void SetX();
    
    //step 4 member function
    void newInteriorS();    //interior stream function at time t+dt
    
    double* getX() const;

private:
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double* v = nullptr;
    double* s = nullptr;
    double* a = nullptr;    //matrix of coefficients
    double* y = nullptr;    //input: vector of vorticities, output: vector of streamfunctions
    double* x = nullptr;    //output of conjugate gradient method - the updated streamfunctions
};

#endif //POISSON_SOLVER_H

/* Poisson Solver class used to solve a set of linear equations in the form y = Ax
 * y = vector of known interior vorticities at time t+dt
 * x = vector of unknown interior streamfunctions at time t+dt
 * A = symmetric banded matrix of constant coefficients*/