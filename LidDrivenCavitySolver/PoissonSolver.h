#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

using namespace std;
#include "LidDrivenCavity.h"
class LidDrivenCavity;

class PoissonSolver
{
public:
    //LidDrivenCavity* temp = new LidDrivenCavity();
    //CONSTRUCTORS
    PoissonSolver();
    ~PoissonSolver();
    
    //MEMBER FUNCTIONS
    void SetA(LidDrivenCavity& LDC, double* a);
    void SetY(LidDrivenCavity& LDC, double* y, double* v);
    
    //step 4 member function
    void newInteriorS(LidDrivenCavity& LDC);    //interior stream function at time t+dt

private:
    double* a = nullptr;
    double* y = nullptr;
};

#endif //POISSON_SOLVER_H

/* Poisson Solver class used to solve a set of linear equations in the form y = Ax
 * y = vector of known interior vorticities at time t+dt
 * x = vector of unknown interior streamfunctions at time t+dt
 * A = symmetric banded matrix of constant coefficients*/