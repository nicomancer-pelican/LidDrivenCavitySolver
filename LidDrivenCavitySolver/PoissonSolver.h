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
    void test(const LidDrivenCavity& LDC);
    void SetA(const LidDrivenCavity& LDC);
    void SetY(const LidDrivenCavity& LDC);
    
    //step 4 member function
    void newInteriorS(const LidDrivenCavity& LDC);    //interior stream function at time t+dt

private:
    double* a = nullptr;    //matrix of coefficients
    double* y = nullptr;    //input: vector of vorticities, output: vector of streamfunctions
};

#endif //POISSON_SOLVER_H

/* Poisson Solver class used to solve a set of linear equations in the form y = Ax
 * y = vector of known interior vorticities at time t+dt
 * x = vector of unknown interior streamfunctions at time t+dt
 * A = symmetric banded matrix of constant coefficients*/