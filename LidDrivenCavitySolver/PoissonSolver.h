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
    void SetA(int Nx, int Ny, double*a);

private:
    double* a = nullptr;
};

#endif //POISSON_SOLVER_H