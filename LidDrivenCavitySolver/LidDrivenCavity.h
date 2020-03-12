#ifndef LIDDRIVENCAVITY_SOLVER_H
#define LIDDRIVENCAVITY_SOLVER_H

#include <string>
#include <iostream>

#include "PoissonSolver.h"
using namespace std;

class LidDrivenCavity
{
public:
    //CONSTRUCTORS
    LidDrivenCavity();
    ~LidDrivenCavity();

    //MEMBER FUNCTIONS
    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    
    int GetNx() const;
    int GetNy() const;
    double GetLx() const;
    double GetLy() const;

    void Initialise();
    void Integrate();
    void FirstPart();

    // Add any other public functions
    
    //Step 1
    void boundaryConditions(); //vorticity boundary conditions at time t
    
    //Step 2
    void interiorV();          //interior voricity at time t
    
    //Step 3
    void newInteriorV();       //interior vorticity at time t+dt
    
    //getter functions for testing
    double* getV() const;
    double* getS() const;

private:
    double* v = nullptr;    //vorticity matrix pointer
    double* s = nullptr;    //streamfunction matrix pointer

    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
    //double deltaX;
    //double deltaY;
    friend class PoissonSolver;
};

//for std::cout
inline std::ostream& operator<<(std::ostream& os, const LidDrivenCavity& a){
    return os << a.getS();
}

#endif //LIDDRIVENCAVITY_SOLVER_H