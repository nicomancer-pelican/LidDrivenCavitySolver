#ifndef LIDDRIVENCAVITY_SOLVER_H
#define LIDDRIVENCAVITY_SOLVER_H

#include <string>
#include "PoissonSolver.h"
#include <iostream>
using namespace std;

class LidDrivenCavity
{
public:
    friend class PoissonSolver;
    
    //CONSTRUCTORS
    LidDrivenCavity();
    ~LidDrivenCavity();

    //MEMBER FUNCTIONS
    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);

    void Initialise();
    void Integrate();
    void FirstPart();

    // Add any other public functions

    
    //Step 1 member functions
    void horizontalBC(); //vorticity boundary conditions along the top and bottom at time t
    void verticalBC();   //vorticity bounday conditions along the left and right at time t
    
    //Step 2 member function
    void interiorV();    //interior voricity at time t
    
    //Step 3 member function
    void newInteriorV(); //interior vorticity at time t+dt
    
    //getter functions for testing
    double* getV() const;
    double* getS() const;

private:
    double* v = nullptr;
    double* s = nullptr;

    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
    double deltaX;
    double deltaY;
};

//for std::cout
inline std::ostream& operator<<(std::ostream& os, const LidDrivenCavity& a){
    return os << a.getV();
}

#endif //LIDDRIVENCAVITY_SOLVER_H