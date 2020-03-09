#ifndef LIDDRIVENCAVITY_SOLVER_H
#define LIDDRIVENCAVITY_SOLVER_H

#include <string>
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

    void Initialise(double xlen, double ylen, int nx, int ny, double deltat, double finalt, double re);
    void Integrate();

    // Add any other public functions
    void SetMatrices(int Nx, int Ny, double* v, double* s);
    
    //Step 1 member functions
    void topBC(double* v, double* s, int Nx, int Ny, double deltaY);        //vorticity boundary conditions along the top and bottom at time t
    void verticalBC(double* v, double* s, int Nx, int Ny, double deltaX);   //vorticity bounday conditions along the left and right at time t
    
    //Step 2 member function
    void interiorV(double* v, double* s, int Nx, int Ny, double deltaX, double deltaY); //interior voricity at time t
    
    //Step 2 member function
    void newInteriorV(double* v, double* s, int Nx, int Ny, double deltaX, double deltaY, double Re, double dt); //interior vorticity at time t+dt

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
};


#endif //LIDDRIVENCAVITY_SOLVER_H