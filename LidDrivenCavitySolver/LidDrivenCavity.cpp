#include "LidDrivenCavity.h"

//CONSTRUCTORS//////////////////////////////////////////////////////
//default constructor
LidDrivenCavity::LidDrivenCavity(){
}

//default destructor
LidDrivenCavity::~LidDrivenCavity(){
}



//MEMBER FUNCTIONS/////////////////////////////////////////////////
void LidDrivenCavity::SetDomainSize(double xlen, double ylen){
    xlen = Lx;
    ylen = Ly;
}

void LidDrivenCavity::SetGridSize(int nx, int ny){
    nx = Nx;
    ny = Ny;
}

void LidDrivenCavity::SetTimeStep(double deltat){
    deltat = dt;
}

void LidDrivenCavity::SetFinalTime(double finalt){
    finalt = T;
}

void LidDrivenCavity::SetReynoldsNumber(double re){
    re = Re;
}

void LidDrivenCavity::Initialise(){
    //pointers to matrices
    v = new double[Nx*Ny];
    s = new double[Nx*Ny];
}

void LidDrivenCavity::Integrate(){
}

void LidDrivenCavity::FirstPart(){
    boundaryConditions();
    interiorV();
    newInteriorV();
}


//STEP 1 MEMBER FUNCTIONS
void LidDrivenCavity::boundaryConditions(){
    const double deltaX = Lx/(Nx-1);
    const double deltaY = Ly/(Ny-1);
    
    for(int i=0; i<Nx; i++){
        //top boundary condition
        *(v + Nx*(Ny-1) + i) = (*(s + (Ny-1)*Nx + i) - *(s+(Ny-2)*Nx + i))*(2/(deltaY*deltaY)) - (2/deltaY);
        //bottom boundary condition
        *(v + i) = (*(s + i) - *(s + Nx + i))*(2/(deltaY*deltaY));
    }
    
    for(int j=0; j<Ny; j++){
        //left boundary condition
        *(v + j*Nx) = (*(s + j*Nx) - *(s + j*Nx + 1))*(2/(deltaX*deltaX));
        //right boundary condition
        *(v + Nx*(j+1) - 1) = (*(s + Nx*(j+1) - 1) - *(s + Nx*(j+1) - 2))*(2/(deltaX*deltaX));
    }
}

//STEP 2 MEMBER FUNCTION
void LidDrivenCavity::interiorV(){
    const double deltaX = Lx/(Nx-1);
    const double deltaY = Ly/(Ny-1);
    for(int j=1; j<Ny-1; j++){
        for(int i=1; i<Nx-1; i++){
            *(v + j*Nx + i) = -(*(s + j*Nx + i + 1) - *(s + j*Nx + i)*2 + *(s + j*Nx + i - 1))*(2/(deltaX*deltaX))
                              -(*(s + Nx*(j+1) + i) - *(s + j*Nx + i)*2 + *(s + Nx*(j-1) + i))*(2/(deltaY*deltaY));
        }
    }
}

//Setp 3 MEMBER FUNCTION
void LidDrivenCavity::newInteriorV(){
    const double deltaX = Lx/(Nx-1);
    const double deltaY = Ly/(Ny-1);
    for(int j=1; j<Ny-1; j++){
        for(int i=1; i<Nx-1; i++){
            *(v + j*Nx + i) = *(v + j*Nx + i) + dt*(
                               (1/Re) * (((*(v + j*Nx + i + 1) - *(v + j*Nx + i)*2 + *(v + j*Nx + i - 1)) / (deltaX*deltaX))
                                        -((*(v + Nx*(j+1) + i) - *(v + j*Nx + i)*2 + *(v + Nx*(j-1) + i)) / (deltaY*deltaY)))
                              -(((*(s + Nx*(j+1) + i) - *(s + Nx*(j-1) + i)) / (2*deltaY))
                               *((*(v + j*Nx + i + 1) - *(v + j*Nx + i - 1)) / (2 *deltaX)))
                              +(((*(s + j*Nx + i + 1) - *(s + j*Nx + i - 1)) / (2*deltaX))
                               *((*(v + Nx*(j+1) + i) - *(v + Nx*(j-1) + i)) / (2*deltaY)))
                                );
        }
    }
}



//getter functions for testing
double* LidDrivenCavity::getV() const{
    return v;
}

double* LidDrivenCavity::getS() const{
    return s;
}


