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

void LidDrivenCavity::Initialise(double xlen, double ylen, int nx, int ny, double deltat, double finalt, double re){
}

void LidDrivenCavity::Integrate(){
}

void LidDrivenCavity::SetMatrices(int Nx, int Ny, double* v, double* s){
    double omega[Nx][Ny];
    v = &omega[0][0];
    double psi[Nx][Ny];
    s = &psi[0][0];
}


//STEP 1 MEMBER FUNCTIONS
void horizontalBC(double* v, double* s, int Nx, int Ny, double deltaY){
    for(int i=0; i<Nx; i++){
        *(v + Nx*(Ny-1) + i) = (*(s + (Ny-1)*Nx + i) - *(s+(Ny-2)*Nx + i))*(2/(deltaY*deltaY)) - (2/deltaY);  //top boundary condition
        *(v + i) = (*(s + i) - *(s + Nx + i))*(2/(deltaY*deltaY));  //bottom boundary condition
    }
}

void verticalBC(double* v, double* s, int Nx, int Ny, double deltaX){
    for(int j=0; j<Ny; j++){
        *(v + j*Nx) = (*(s + j*Nx) - *(s + j*Nx + 1))*(2/(deltaX*deltaX));  //left boundary condition
        *(v + Nx*(j+1) - 1) = (*(s + Nx*(j+1) - 1) - *(s + Nx*(j+1) - 2))*(2/(deltaX*deltaX));  //right boundary condition
    }
}

//STEP 2 MEMBER FUNCTION
void interiorV(double* v, double* s, int Nx, int Ny, double deltaX, double deltaY){
    for(int j=1; j<Ny-1; j++){
        for(int i=1; i<Nx-1; i++){
            *(v + j*Nx + i) = -(*(s + j*Nx + i + 1) - *(s + j*Nx + i)*2 + *(s + j*Nx + i - 1))*(2/(deltaX*deltaX))
                              -(*(s + Nx*(j+1) + i) - *(s + j*Nx + i)*2 + *(s + Nx*(j-1) + i))*(2/(deltaY*deltaY));
        }
    }
}