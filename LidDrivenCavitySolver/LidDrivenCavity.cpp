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
    Lx = xlen;
    Ly = ylen;
}

void LidDrivenCavity::SetGridSize(int nx, int ny){
    Nx = nx;
    Ny = ny;
}

void LidDrivenCavity::SetTimeStep(double deltat){
    dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt){
    T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re){
    Re = re;
}

/*//getters
int LidDrivenCavity::GetNx() const{
    return Nx;
}

int LidDrivenCavity::GetNy() const{
    return Ny;
}

double LidDrivenCavity::GetLx() const{
    return Lx;
}

double LidDrivenCavity::GetLy() const{
    return Ly;
}*/

//other member functions
void LidDrivenCavity::Initialise(){
    //pointers to matrices
    v = new double[Nx*Ny];
    s = new double[Nx*Ny];
}

void LidDrivenCavity::Integrate(){
    PoissonSolver* poisson = new PoissonSolver();
    
    for(double t=0.0; t<T; t+=dt){
        boundaryConditions();
        interiorV();
        newInteriorV();
        updateS(poisson->Execute(Lx, Ly, Nx, Ny, v, s));
    }
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

//STEP 4 MEMBER FUNCTION
void LidDrivenCavity::updateS(double* x){
    int k = 0;
    for(int j=1; j<Ny-1; j++){
        for(int i=1; i<Nx-1; i++){
            *(s + j*Nx + i) = *(x + k);//omega[i][j];
            k ++;
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


