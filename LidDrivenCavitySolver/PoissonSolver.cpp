#include "PoissonSolver.h"

//CONSTRUCTORS
PoissonSolver::PoissonSolver(){
}

PoissonSolver::~PoissonSolver(){
}



//MEMBER FUNCTIONS
void SetA(int Nx, int Ny, double*a, double deltaX, double deltaY){
    double A[4][(Nx-2)*(Ny-2)];
    a = &A[0][0];
    
    for(int i=0; i<(Nx-2)*(Ny-2); i++){
        *(a + 9*(Ny - 2) + i) = -1/(deltaY*deltaY);
        *(a + 6*(Ny - 2) + i) = 0.0;
        *(a + 3*(Ny - 2) + i) = -1/(deltaX*deltaX);
        *(a + i) = 2*((1/deltaX)+(1/deltaY));
    }
    for(int i=Nx-2; i<(Nx-2)*(Ny-2); i+=Nx-2){
        *(a + 3*(Ny - 2) + i) = 0.0;
    }
}