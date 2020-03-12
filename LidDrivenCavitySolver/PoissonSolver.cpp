#include "PoissonSolver.h"

//LAPACK stuff
#define F77NAME(x) x##_
extern "C" {
    void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                        const int& lda, int * ipiv, double * B,
                        const int& ldb, int& info);
}

//CONSTRUCTORS
PoissonSolver::PoissonSolver(){
}

PoissonSolver::~PoissonSolver(){
}



//MEMBER FUNCTIONS
void PoissonSolver::test(const LidDrivenCavity& LDC){
    
}

void PoissonSolver::SetA(const LidDrivenCavity& LDC){
    const double deltaX = LDC.Lx/(LDC.Nx-1);
    const double deltaY = LDC.Ly/(LDC.Ny-1);
    /*//create matrix A - needs to be transposed for LAPACK
    double A[7][(Nx-2)*(Ny-2)];
    a = &A[0][0];
    for(int i=0; i<(Nx-2)*(Ny-2); i++){
        *(a + i) = -1/(deltaY*deltaY);
        *(a + 3*(Ny - 2) + i) = 0.0;
        *(a + 6*(Ny - 2) + i) = -1/(deltaX*deltaX);
        *(a + 9*(Ny - 2) + i) = 2*((1/deltaX)+(1/deltaY));
        //symmetric bit
        *(a + 12*(Ny - 2) + i) = -1/(deltaX*deltaX);
        *(a + 15*(Ny - 2) + i) = 0.0;
        *(a + 18*(Ny - 2) + i) = -1/(deltaY*deltaY);
    }*/
    //non-banded version for now
    //double A[(LDC.Nx-2)*(LDC.Ny-2)][(LDC.Nx-2)*(LDC.Ny-2)];
    a = new double[(LDC.Nx-2)*(LDC.Ny-2)*(LDC.Nx-2)*(LDC.Ny-2)]; //&A[0][0];
    for(int i=0; i<(LDC.Nx-2)*(LDC.Ny-2)*(LDC.Nx-2)*(LDC.Ny-2); i++){
        *(a+i) = 0.0;
    }
    for(int i=0; i<(LDC.Nx-2)*(LDC.Ny-2); i++){
        *(a + i + i*(LDC.Nx-2)*(LDC.Ny-2)) = 2*((1/deltaX)+(1/deltaY));
    }
    for(int i=0; i<(LDC.Nx-2)*(LDC.Ny-3); i++){
        *(a + 3 + i + i*(LDC.Nx-2)*(LDC.Ny-2)) = -1/(deltaY*deltaY);
        *(a + 27 + i + i*(LDC.Nx-2)*(LDC.Ny-2)) = -1/(deltaY*deltaY);
    }
    for(int i=0; i<(LDC.Nx-2)*(LDC.Ny-2)-1; i++){
        *(a + 1 + i + i*(LDC.Nx-2)*(LDC.Ny-2)) = -1/(deltaX*deltaX);
        *(a + 9 + i + i*(LDC.Nx-2)*(LDC.Ny-2)) = -1/(deltaX*deltaX);
    }
    for(int i=0; i<(LDC.Nx-2); i++){
        *(a + (i+1)*(LDC.Nx-2)*(1 +  (LDC.Ny-2)*(LDC.Nx-3)) + i*(LDC.Nx-2)*(LDC.Ny-2)) = 0.0;
        *(a + (i+1)*((LDC.Nx-2)*(LDC.Nx-2)*(LDC.Ny-2) + (LDC.Nx-3)) + 1) = 0.0;
    }
}

void PoissonSolver::SetY(const LidDrivenCavity& LDC){
    double Y[(LDC.Nx-2)*(LDC.Ny-2)];
    y = &Y[0];
    int k = 0;
    for(int j=1; j<LDC.Ny-1; j++){
        for(int i=1; i<LDC.Nx-1; i++){
            *(y + k) = *(LDC.v + j*LDC.Nx + i);//omega[i][j];
            k ++;
        }
    }
}


void PoissonSolver::newInteriorS(const LidDrivenCavity& LDC){
    //use LAPACK to solve for x
    const int n = (LDC.Nx-2)*(LDC.Ny-2);
    //const int kl = 3;
    //const int ku = 3;
    const int nrhs = 1;
    //const int ldab = 1 + 2*kl + ku;
    //int ldb = n;
    int info;
    int * ipiv = new int[n];
    
    //general
    F77NAME(dgesv)(n, nrhs, a, n, ipiv, y, n, info);
    
    //banded
    //F77NAME(dgbsv)(n, kl, ku, nrhs, w, ldab, ipiv, y, ldb, info);
    
    //set y back into s
    int k = 0;
    for(int j=1; j<LDC.Ny-1; j++){
        for(int i=1; i<LDC.Nx-1; i++){
            *(LDC.s + k) = *(y + j*LDC.Nx + i);//omega[i][j];
            k ++;
        }
    }
    
    delete ipiv;
}