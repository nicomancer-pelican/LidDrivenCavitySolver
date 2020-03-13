#include "PoissonSolver.h"
#include "cblas.h"

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
    
    //create matrix A - symmetric upper triangle
    int dim = (LDC.Nx-2)*(LDC.Ny-2);
    double A[dim*dim];
    double* a = &A[0];
    for(int j=0; j<dim; j++){
        for(int i=0; i<dim; i++){
            if(i==j){
                *(a + j*dim + i) = 2*((1/deltaX)+(1/deltaY));
                *(a + j*dim + i + (LDC.Nx-2)) = -1/(deltaY*deltaY);
                if(i%(LDC.Nx-2) != 0){
                    *(a + (j-1)*dim + i) = -1/(deltaX*deltaX);
                }
            }
        }
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

void PoissonSolver::SetX(const LidDrivenCavity& LDC){
    double X[(LDC.Nx-2)*(LDC.Ny-2)];
    x = &X[0];
}

void PoissonSolver::newInteriorS(const LidDrivenCavity& LDC){
    int n = (LDC.Nx-2)*(LDC.Ny-2);
    double R[n] = {0.0};
    double* r = &R[0];
    double P[n] = {0.0};
    double* p = &P[0];
    double T[n] = {0.0};
    double* t = &T[0];  //temp vector for r vector
    double X[n] = {0.0};
    double* x = &X[0];  //x_0
    int k = 0;
    double alpha;
    double beta;
    double eps;
    double tol = 0.00001;
    
    //setup things
    cblas_dcopy(n, y, 1, r, 1); //r_0 = b (i.e. y)
    cblas_dsymv(CblasRowMajor, CblasUpper, n, -1.0, a, n, x, 1, 1.0, r, 1); //r_0 = y - Ax_0
    cblas_dcopy(n, r, 1, p, 1); //p_0 = r_0
    //loop
    do{
        cblas_dsymv(CblasRowMajor, CblasUpper, n, 1.0, a, n, p, 1, 0.0, t, 1); //t = A*p_k
        alpha = cblas_ddot(n, t, 1, p, 1);          //alpha = trans(p_k) * A * p_k
        alpha = cblas_ddot(n, r, 1, r, 1) / alpha;  //alpha_k = trans(r_k) * r_k / trans(p_k) * A * p_k
        beta  = cblas_ddot(n, r, 1, r, 1);          //trans(r_k) * r_k
        
        cblas_daxpy(n,  alpha, p, 1, x, 1);         //x_k+1 = x_k + alpha_k * p_k
        cblas_daxpy(n, -alpha, t, 1, r, 1);         //r_k+1 = r_k - alpha_k * A * p_k
        
        eps = cblas_dnrm2(n, r, 1);     //euclidean norm
        if(eps < tol*tol){
            break;
        }
        
        beta = cblas_ddot(n, r, 1, r, 1) / beta;    //trans(r_k+1) * r_k+1 / trans(r_k) * r_k
        
        cblas_dcopy(n, r, 1, t, 1);                 //copy r into temp, t
        cblas_daxpy(n, beta, p, 1, t, 1);           //P_k+1 = r_k+1 + beta * P_k
        cblas_dcopy(n, t, 1, p, 1);                 //copy temp, t into p
        
        k++;
    } while(k<5000);
    
    delete r;
    delete p;
    delete t;
}