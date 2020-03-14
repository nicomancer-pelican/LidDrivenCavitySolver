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
double* PoissonSolver::SetA(double Lx, double Ly, int Nx, int Ny){
    const double deltaX = Lx/(Nx-1);
    const double deltaY = Ly/(Ny-1);
    
    //create matrix A - symmetric upper triangle
    int dim = (Nx-2)*(Ny-2);
    double A[dim*dim] = {0.0};
    a = &A[0];
    for(int j=0; j<dim; j++){
        for(int i=0; i<dim; i++){
            if(i==j){
                *(a + j*dim + i) = 2*((1/deltaX)+(1/deltaY));
                if(j*dim + i + (Nx-2) < dim*dim){
                    *(a + j*dim + i + (Nx-2)) = -1/(deltaY*deltaY);
                }
                if(i%(Nx-2) != 0 && (j-1)*dim + i < dim*dim){
                    *(a + (j-1)*dim + i) = -1/(deltaX*deltaX);
                }
            }
        }
    }
    return a;
}

double* PoissonSolver::SetY(int Nx, int Ny, double* v){
    double Y[(Nx-2)*(Ny-2)];
    y = &Y[0];
    int k = 0;
    for(int j=1; j<Ny-1; j++){
        for(int i=1; i<Nx-1; i++){
            *(y + k) = *(v + j*Nx + i);//omega[i][j];
            k ++;
        }
    }
    return y;
}

double* PoissonSolver::SetX(int Nx, int Ny){
    double X[(Nx-2)*(Ny-2)];
    x = &X[0];
    return x;
}

double* PoissonSolver::Execute(double Lx, double Ly, int Nx, int Ny, double* v, double* s){
    a = SetA(Lx, Ly, Nx, Ny);
    y = SetY(Nx, Ny, v);
    x = SetX(Nx, Ny);
    
    int n = (Nx-2)*(Ny-2);
    double R[n] = {0.0};
    double* r = &R[0];
    double P[n] = {0.0};
    double* p = &P[0];
    double T[n] = {0.0};
    double* t = &T[0];  //temp vector for r vector
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
    
    return x;
}
