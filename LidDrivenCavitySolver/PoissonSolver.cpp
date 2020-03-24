#include "PoissonSolver.h"
#include "cblas.h"

//CONSTRUCTORS
PoissonSolver::PoissonSolver(){
}


PoissonSolver::~PoissonSolver(){
}



//MEMBER FUNCTIONS

//create the coefficient matrix A for the global problem
double* PoissonSolver::SetGlobalA(double Lx, double Ly, int Nx, int Ny){
    const double deltaX = Lx/(Nx-1);
    const double deltaY = Ly/(Ny-1);
    const double X2 = deltaX*deltaX;
    const double Y2 = deltaY*deltaY;
    
    //create matrix A - symmetric upper triangle
    int Ax = Nx - 2;
    int Ay = Ny - 2;
    int Dim = Ax*Ay;
    
    A = new double[Dim*Dim];
    for(int j=0; j<Dim; j++){
        for(int i=0; i<Dim; i++){
            if(i==j){
                *(A + j*Dim + i) = 2*((1/X2)+(1/Y2));
                if((j*Dim + i + Ax) < (Dim*Dim - Ay*Dim)){
                    *(A + j*Dim + i + Ax) = -1/Y2;
                }
                if((j-(Ax-1))%Ax != 0){
                    *(A + j*Dim + i + 1) = -1/X2;
                }
            }
        }
    }
    
    return A;
}

//create the coefficient matrix a for the local problem
double* PoissonSolver::SetLocalA(int Nx, int Ny, int Px, int Py, int startCol, int endCol, int startRow, int endRow, int rank){
    int nx = Nx/Px;
    int ny = Ny/Py;
    
    int Ax = Nx - 2;
    int Ay = Ny - 2;
    int Dim = Ax*Ay;
    
    //relate local start and end points to the global
    int rankCol = rank%Px;
    int globalStartCol = (startCol - 1) + nx*rankCol;
    int globalEndCol   = globalStartCol + endCol - startCol;
    
    int rankRow = rank%Py;
    int globalStartRow = (startRow - 1) + ny*rankRow;
    int globalEndRow   = globalStartRow + endRow - startRow;
    
    int ax = endCol - startCol + 1;
    int ay = endRow - startRow + 1;
    int dim = ax*ay;
    
    //create matrix a (symmetric upper triangle) - coefficient matrix for the 'local' case
    int k = 0;
    
    a = new double[dim*dim];
    for(int R1 = globalStartRow; R1 <= globalEndRow; R1++){
        for(int C1 = globalStartCol; C1 <= globalEndCol; C1++){
            for(int R2 = globalStartRow; R2 <= globalEndRow; R2++){
                for(int C2 = globalStartCol; C2 <= globalEndCol; C2++){
                    *(a + k) = *(A + Dim*(C1-1 + Ay*(R1-1)) + (C2-1 + Ax*(R2-1)));
                    k++;
                }
            }
        }
    }
    return a;
}

double* PoissonSolver::SetY(int Nx, int Ny, int Px, int Py, int startCol, int endCol, int startRow, int endRow, double* v){
    int nx = Nx/Px;
    int ny = Ny/Py;
    int augX = nx + 2;
    
    y = new double[(endCol - startCol + 1)*(endRow - startRow + 1)];
    int k = 0;
    for(int j=startRow; j<=endRow; j++){
        for(int i=startCol; i<=endCol; i++){
            *(y + k) = *(v + j*augX + i);//omega[i][j];
            k ++;
        }
    }
    return y;
}

double* PoissonSolver::SetX(int Nx, int Ny, int Px, int Py, int startCol, int endCol, int startRow, int endRow){
    int nx = Nx/Px;
    int ny = Ny/Py;
    
    x = new double[(endCol - startCol + 1)*(endRow - startRow + 1)];
    return x;
}


double* PoissonSolver::Execute(double Lx, double Ly, int Nx, int Ny, int Px, int Py, int startCol, int endCol, int startRow, int endRow, double* v, double* s, int rank){
    SetGlobalA(Lx, Ly, Nx, Ny);
    SetLocalA(Nx, Ny, Px, Py, startCol, endCol, startRow, endRow, rank);
    SetY(Nx, Ny, Px, Py, startCol, endCol, startRow, endRow, v);
    SetX(Nx, Ny, Px, Py, startCol, endCol, startRow, endRow);
    
    int nx = Nx/Px;
    int ny = Ny/Py;
    
    int n = (endCol - startCol + 1)*(endRow - startRow + 1);
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
    double tol = 0.0001;
    
    //setup things
    cblas_dcopy(n, y, 1, r, 1); //r_0 = b (i.e. y)
    cblas_dsymv(CblasRowMajor, CblasUpper, n, -1.0, a, n, x, 1, 1.0, r, 1); //r_0 = y - Ax_0
    cblas_dcopy(n, r, 1, p, 1); //p_0 = r_0

    //loop
    do{
        cblas_dsymv(CblasRowMajor, CblasUpper, n, 1.0, a, n, p, 1, 0.0, t, 1); //t = A*p_k
        
        alpha = cblas_ddot(n, t, 1, p, 1);          //alpha = trans(p_k) * A * p_k
        
        /*for(int i=0; i<8; i++){
            *(x + i) = alpha;//*(t + i);
        }*/
        alpha = cblas_ddot(n, r, 1, r, 1) / alpha;  //alpha_k = trans(r_k) * r_k / trans(p_k) * A * p_k
        beta  = cblas_ddot(n, r, 1, r, 1);          //trans(r_k) * r_k
        
        cblas_daxpy(n,  alpha, p, 1, x, 1);         //x_k+1 = x_k + alpha_k * p_k
        cblas_daxpy(n, -alpha, t, 1, r, 1);         //r_k+1 = r_k - alpha_k * A * p_k
        
        eps = cblas_dnrm2(n, r, 1);                 //euclidean norm
        if(eps < tol*tol){
            break;
        }
        
        beta = cblas_ddot(n, r, 1, r, 1) / beta;    //trans(r_k+1) * r_k+1 / trans(r_k) * r_k
        
        cblas_dcopy(n, r, 1, t, 1);                 //copy r into temp, t
        cblas_daxpy(n, beta, p, 1, t, 1);           //P_k+1 = r_k+1 + beta * P_k
        cblas_dcopy(n, t, 1, p, 1);                 //copy temp, t into p
        
        k++;
    } while(k<5000);
    
    /*for(int i=0; i<16; i++){
        *(x + i) = *(a + 3*(Nx-2)*(Ny-2) + i);
    }*/

    return x;
}
