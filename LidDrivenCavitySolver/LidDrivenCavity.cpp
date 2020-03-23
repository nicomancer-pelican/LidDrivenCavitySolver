#include "LidDrivenCavity.h"
#include <mpi.h>

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

void LidDrivenCavity::SetPartitions(int px, int py){
    Px = px;
    Py = py;
}

void LidDrivenCavity::SetRank(int Rank){
    rank = Rank;
}



//other member functions
void LidDrivenCavity::Initialise(){
    //set size local to each processor
    nx = Nx/Px;
    ny = Ny/Py;
    
    //set augmented size to account for guard cells
    augX = nx + 2;
    augY = ny + 2;
    
    //pointers to matrices
    v = new double[augX*augY];
    s = new double[augX*augY];
    
    //mark the start of the 'internal' section - assume processor is in the middle somewhere and not subject to wall conditions
    startRow = 1;
    startCol = 1;
    endRow = ny;
    endCol = nx;
}

void LidDrivenCavity::Integrate(){
    PoissonSolver* poisson = new PoissonSolver();
    poisson->SetGlobalA(Lx, Ly, Nx, Ny);
    poisson->SetLocalA(Nx, Ny, Px, Py, startCol, endCol, startRow, endRow, rank);
    poisson->SetY(Nx, Ny, Px, Py, startCol, endCol, startRow, endRow, v);
    //poisson->InitialisePoisson(Nx, Ny, Lx, Ly, Px, Py, startCol, endCol, startRow, endRow, v, rank);
    
    //for(double t=0.0; t<T; t+=dt){
        guardCells();
        boundaryConditions();
        interiorV();
        guardCells();
        newInteriorV();
        //updateS(poisson->Execute(Lx, Ly, Nx, Ny, v, s));
    //}
}


//GUARD CELLS
void LidDrivenCavity::guardCells(){
    //setup Cartesian topology - 2x2 grid, 1 processor for each grid
    MPI_Comm cartesianGrid;
    const int ndims = 2;           //number of dimentions of cartesian array
    const int dims[2] = {Py, Px};  //integer array of size ndims with number of processors
    const int periods[2] = {0, 0}; //logical array of size ndims specifygin if grid is periodic
    
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 0, &cartesianGrid);
    
    //find the rank of processes on each side
    int shift[4] = {rank, rank, rank, rank};                    //{left, right, up, down}
    MPI_Cart_shift(cartesianGrid, 1, 1, &shift[0], &shift[1]);  //left and right
    MPI_Cart_shift(cartesianGrid, 0, 1, &shift[3], &shift[2]);  //up and down
    
    //buffer vectors
    double* BuffCol = new double[ny];
    double* BuffRow = new double[nx];
    
    //left guard cells
    if(shift[0] < 0){
        startCol = 2;
    } else{
        //fill buffer to send
        for(int j=1; j<augY-1; j++){
            *(BuffCol + j - 1) = *(v + j*augX + 1);
        }
        MPI_Sendrecv_replace(BuffCol, ny, MPI_DOUBLE, shift[0], rank, shift[0], rank-1, cartesianGrid, MPI_STATUS_IGNORE);
        //put buffer in matrix
        for(int j=1; j<augY-1; j++){
            *(v + j*augX) = *(BuffCol + j - 1);
        }
    }
    
    //right guard cells
    if(shift[1] < 0){
        endCol = nx - 1;
    } else{
        //fill buffer to send
        for(int j=1; j<augY-1; j++){
            *(BuffCol + j - 1) = *(v + j*augX + nx);
        }
        MPI_Sendrecv_replace(BuffCol, ny, MPI_DOUBLE, shift[1], rank, shift[1], rank+1, cartesianGrid, MPI_STATUS_IGNORE);
        //put buffer in matrix
        for(int j=1; j<augY-1; j++){
            *(v+ j*augX + nx + 1) = *(BuffCol + j - 1);
        }
    }
    
    //up guard cells
    if(shift[2] < 0){
        endRow = ny - 1;
    } else{
        //fill buffer to send
        for(int i=1; i<augX-1; i++){
            *(BuffRow + i - 1) = *(v + augX*ny + i);
        }
        MPI_Sendrecv_replace(BuffRow, nx, MPI_DOUBLE, shift[2], rank, shift[2], rank+Px, cartesianGrid, MPI_STATUS_IGNORE);
        //put buffer in matrix
        for(int i=1; i<augX-1; i++){
            *(v + augX*(ny+1) + i) = *(BuffRow + i - 1);
        }
    }
    
    //down guard cells
    if(shift[3] < 0){
        startRow = 2;
    } else{
        //fill buffer to send
        for(int i=1; i<augX-1; i++){
            *(BuffRow + i - 1) = *(v + augX + i);
        }
        MPI_Sendrecv_replace(BuffRow, nx, MPI_DOUBLE, shift[3], rank, shift[3], rank-Px, cartesianGrid, MPI_STATUS_IGNORE);
        //put buffer in matrix
        for(int i=1; i<augX-1; i++){
            *(v + i) = *(BuffRow + i - 1);
        }
    }
    
    
    delete BuffCol; BuffCol = nullptr;
    delete BuffRow; BuffRow = nullptr;
}


//STEP 1 MEMBER FUNCTIONS
void LidDrivenCavity::boundaryConditions(){
    const double deltaX = Lx / (Nx-1);
    const double deltaY = Ly / (Ny-1);
    const double X2     = deltaX * deltaX;
    const double Y2     = deltaY * deltaY;
    
    //top wall condition
    if(endRow == ny-1){
        for(int i=1; i<augX-1; i++){
            *(v + ny*augX + i) = (*(s + augX*ny + i) - *(s + (ny-1)*augX + i))*(2/Y2) - (2/deltaY);
        }
    }
    
    //bottom wall condition
    if(startRow == 2){
        for(int i=1; i<augX-1; i++){
            *(v + augX + i) = (*(s + augX + i) - *(s + 2*augX + i))*(2/Y2);
        }
    }
    
    //left wall condition
    if(startCol == 2){
        for(int j=1; j<augY-1; j++){
            *(v + j*augX + 1) = (*(s + j*augX + 1) - *(s + j*augX + 2))*(2/X2);
        }
    }
    
    //right wall condition
    if(endCol == nx-1){
        for(int j=1; j<augY-1; j++){
            *(v + augX*(j+1) - 2) = (*(s + augX*(j+1) - 2) - *(s + augX*(j+1) - 3))*(2/X2);
        }
    }
}

//STEP 2 MEMBER FUNCTION
void LidDrivenCavity::interiorV(){
    const double deltaX = Lx / (Nx-1);
    const double deltaY = Ly / (Ny-1);
    const double X2     = deltaX * deltaX;
    const double Y2     = deltaY * deltaY;
    for(int j=startRow; j<=endRow; j++){
        for(int i=startCol; i<endCol; i++){
            *(v + j*augX + i) = -((*(s + j*augX + i + 1) - *(s + j*augX + i)*2 + *(s + j*augX + i - 1)) / X2)
                              -((*(s + (j+1)*augX + i) - *(s + j*augX + i)*2 + *(s + (j-1)*augX + i)) / Y2);
        }
    }
}

//Setp 3 MEMBER FUNCTION
void LidDrivenCavity::newInteriorV(){
    const double deltaX = Lx / (Nx-1);
    const double deltaY = Ly / (Ny-1);
    const double X2     = deltaX * deltaX;
    const double Y2     = deltaY * deltaY;
    for(int j=startRow; j<=endRow; j++){
        for(int i=startCol; i<=endCol; i++){
            double a = ((*(v + j*augX + i + 1) - *(v + j*augX + i)*2 + *(v + j*augX + i -1))/X2) + ((*(v + (j+1)*augX + i) - *(v + j*augX + i)*2 + *(v + (j-1)*augX + i))/Y2);
            double b = ((*(v + (j+1)*augX + i) - *(v + (j-1)*augX + i)) / (2*deltaY)) * ((*(v + j*augX + i + 1) - *(v + j*augX + i - 1)) / (2*deltaX));
            double c = ((*(v + j*augX + i + 1) - *(v + j*augX + i - 1)) / (2*deltaX)) * ((*(v + (j+1)*augX + i) - *(v + (j-1)*augY + i)) / (2*deltaY));
            *(v + j*augX + i) = ((1/Re)*a - b + c)*dt + *(v + j*augX + i);
            
            /**(v + j*augX + i) = *(v + j*augX + i) + dt*(
                               (1/Re) * (((*(v + j*augX + i + 1) - *(v + j*augX + i)*2 + *(v + j*augX + i - 1)) / X2)
                                        +((*(v + augX*(j+1) + i) - *(v + j*augX + i)*2 + *(v + augX*(j-1) + i)) / Y2))
                              -(((*(s + augX*(j+1) + i) - *(s + augX*(j-1) + i)) / (2*deltaY))
                               *((*(v + j*augX + i + 1) - *(v + j*augX + i - 1)) / (2 *deltaX)))
                              +(((*(s + j*augX + i + 1) - *(s + j*augX + i - 1)) / (2*deltaX))
                               *((*(v + augX*(j+1) + i) - *(v + augX*(j-1) + i)) / (2*deltaY)))
                                );*/
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


