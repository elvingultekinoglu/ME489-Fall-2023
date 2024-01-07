# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <string.h>

#define BUFSIZE 512

// Include MPI header
# include "mpi.h"


// Function definitions
/******************************************************************************/
int main ( int argc, char *argv[] );
double exactSoln( double c, double x, double y, double t );
void applyBC(double *data,  double *x, double *y, double c, double time, int nx, int ny,int rank, int size);
void solverPlot(char *fileName, double *x, double *y, int nx, int ny, double *data); 
double readInputFile(char *fileName, char* tag); 

// Solver Info
/******************************************************************************
  Purpose:
    wave2d solves the wave equation in parallel using MPI.
  Discussion:
    Discretize the equation for u(x,t):
      d^2 u/dt^2  =  c^2 * (d^2 u/dx^2 + d^2 u/dy^2)  
      for 0 < x < 1, 0 < y < 1, t>0
    with boundary conditions and Initial conditions obtained from the exact solutions:
      u(x,y, t) = sin ( 2 * pi * ( x - c * t ) )
   Usage: serial -> ./wave input.dat  parallel> mpirun -np 4 ./wave input.dat 
******************************************************************************/

int main ( int argc, char *argv[] ){
  
  // Read input file for solution parameters
  double tstart = readInputFile(argv[1], "TSART"); // Start time
  double tend   = readInputFile(argv[1], "TEND");  // End time
  double dt     = readInputFile(argv[1], "DT");    // Time step size

  // Global node number in x and y
  int NX        = (int) readInputFile(argv[1], "NX"); // Global node numbers in x direction
  int NY        = (int) readInputFile(argv[1], "NY"); // Global node numbers in y direction

  double xmax = readInputFile(argv[1], "XMAX"); // domain boundaries
  double xmin = readInputFile(argv[1], "XMIN"); // domain boundaries
  double ymax = readInputFile(argv[1], "YMAX"); // domain boundaries
  double ymin = readInputFile(argv[1], "YMIN"); // domain boundaries
  double c = readInputFile(argv[1], "CONSTANT_C");



  // Paralelization start 
  int rank, size;
  // Initialize MPI, get size and rank
  MPI_Init ( &argc, &argv );
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &size );

  // Record start time
  double startTime = MPI_Wtime();


  double *qn, *q0, *q1;               // Solution field at t+dt, t and t-dt
  static int frame=0; 
 
  // DOMAIN DECOMPOSITION
  // For serial implementation nx = NX and ny = NY; 
  int nx = NX;      // local number of nodes in x direction
  int ny = NY/size;      // local number of nodes in y direction

  // ALLOCATE MEMORY for COORDINATES (x, y) and compute them
  double *x = ( double * ) malloc ( nx * (ny+2) * sizeof ( double ) );
  double *y = ( double * ) malloc ( nx * (ny+2) * sizeof ( double ) );
  // find uniform spacing in x and y directions
  double hx = (xmax - xmin)/(NX-1.0); 
  double hy = (ymax - ymin)/(NY-1.0); 
  // Compute coordinates of the nodes
  for(int j=0; j < ny+2; ++j){
    for(int i=0; i < nx;++i){
      x[i+j*nx]= ( i*hx); 
      y[i+j*nx]= -hy + ( j*hy) + (rank * hy * ny) ;
    }
  }


  // ALLOCATE MEMORY for SOLUTION and its HISTORY
  // Solution at time (t+dt)
  qn = ( double * ) malloc ( nx * (ny+2) * sizeof ( double ) );
  // Solution at time (t)
  q0 = ( double * ) malloc ( nx * (ny+2) * sizeof ( double ) );
  // Solution at time t-dt
  q1 = ( double * ) malloc ( nx * (ny+2) * sizeof ( double ) );

  // USE EXACT SOLUTION TO FILL HISTORY
   for(int i=0; i<nx; i++){
      for(int j=0; j<ny+2; j++){
      const double xn = x[i+ j*nx]; 
      const double yn = y[i+ j*nx]; 
      // Exact solutions at history tstart and tstart+dt
      q0[i + j*nx] = exactSoln(c, xn, yn, tstart + dt);  
      q1[i + j*nx] = exactSoln(c, xn, yn, tstart);  
    }
  }

 
  // Write the initial solution 
  {
    char fname[BUFSIZ];
    sprintf(fname, "Core_%d_test_%04d.csv", rank,frame++);
    solverPlot(fname, x, y, nx, ny, q1);
  }

// RUN SOLVER 
  int Noutput = 10000; 
  int Nsteps=(tend - tstart)/dt;     // Assume  dt divides (tend- tstart)
  double alphax2 = pow((c*dt/hx),2); 
  double alphay2 = pow((c*dt/hy),2);
  
  // We already have 2 steps computed with exact solution
  double time = dt; 


  // for every time step
  for(int tstep = 2; tstep<=Nsteps+1; ++tstep){
  //for(int tstep = 2; tstep<=3; ++tstep){
    // increase  time
    time = tstart + tstep*dt; 
    
    // Apply Boundary Conditions i.e. at i, j = 0, i,j = nx-1, ny-1
    applyBC(q0, x, y, c, time, nx, ny,rank,size);

//Communication Part 
if(size>1){

if (rank != size - 1) {
  MPI_Request send_upper_request[2];
    MPI_Isend(&q0[nx*ny], nx, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD, &send_upper_request[0]);
    MPI_Isend(&q1[nx*ny], nx, MPI_DOUBLE, rank + 1, rank+nx, MPI_COMM_WORLD, &send_upper_request[1]);
}
if (rank != 0) {
    //double rec;
    MPI_Request recv_lower_request[2];
    MPI_Irecv(&q0[0], nx, MPI_DOUBLE, rank - 1, rank-1, MPI_COMM_WORLD, &recv_lower_request[0]);
    MPI_Irecv(&q1[0], nx, MPI_DOUBLE, rank - 1, rank-1 + nx, MPI_COMM_WORLD, &recv_lower_request[1]);
    MPI_Waitall(2, recv_lower_request, MPI_STATUSES_IGNORE);
}
if (rank != 0) {
  MPI_Request send_lower_request[2];
    MPI_Isend(&q0[nx], nx, MPI_DOUBLE, rank - 1, rank+nx+nx-1, MPI_COMM_WORLD, &send_lower_request[0]);
    MPI_Isend(&q1[nx], nx, MPI_DOUBLE, rank - 1, rank+nx+nx+nx-1, MPI_COMM_WORLD, &send_lower_request[1]);
}
if (rank != size - 1) {
    MPI_Request recv_upper_request[2];
    MPI_Irecv(&q0[nx*(ny+1)], nx, MPI_DOUBLE, rank + 1, rank+nx+nx, MPI_COMM_WORLD, &recv_upper_request[0]);
    MPI_Irecv(&q1[nx*(ny+1)], nx, MPI_DOUBLE, rank + 1, rank+nx+nx+nx, MPI_COMM_WORLD, &recv_upper_request[1]);
    MPI_Waitall(2, recv_upper_request, MPI_STATUSES_IGNORE);
}

}




// Update solution using second order central differencing in time and space
    
if(rank==0){
    for(int i=1; i<nx-1; i++){ // exclude left right boundaries
      for(int j= 2; j<ny+1 ; j++){ // exclude top and bottom boundaries
        const int n0   = i + j*nx; 
        const int nim1 = i - 1 + j*nx; // node i-1,j
        const int nip1 = i + 1 + j*nx; // node i+1,j
        const int njm1 = i + (j-1)*nx; // node i, j-1
        const int njp1 = i + (j+1)*nx; // node i, j+1
        // update solution 
        qn[n0] = 2.0*q0[n0] - q1[n0] + alphax2*(q0[nip1]- 2.0*q0[n0] + q0[nim1])
                                     + alphay2*(q0[njp1] -2.0*q0[n0] + q0[njm1]); 
      }
    }
}else if(rank==(size-1)) {
    for(int i=1; i<nx-1; i++){ // exclude left right boundaries
      for(int j= 1; j<ny ; j++){ // exclude top and bottom boundaries
        const int n0   = i + j*nx; 
        const int nim1 = i - 1 + j*nx; // node i-1,j
        const int nip1 = i + 1 + j*nx; // node i+1,j
        const int njm1 = i + (j-1)*nx; // node i, j-1
        const int njp1 = i + (j+1)*nx; // node i, j+1
        // update solution 
        qn[n0] = 2.0*q0[n0] - q1[n0] + alphax2*(q0[nip1]- 2.0*q0[n0] + q0[nim1])
                                     + alphay2*(q0[njp1] -2.0*q0[n0] + q0[njm1]); 
      }
    }
}else {
    for(int i=1; i<nx-1; i++){ // exclude left right boundaries
      for(int j= 1; j<ny+1 ; j++){ // exclude top and bottom boundaries
        const int n0   = i + j*nx; 
        const int nim1 = i - 1 + j*nx; // node i-1,j
        const int nip1 = i + 1 + j*nx; // node i+1,j
        const int njm1 = i + (j-1)*nx; // node i, j-1
        const int njp1 = i + (j+1)*nx; // node i, j+1
        // update solution 
        qn[n0] = 2.0*q0[n0] - q1[n0] + alphax2*(q0[nip1]- 2.0*q0[n0] + q0[nim1])
                                     + alphay2*(q0[njp1] -2.0*q0[n0] + q0[njm1]); 
      }
    }
}

// Update history q1 = q0; q0 = qn, except the boundaries
    if(rank==0){
    for(int i=1; i<nx-1; i++){
      for(int j=2; j<ny+1; j++){
        q1[i + j*nx] = q0[i + j*nx]; 
        q0[i + j*nx] = qn[i + j*nx]; 
      }
    }
  } else if (rank==(size-1)) {
    for(int i=1; i<nx-1; i++){
      for(int j=1; j<ny; j++){
        q1[i + j*nx] = q0[i + j*nx]; 
        q0[i + j*nx] = qn[i + j*nx]; 
      }
    }
  } else {
    for(int i=1; i<nx-1; i++){
      for(int j=1; j<ny+1; j++){
        q1[i + j*nx] = q0[i + j*nx]; 
        q0[i + j*nx] = qn[i + j*nx]; 
      }
    }
  }
    // Dampout a csv file for postprocessing
    if(tstep%Noutput == 0){
    //if(time==1){ 
      char fname[BUFSIZ];
      sprintf(fname, "Core_%d_test_%04d.csv", rank,frame++);
      solverPlot(fname, x, y, nx, ny, q0);
    }
  }

  // Compute Linf norm of error at tend
    double linf = 0.0; 
    //int xx= 0;
    //int yy= 0;
    for(int i=0; i<nx; i++){
      for(int j=1; j<ny+1; j++){
         double xn = x[i+ j*nx]; 
         double yn = y[i+ j*nx]; 
         // solution and the exact one
         double qn = q0[i+ j*nx]; 
         double qe = exactSoln(c, xn, yn, time); 
         //if (fabs(qn-qe)>linf) {
         //xx  = i;
         //yy  = j;
         //} 
         linf  = fabs(qn-qe)>linf ? fabs(qn -qe):linf; 
      }
    }

    // Compute biggest Linf norm of error at tend
    double globalLinf = 0.0;
    MPI_Reduce(&linf, &globalLinf, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Record end time
    double endTime = MPI_Wtime();

    // Only print the result on rank 0
    if (rank == 0) {
    double executionTime = endTime - startTime;
    printf("Global Infinity norm of the error: %.4e, total calculation time: %.4f, Total core number: %d\n", globalLinf,executionTime,size);
    }

  // Terminate MPI.
  MPI_Finalize ( );
free(x);
free(y);
free(qn);
free(q0);
free(q1);

  return 0;
}




/***************************************************************************************/
double exactSoln( double c, double x, double y, double t){
  const double pi = 3.141592653589793;
  double value = sin( 2.0*pi*( x - c*t));
  return value;
}

/***************************************************************************************/
void applyBC(double *data,  double *x, double *y, double c, double time, int nx, int ny,int rank, int size){

  // Apply Boundary Conditions
  double xn, yn; 

  for(int j=0; j<ny+2;++j){ // left right boundaries i.e. i=0 and i=nx-1
    xn = x[0 + j*nx]; 
    yn = y[0 + j*nx];    
    data[0 + j*nx] = exactSoln(c, xn, yn, time); 

    xn = x[nx-1 + j*nx]; 
    yn = y[nx-1 + j*nx];    
    data[nx-1 + j*nx] = exactSoln(c, xn, yn, time); 
  }

  if (rank ==(size-1)){
  for(int i=0; i< nx; ++i){ // top and  bottom boundaries i.e. j=0 and j=ny-1
    xn = x[i+ (ny)*nx]; 
    yn = y[i+ (ny)*nx];       
    data[i +  (ny)*nx] = exactSoln(c, xn, yn, time);
  }
}

  if(rank==0){
  for(int i=0; i< nx; ++i){ // top and  bottom boundaries i.e. j=0 and j=ny-1
    xn = x[i+ 1*nx]; 
    yn = y[i+ 1*nx]; 
    data[i + 1*nx] = exactSoln(c, xn, yn, time); 
  }
}
}

/* ************************************************************************** */
void solverPlot(char *fileName, double *x, double *y, int nx, int ny, double *Q){
    FILE *fp = fopen(fileName, "w");
    if (fp == NULL) {
        printf("Error opening file\n");
        return;
    }

    fprintf(fp, "X,Y,Z,Q \n");
     for(int i=0; i<nx; i++){
      for(int j=1; j<ny+1; j++){ // not printing ghost lines
        const double xn = x[i + j*nx]; 
        const double yn = y[i + j*nx]; 
        fprintf(fp, "%.8f, %.8f,%.8f,%.8f\n", xn, yn, 0.0, Q[i + j*nx]);
      }
    }
}


/* ************************************************************************** */
double readInputFile(char *fileName, char* tag){
  FILE *fp = fopen(fileName, "r");
  if (fp == NULL) {
    printf("Error opening the input file\n");
    return -1;
  }

  int sk = 0; 
  double result; 
  char buffer[BUFSIZE];
  char fileTag[BUFSIZE]; 
  while(fgets(buffer, BUFSIZE, fp) != NULL){
    sscanf(buffer, "%s", fileTag);
    if(strstr(fileTag, tag)){
      fgets(buffer, BUFSIZE, fp);
      sscanf(buffer, "%lf", &result); 
      return result;
    }
    sk++;
  }

  if(sk==0){
    printf("could not find the tag: %s in the file %s\n", tag, fileName);
  }
}