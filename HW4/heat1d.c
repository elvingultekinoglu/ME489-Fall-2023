# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# define OUT 0

// Include MPI header
# include "mpi.h"

// Function definitions
int main ( int argc, char *argv[] );
double boundary_condition ( double x, double time );
double initial_condition ( double x, double time );
double source ( double x, double time );
void runSolver( int n, int rank, int size );



/*-------------------------------------------------------------
  Purpose: Compute number of primes from 1 to N with naive way
 -------------------------------------------------------------*/
// This function is fully implemented for you!!!!!!
// usage: mpirun -n 4 heat1d N
// N    : Number of nodes per processor
int main ( int argc, char *argv[] ){
  int rank, size;
  double wtime;

  // Initialize MPI, get size and rank
  MPI_Init ( &argc, &argv );
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &size );

  // get number of nodes per processor
  int N = strtol(argv[1], NULL, 10);


  // Solve and update the solution in time
  runSolver(N, rank, size);

  // Terminate MPI.
  MPI_Finalize ( );
  // Terminate.
  return 0;
}

/*-------------------------------------------------------------
  Purpose: computes the solution of the heat equation.
 -------------------------------------------------------------*/
void runSolver( int n, int rank, int size ){
  // CFL Condition is fixed
  double cfl = 0.5; 
  // Domain boundaries are fixed
  double x_min=0.0, x_max=1.0;
  // Diffusion coefficient is fixed
  double k   = 0.002;
  // Start time and end time are fixed
  double tstart = 0.0, tend = 10.0;  

  // Storage for node coordinates, solution field and next time level values
  double *x, *q, *qn;
  // Set the x coordinates of the n nodes padded with +2 ghost nodes. 
  x  = ( double*)malloc((n+2)*sizeof(double));
  q  = ( double*)malloc((n+2)*sizeof(double));
  qn = ( double*)malloc((n+2)*sizeof(double));

  // Write solution field to text file if size==1 only
  FILE *qfile, *xfile;

  // uniform grid spacing
  double dx = ( x_max - x_min ) / ( double ) ( size * n - 1 );

  // Set time step size dt <= CFL*h^2/k
  // and then modify dt to get integer number of steps for tend
  double dt  = cfl*dx*dx/k; 
  int Nsteps = ceil(( tend - tstart )/dt);
  dt =  ( tend - tstart )/(( double )(Nsteps)); 
    //printf("nsteps: %d dt: %f \n ",Nsteps, dt);
    //fflush(stdout);

  int tag;
  MPI_Status status;
  double time, time_new, wtime;  

  // find the coordinates for uniform spacing 
  for ( int i = 0; i <= n + 1; i++ ){
    // COMPLETE THIS PART
    // x[i] = ....
    x[i]= -dx + ( i*dx) + (rank * dx * n) ; 
    //printf("rank:%d, i:%d, position:%f \n",rank,i,x[i]); // control position

  }

  // Set the values of q at the initial time.
  time = tstart; q[0] = 0.0; q[n+1] = 0.0;
  for (int i = 1; i <= n; i++ ){
    q[i] = initial_condition(x[i],time);
  }


 // Record the starting time.
  wtime = MPI_Wtime();

     
  // Compute the values of H at the next time, based on current data.
  for ( int step = 1; step <= Nsteps; step++ ){

    time_new = time + step*dt; 

    
// Perform point to point communications here!!!!

if(size>1){

if (rank != size - 1) {
  MPI_Request send_right_request;
    MPI_Isend(&q[n], 1, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD, &send_right_request);
    //printf("rank : %d,send to %d, massage %f \n",rank,rank+1,q[n]);
    //fflush(stdout);
}

if (rank != 0) {
    //double rec;
    MPI_Request recv_left_request;
    MPI_Irecv(&q[0], 1, MPI_DOUBLE, rank - 1, rank-1, MPI_COMM_WORLD, &recv_left_request);
    MPI_Wait(&recv_left_request, MPI_STATUS_IGNORE); // Wait for the receive operation to complete
    //printf("rank : %d,recito %d massage %f \n",rank,rank-1,rec);
    //fflush(stdout);
    //q[0]=rec;

}

if (rank != 0) {
  MPI_Request send_left_request;
    MPI_Isend(&q[1], 1, MPI_DOUBLE, rank - 1, rank+5, MPI_COMM_WORLD, &send_left_request);
    //printf("rank : %d,send to %d \n",rank,rank+1);
    //printf("aaa %d",rank);
}

if (rank != size - 1) {

    MPI_Request recv_right_request;
    MPI_Irecv(&q[n + 1], 1, MPI_DOUBLE, rank + 1, rank+1+5, MPI_COMM_WORLD, &recv_right_request);
    MPI_Wait(&recv_right_request, MPI_STATUS_IGNORE); // Wait for the receive operation to complete
    //printf("rank : %d,recito %d \n",rank,rank-1);
}

/*    // Wait for the completion of the first round of communication
    MPI_Wait(&send_right_request, MPI_STATUS_IGNORE);
    MPI_Wait(&recv_left_request, &recv_left_status);

    // Wait for the completion of the second round of communication
    MPI_Wait(&send_left_request, MPI_STATUS_IGNORE);
    MPI_Wait(&recv_right_request, &recv_right_status);
*/
}

    double rhs;
    // UPDATE the solution based on central differantiation.
    // qn[i] = q[i] + dt*rhs(q,t)
    // For OpenMP make this loop parallel also
    for ( int i = 1; i <= n; i++ ){
      // COMPLETE THIS PART
      rhs =  k * (q[i-1] - 2*q[i] + q[i+1] + source(x[i],time)) / (dx*dx); 
      qn[i] = q[i] + dt*rhs;
    }

  
    // q at the extreme left and right boundaries was incorrectly computed
    // using the differential equation.  
    // Replace that calculation by the boundary conditions.
    // global left endpoint 
    if (rank==0){
      qn[1] = boundary_condition ( x[1], time_new );
    }
    // global right endpoint 
    if (rank == size - 1 ){
      qn[n] = boundary_condition ( x[n], time_new );
    }


  // Update time and field.
    time = time_new;
    // For OpenMP make this loop parallel also
    for ( int i = 1; i <= n; i++ ){
      q[i] = qn[i];
    }

  // In single processor mode, add current solution data to output file.
    if (size == 1 && OUT==1){
      for ( int i = 1; i <= n; i++ ){
        fprintf ( qfile, "  %f", q[i] );
      }
      fprintf ( qfile, "\n" );
    }




 if (step==Nsteps){
    char x_filename[40];
    char q_filename[40];
    sprintf(x_filename, "mpi_x_data_size%d_rank%d.txt",size, rank);
    sprintf(q_filename, "mpi_q_data_size%d_rank%d.txt",size, rank);
    // write out the x coordinates for display.
    xfile = fopen ( x_filename, "w" );
    for (int i = 1; i<(n+1); i++ ){
      fprintf ( xfile, "  %f", x[i] );
    }
    fprintf ( xfile, "\n" );
    fclose ( xfile );
    // write out the initial solution for display.
    qfile = fopen ( q_filename, "w" );
    for ( int i = 1; i < (n+1); i++ ){
      fprintf ( qfile, "  %f", q[i] );
    }
    fprintf ( qfile, "\n" );
    fclose(qfile);
  }  
  /*
      if (step==1){
    char x_filename[40];
    char q_filename[40];
    sprintf(x_filename, "x_data_rank%d_step%d.txt", rank,step);
    sprintf(q_filename, "q_data_rank%d_step%d.txt", rank,step);
    // write out the x coordinates for display.
    xfile = fopen ( x_filename, "w" );
    for (int i = 1; i<(n+1); i++ ){
      fprintf ( xfile, "  %f", x[i] );
    }
    fprintf ( xfile, "\n" );
    fclose ( xfile );
    // write out the initial solution for display.
    qfile = fopen ( q_filename, "w" );
    for ( int i = 1; i < (n+1); i++ ){
      fprintf ( qfile, "  %f", q[i] );
    }
    fprintf ( qfile, "\n" );
    fclose(qfile);
  }
 */


  }

 /*
    char x_filename[40];
    char q_filename[40];
    sprintf(x_filename, "mpi_x_data_rank%d_step%d.txt", rank,10);
    sprintf(q_filename, "mpi_q_data_rank%d_step%d.txt", rank,10);
    // write out the x coordinates for display.
    xfile = fopen ( x_filename, "w" );
    for (int i = 1; i<(n+1); i++ ){
      fprintf ( xfile, "  %f", x[i] );
    }
    fprintf ( xfile, "\n" );
    fclose ( xfile );
    // write out the initial solution for display.
    qfile = fopen ( q_filename, "w" );
    for ( int i = 1; i < (n+1); i++ ){
      fprintf ( qfile, "  %f", q[i] );
    }
    fprintf ( qfile, "\n" );
    fclose(qfile);
  */
  
  // Record the final time.
  // if (rank == 0 ){
  wtime = MPI_Wtime( )-wtime;

  // Add local number of primes
  double global_time = 0.0; 
  MPI_Reduce( &wtime, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

 if(rank==0)
   printf ( " SİZE: %d RANK: %d Wall clock elapsed seconds = %f\n",size, rank,global_time );      


  if( size == 1 && OUT==1)
    fclose ( qfile ); 

  free(q); free(qn); free(x);

  return;
}
/*-----------------------------------------------------------*/
double boundary_condition ( double x, double time ){
  double value;

  // Left condition:
  if ( x < 0.5 ){
    value = 100.0 + 10.0 * sin ( time );
  }else{
    value = 75.0;
  }
  return value;
}
/*-----------------------------------------------------------*/
double initial_condition ( double x, double time ){
  double value;
  value = 95.0;

  return value;
}
/*-----------------------------------------------------------*/
double source ( double x, double time ){
  double value;

  value = 0.0;

  return value;
}
