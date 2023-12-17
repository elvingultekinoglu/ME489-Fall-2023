# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <time.h>
#include <omp.h>
# define OUT 0




// Function definitions
int main ( int argc, char *argv[] );
double boundary_condition ( double x, double time );
double initial_condition ( double x, double time );
double source ( double x, double time );
void runSolver( int n, int size );



/*-------------------------------------------------------------
  Purpose: Compute number of primes from 1 to N with naive way
 -------------------------------------------------------------*/
// This function is fully implemented for you!!!!!!
// usage: mpirun -n 4 heat1d N
// N    : Number of nodes per processor
int main ( int argc, char *argv[] ){
  int rank, size;

  // get number of nodes per processor
  int N = strtol(argv[1], NULL, 10);
  //printf("number of nodes per core %d\n",N);

  size = strtol(argv[2], NULL, 10);
  //printf("number of core %d\n",size);

  omp_set_num_threads(size);

  // Solve and update the solution in time
  runSolver(N, size);


  return 0;
}

/*-------------------------------------------------------------
  Purpose: computes the solution of the heat equation.
 -------------------------------------------------------------*/
void runSolver( int n, int size ){
  // CFL Condition is fixed
  double cfl = 0.5; 
  // Domain boundaries are fixed
  double x_min=0.0, x_max=1.0;
  // Diffusion coefficient is fixed
  double k   = 0.002;
  // Start time and end time are fixed
  double tstart = 0.0, tend = 10.0;  


  double *x, *q, *qn;

  x  = ( double*)malloc(((n*size))*sizeof(double));
  q  = ( double*)malloc(((n*size))*sizeof(double));
  qn = ( double*)malloc(((n*size))*sizeof(double));
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
  double time, time_new, wtime;  

  // find the coordinates for uniform spacing 
  for ( int i = 0; i <= (size*n) - 1; i++ ){
    // COMPLETE THIS PART
    // x[i] = ....
    x[i]=  i*dx; 
    //printf("rank:%d, i:%d, position:%f \n",rank,i,x[i]); // control position

  }

  time = tstart; 
  for (int i = 0; i <= size*n-1; i++ ){
    q[i] = initial_condition(x[i],time);
  }



     
  // Compute the values of H at the next time, based on current data.
  double start_time = omp_get_wtime();

 
  for ( int step = 1; step <= Nsteps; step++ ){
    time_new = time + step*dt;
    double rhs;
    // UPDATE the solution based on central differantiation.
    // qn[i] = q[i] + dt*rhs(q,t)
    // For OpenMP make this loop parallel also
    int i;
#pragma omp parallel default(none) private(i, rhs) shared(size, n, q, k, x, dx, qn, dt, time)
{
    #pragma omp for
    for (i = 1; i <= (size * n) - 2; i++) {
        // COMPLETE THIS PART
        rhs = k * (q[i - 1] - 2 * q[i] + q[i + 1] + source(x[i], time)) / (dx * dx);
        qn[i] = q[i] + dt * rhs;
    }
}

    // q at the extreme left and right boundaries was incorrectly computed
    // using the differential equation.  
    // Replace that calculation by the boundary conditions.
    // global left endpoint 
      qn[0] = boundary_condition ( x[0], time_new );
    // global right endpoint 

      qn[size*n-1] = boundary_condition ( x[size*n-1], time_new );
  // Update time and field.
    time = time_new;
    // For OpenMP make this loop parallel also
#pragma omp parallel default(none) private(i) shared(n, size, q, qn)
{
    #pragma omp for
    for (i = 0; i <= n * size - 1; i++) {
        q[i] = qn[i];
    }
}


  }


      // Stop the timer
    double end_time = omp_get_wtime();

    // Calculate and print the elapsed time
    double elapsed_time = end_time - start_time;
    printf("Elapsed time of parallel OpenMp solver: %f seconds\n", elapsed_time);


    // write out the x coordinates for display.
    xfile = fopen ( "openmp_x_data.txt", "w" );
    for (int i = 0; i<(size*n); i++ ){
      fprintf ( xfile, "  %f", x[i] );
    }
    fprintf ( xfile, "\n" );
    fclose ( xfile );
    // write out the initial solution for display.
    qfile = fopen ( "openmp_q_data.txt", "w" );
     for (int i = 0; i<(size*n); i++ ){
      fprintf ( qfile, "  %f", q[i] );
    }
    fprintf ( qfile, "\n" );


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
