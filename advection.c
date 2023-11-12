/* This is a sample Advection solver in C 
The advection equation-> \partial q / \partial t - u \cdot \nabla q(x,y) = 0
The grid of NX by NX evenly spaced points are used for discretization.  
The first and last points in each direction are boundary points. 
Approximating the advection operator by 1st order finite difference. 
*/
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include "advection.h"
#include <string.h>
#include <stdbool.h>

#define BUFSIZE 512
/* ************************************************************************** */
int main ( int argc, char *argv[] ){
  if(argc!=2){
    printf("Usage: ./levelSet input.dat\n");
    return -1;  
  }
  static int frame=0;

  // Create an advection solver
  solver_t advc; 
  // Create uniform rectangular (Cartesian) mesh
  advc.msh = createMesh(argv[1]); 
  // Create time stepper 
  tstep_t tstep = createTimeStepper(advc.msh.Nnodes); 
  // Create Initial Field
  initialCondition(&advc);

  // Read input file for time variables 
  tstep.tstart = readInputFile(argv[1], "TSART");
  tstep.tend   = readInputFile(argv[1], "TEND");
  tstep.dt     = readInputFile(argv[1], "DT");
  tstep.time = 0.0; 

  // adjust time step size 
  int Nsteps = ceil( (tstep.tend - tstep.tstart)/tstep.dt);
  tstep.dt = (tstep.tend - tstep.tstart)/Nsteps;

  // Read input file for OUTPUT FREQUENCY i.e. in every 1000 steps
  int Noutput = readInputFile(argv[1], "OUTPUT_FREQUENCY");


  // write the initial solution i.e. q at t = tstart
  {
    char fname[BUFSIZ];
    sprintf(fname, "test_%04d.csv", frame++);
    solverPlot(fname, &advc.msh, advc.q);
  }


  // ********************Time integration***************************************/
  // for every steps
  for(int step = 0; step<Nsteps; step++){
    // for every stage
    for(int stage=0; stage<tstep.Nstage; stage++){
      // Call integration function
      RhsQ(&advc, &tstep, stage); 
    }

    tstep.time = tstep.time+tstep.dt;

    if(step%Noutput == 0){
      char fname[BUFSIZ];
      sprintf(fname, "test_%04d.csv", frame++);
      solverPlot(fname, &advc.msh, advc.q);
    }
  }
}

/* ************************************************************************** */
void RhsQ(solver_t *solver, tstep_t *tstep, int stage){

mesh_t *msh = &solver->msh;
for(int j=0; j<msh->NY; j++){
    for(int i=0; i<msh->NX; i++){
      /*
      fill this part to compute rhsq at every n = j*msh.NX + i
      */


      /* 
      Time integration in 2 steps 
      Step 1: Update residual
      resq = rk4a(stage)* resq + dt*rhsq
      Step:2 Update solution and store
       q = q + rk4b(stage)*resq
      */

      // Store updated solutions in solver->q and solver->resq 

    }
  }
}

/* ************************************************************************** */
void initialCondition(solver_t *solver){
  mesh_t *msh = &(solver->msh); 

  solver->q = (double *)malloc(msh->Nnodes*sizeof(double)); 
  solver->u = (double *)malloc(2*msh->Nnodes*sizeof(double));

  for(int j=0; j<msh->NY; j++){
    for(int i=0; i<msh->NX; i++){
     /*
     Create initial condition and velocity field
     */

    }
  }

}



/* ************************************************************************** */
// void createMesh(struct mesh *msh){
mesh_t createMesh(char* inputFile){

  mesh_t msh; 

  // Read required fields i.e. NX, NY, XMIN, XMAX, YMIN, YMAX

  msh.NX   = readInputFile(inputFile, "NX");
  msh.NY   = readInputFile(inputFile, "NY");
  msh.xmin   = readInputFile(inputFile, "XMIN");
  msh.xmax   = readInputFile(inputFile, "XMAX");
  msh.ymin   = readInputFile(inputFile, "YMIN");
  msh.ymax   = readInputFile(inputFile, "YMAX");
  /* 
  Continue with other required fields DONE!
  */


  msh.Nnodes = msh.NX*msh.NY;
  msh.x = (double *) malloc(msh.Nnodes*sizeof(double));
  msh.y = (double *) malloc(msh.Nnodes*sizeof(double));

  /*
  Compute Coordinates of the nodes
  */

  double x_interval=(msh.xmax-msh.xmin)/(msh.NX-1); //horizontal interval between two nodes 
  double y_interval=(msh.ymax-msh.ymin)/(msh.NY-1); //vertical interval between two nodes 

  double grid_array[msh.NX*msh.NY][2]; //creating array for calculating positions of nodes , x-->row , y--> column 
  int n;

  for(int j=0; j<msh.NY; j++){
    for(int i=0; i<msh.NX; i++){
     /*
      Complete this part
    */
    n=j*msh.NX+i; //calculation of n to get rid of (i,j) representation
    grid_array[n][0]=i*x_interval; //x coordinates of nodes
    grid_array[n][1]=j*y_interval; //y coordinates of nodes 
    }
  }

  //printf("x:%f\n",grid_array[(401*401)-1][1]); //to check the result correctness 

  // Create connectivity and periodic connectivity
  /* 
  for every node 4 connections east north west and south
  Nothe that periodic connections require specific treatment
  */
  msh.N2N = (int *)malloc(4*msh.Nnodes*sizeof(int)); // yukarı-->0, sağ-->1, aşağı-->2, sol-->3 

  for(int j=0; j<msh.NY; j++){
    for(int i=0; i<msh.NX; i++){
    /*
     Complete this part
    */
    n=j*msh.NX+i;
      if(i==0 && j!=0 && j!=msh.NY-1){
        //left side
        msh.N2N[4*n]=n+msh.NX; 
        msh.N2N[4*n+1]=n+1;
        msh.N2N[4*n+2]=n-msh.NX; 
        msh.N2N[4*n+3]=n+(n-1)-1; 
      }else if (j==0 && i!=0 && i!=msh.NX-1){
        //bottom side
        msh.N2N[4*n]=n+msh.NX; 
        msh.N2N[4*n+1]=n+1;
        msh.N2N[4*n+2]=((msh.NY-1)*4)+i; 
        msh.N2N[4*n+3]=n-1; 
      }else if (i==msh.NX-1 && j!=0 && j!=msh.NY-1){
        //right side 
        msh.N2N[4*n]=n+msh.NX; 
        msh.N2N[4*n+1]=(n-msh.NX)+1;
        msh.N2N[4*n+2]=n-msh.NX; 
        msh.N2N[4*n+3]=n-1; 
      }else if (j==msh.NY-1 && i!=0 && i!=msh.NX-1){
        //upper side 
        msh.N2N[4*n]=i; 
        msh.N2N[4*n+1]=n+1;
        msh.N2N[4*n+2]=n-msh.NX; 
        msh.N2N[4*n+3]=n-1; 
      }else if(i==0&&jj==0){
        //bottom left corner 
        msh.N2N[4*n]=n+msh.NX; 
        msh.N2N[4*n+1]=n+1;
        msh.N2N[4*n+2]=((msh.NY-1)*4); 
        msh.N2N[4*n+3]=(n+(msh.NX))-1; 
      }else if (i==msh.NX-1 && j==0){
        //bottom right corner 
        msh.N2N[4*n]=n+msh.NX; 
        msh.N2N[4*n+1]=(n-msh.NX)+1;
        msh.N2N[4*n+2]=((msh.NY-1)*4)+(msh.NX-1); 
        msh.N2N[4*n+3]=n-1; 
      }else if (j==msh.NY-1 && i==0){
        //top left corner 
        msh.N2N[4*n]=0; 
        msh.N2N[4*n+1]=n+1;
        msh.N2N[4*n+2]=n-msh.NX; 
        msh.N2N[4*n+3]=(msh.NX*msh.NY)-1;
      }else if (j==msh.NY-1 && i==msh.NX-1){
        //top rigth corner 
        msh.N2N[4*n]=i; 
        msh.N2N[4*n+1]=(n-msh.NX)+1;
        msh.N2N[4*n+2]=n-msh.NX; 
        msh.N2N[4*n+3]=n-1;
      } else {
        //inner
        msh.N2N[4*n]=n+msh.NX; 
        msh.N2N[4*n+1]=n+1;
        msh.N2N[4*n+2]=n-(msh.NX); 
        msh.N2N[4*n+3]=n-1;
      }
    }
  }

  return msh; 
}

/* ************************************************************************** */
void solverPlot(char *fileName, mesh_t *msh, double *Q){
    FILE *fp = fopen(fileName, "w");
    if (fp == NULL) {
        printf("Error opening file\n");
        return;
    }

    fprintf(fp, "X,Y,Z,Q \n");
    for(int n=0; n< msh->Nnodes; n++){
      fprintf(fp, "%.8f, %.8f,%.8f,%.8f\n", msh->x[n], msh->y[n], 0.0, Q[n]);
    } 
}

/* ************************************************************************** */
double readInputFile(char *fileName, char* tag){
  FILE *fp = fopen(fileName, "r");
  if (fp == NULL) {
    printf("Error opening the input file\n");
    return -1;
  }

  size_t original_size=strlen(tag); //size of original string
  char str[original_size+3]; //new string to copy old one 
  snprintf(str,sizeof(str),"[%s]",tag); //creating a new string with brackets 
  //printf("str:%s\n",str); 
  size_t size=strlen(str); //length of final string 

  char line[5000]; //to read line by line 
  bool found=false; 
  char *conversion_variable; //to convert it double at the end 
  double variable; //to hold variable in double type 

  while (fgets(line, sizeof(line), fp)) { //reading line by line 
 
    char str_2[size]; //new string to copy read line 
    if(found==true){
      variable=strtod(line,&conversion_variable); //conversion to double 
      //printf("next read line:%s",line);
      found=false;
    }
    strncpy(str_2,line,size); //holding new string 
    if(strcmp(str,str_2)==0){
      //printf("current read line:%s",line); 
      found=true; 
    }
  }
  //printf("Variable:%f\n",variable); 

  fclose(fp);
  return variable; 
  /* 
  Complete this function to read the file for the given tag
  */
}


/* ************************************************************************** */
// Time stepper clas RK(4-5)
// resq = rk4a(stage)* resq + dt*rhsq
//  q = q + rk4b(stage)*resq
tstep_t createTimeStepper(int Nnodes){
  tstep_t tstep; 
  tstep.Nstage = 5; 
  tstep.resq = (double *)calloc(Nnodes,sizeof(double)); 
  tstep.rhsq = (double *)calloc(Nnodes,sizeof(double));
  tstep.rk4a = (double *)malloc(tstep.Nstage*sizeof(double));
  tstep.rk4b = (double *)malloc(tstep.Nstage*sizeof(double));
  tstep.rk4c = (double *)malloc(tstep.Nstage*sizeof(double));

  tstep.rk4a[0] = 0.0; 
  tstep.rk4a[1] = -567301805773.0/1357537059087.0; 
  tstep.rk4a[2] = -2404267990393.0/2016746695238.0;
  tstep.rk4a[3] = -3550918686646.0/2091501179385.0;
  tstep.rk4a[4] = -1275806237668.0/842570457699.0;
        
  tstep.rk4b[0] = 1432997174477.0/9575080441755.0;
  tstep.rk4b[1] = 5161836677717.0/13612068292357.0; 
  tstep.rk4b[2] = 1720146321549.0/2090206949498.0;
  tstep.rk4b[3] = 3134564353537.0/4481467310338.0;
  tstep.rk4b[4] = 2277821191437.0/14882151754819.0;
             
  tstep.rk4c[0] = 0.0;
  tstep.rk4c[1] = 1432997174477.0/9575080441755.0;
  tstep.rk4c[2] = 2526269341429.0/6820363962896.0;
  tstep.rk4c[3] = 2006345519317.0/3224310063776.0;
  tstep.rk4c[4] = 2802321613138.0/2924317926251.0;
  return tstep; 
}


