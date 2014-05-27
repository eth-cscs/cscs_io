/*
// Test for parallel I/O benchmark
// Written by Jean M. Favre, Swiss National Supercomputer Center
// Last update: Wed Jun 8 2011
//
// tested and compiled on palu and rosa with the following modules
module load PrgEnv-gnu
module load netcdf-hdf5parallel/4.1.1.0
module load adios
module load mxml
module load python/2.6.4
//
*/
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/times.h>
#include <stdlib.h>
#include <mpi.h>
#include <zlib.h>

int get_procmem(double *bytes) 
{
  FILE *fh;
  int proc_ret;
  char proc_var[80];
  char *cp;
  long long int ibytes;
 
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

  
  *bytes=0.0;
  fh = fopen("/proc/self/status","r");
  while(!feof(fh)) {
    fgets(proc_var,80,fh);
    cp = strstr(proc_var,"VmHWM:");
    if (cp) {sscanf(cp, "VmHWM:"" %llu",&ibytes );
      *bytes=max(*bytes,ibytes);
    }
  }
  fclose(fh);
  *bytes *= 1024.0;
   return 0;
}

#define ZONAL_DATA 1

typedef float mytype;

//#define GX 256 // size of MPI grid in X 
//#define GY 128 // size of MPI grid in Y 
const int NFIELDS = 10;   // number of arrays of scalar variables (all float)
int BlockSize = 64;
//int LONGITUDE_SIZE = 64;    // number of cells in X direction
//int LATITUDE_SIZE = 64;     // number of cells in Y direction
//int LEVELS_SIZE = 64;         // number of cells in Z direction
int NUMBER_OF_ITERATIONS = 1;
int calculate = 0; // set to 0 to skip all computation and do I/O only; otherwise 1
int gzipped = 0; // set to 1 to gzip the file while writing; otherwise 0

typedef struct datagrid
{
    int    Nrows;     /* local size of data array */
    int    Ncolumns;  /* local size of data array */
    int    Nlevels;   /* local size of data array */
    float *rmesh_x;
    float *rmesh_y;
    float *rmesh_z;
    int    rmesh_dims[3];
    int    rmesh_ndims;
    mytype **data;
} datagrid;

typedef struct simulation_data
{
    int       cycle;
    double    time;
    int       runMode;
    int       updateplots;
    int       done;
    int       savingFiles;
    int       saveCounter;
    int       par_rank;
    int       par_size;
    int       global_dims[3];
    int       coords[2]; // MPI topological grid of size 2
    int       cart_dims[2]; // MPI topological grid of size 2
    char      filename[256];
    MPI_Comm  comm_cart;
    datagrid grid;
} simulation_data;

#ifdef ADIOS
#include "adios.h"
#include "adios_types.h"
#include "ADIOS_WriteData.c"

#elif NETCDF4

#ifdef PARIO
#include <netcdf_par.h>
#include <netcdf.h>
#else
#include <netcdf.h>
#endif

#include "NETCDF4_WriteData.c"

#elif HDF5
#include <hdf5.h>
#include "HDF5_WriteData.c"
#endif

struct tms mpitms;
double mpi_elapsed;

#define SIM_STOPPED       0
#define SIM_RUNNING       1

#ifdef VISIT_INSITU
// used in in-situ visualization only
#define SHOW_GHOST_ARRAY 1
// number of grid points would be 2+ in the X and Y direction

#include <VisItControlInterface_V2.h>
#include <VisItDataInterface_V2.h>
#define PARALLEL 1
#include "SimulationExample.h"

/* Data Access Function prototypes */
visit_handle SimGetMetaData(void *);
visit_handle SimGetMesh(int, const char *, void *);
visit_handle SimGetVariable(int, const char *, void *);
visit_handle SimGetDomainList(const char *, void *);
void ControlCommandCallback(const char *cmd, const char *args, void *cbdata);
void SlaveProcessCallback();
static int visit_broadcast_int_callback(int *value, int sender, void *cbdata);
static int visit_broadcast_string_callback(char *str, int len, int sender, void *cbdata);
#endif

void MPIIO_WriteData(simulation_data *sim, char *Filename);
void BOV_WriteData(simulation_data *sim, char *Filename, int tindex, const char *postfix);

void
grid_data_ctor(int par_size, datagrid *grid)
{
  int f;
  grid->Ncolumns = BlockSize; //LONGITUDE_SIZE;
  grid->Nrows    = BlockSize; //LATITUDE_SIZE;
  grid->Nlevels  = BlockSize; //LEVELS_SIZE;
  grid->rmesh_x = NULL;
  grid->rmesh_y = NULL;
  grid->rmesh_z = NULL;
  grid->rmesh_dims[0] = 1; /* shall be redefined later after Parallel init*/
  grid->rmesh_dims[1] = 1;
  grid->rmesh_dims[2] = 1;
  grid->rmesh_ndims = 3;

  grid->data = (float**)malloc(sizeof(float*) * NFIELDS);
  for(f=0; f < NFIELDS; f++)
    {
    grid->data[f] = NULL;
    }
}

void grid_data_dtor(datagrid *grid)
{
  int f;
  if(grid->rmesh_x != NULL)
    {
    free(grid->rmesh_x);
    grid->rmesh_x = NULL;
    }
  if(grid->rmesh_y != NULL)
    {
    free(grid->rmesh_y);
    grid->rmesh_y = NULL;
    }
  if(grid->rmesh_z != NULL)
    {
    free(grid->rmesh_z);
    grid->rmesh_z = NULL;
    }
  for(f=0; f < NFIELDS; f++)
    {
    if(grid->data[f] != NULL)
      {
      free(grid->data[f]);
      grid->data[f] = NULL;
      }
    }
    free(grid->data);
    grid->data = NULL;
}

/*
A 3D Rectilinear mesh of size sim->global_dims[0] * sim->global_dims[1] * sim->global_dims[2] is built and is partitioned among MPI tasks
*/
void grid_data_allocate(simulation_data *sim)
{
    int i, size;
    datagrid *grid = &(sim->grid);   
    int *coords = sim->coords;
    float offset, zone_width;

#ifdef SHOW_GHOST_ARRAY
    grid->Ncolumns += 2; // add one cell layer to the east and to the west
    grid->Nrows += 2;   // add one cell layer to the south and to the north
#endif

    grid->rmesh_dims[0] = grid->Ncolumns+1; // size of X local coordinate array
    grid->rmesh_dims[1] = grid->Nrows+1;    // size of Y local coordinate array
    grid->rmesh_dims[2] = grid->Nlevels+1;  // size of Z local coordinate array

    grid->rmesh_x = (float*)malloc(sizeof(float) * grid->rmesh_dims[0]);
    grid->rmesh_y = (float*)malloc(sizeof(float) * grid->rmesh_dims[1]);
    grid->rmesh_z = (float*)malloc(sizeof(float) * grid->rmesh_dims[2]);

    //offset = coords[0] * (grid->rmesh_dims[0]-2);
#ifdef SHOW_GHOST_ARRAY
    offset = (grid->Ncolumns-2)*(coords[0])-1; // substract two degrees
#else
    offset = grid->Ncolumns*coords[0];
#endif
    zone_width = 1.0/sim->global_dims[0];
    for(i=0; i < grid->rmesh_dims[0]; i++)
    {
        grid->rmesh_x[i] = (offset + i)*zone_width;
    }

    //offset = coords[1] * (grid->rmesh_dims[1]-2);
#ifdef SHOW_GHOST_ARRAY
    offset = (grid->Nrows-2)*(coords[1])-1; // substract two degrees
#else
    offset = grid->Nrows*coords[1];
#endif
    zone_width = 1.0/sim->global_dims[1];
    for(i=0; i < grid->rmesh_dims[1]; i++)
    {
        grid->rmesh_y[i] = (offset + i)*zone_width; // was -90
    }
    if((coords[1] == 0) && (coords[0] >=0 && coords[0] <=2)){
      for(i=0; i < grid->rmesh_dims[0]; i++)
        {
        //fprintf(stderr,"%f ", grid->rmesh_x[i]);
        }
      //fprintf(stderr,"\n");
    }
    zone_width = 1.0/sim->global_dims[2];
    for(i=0; i < grid->rmesh_dims[2]; i++)
    {
        grid->rmesh_z[i] = i * zone_width;
    }

#ifdef ZONAL_DATA
    /* 3D array of data items exposed as solution. Data centered at cells*/
    size = sizeof(float) * (grid->Ncolumns) * (grid->Nrows) * (grid->Nlevels);
#else
    size = sizeof(float) * (grid->Ncolumns+1) * (grid->Nrows+1) * (grid->Nlevels+1);
#endif
    for(i=0; i < NFIELDS; i++)
      {
      grid->data[i] = (float*)malloc(size);
      }

    //fprintf(stderr,"p=[%d*%d], cell array alloc [%d*%d*%d]\n", coords[0],coords[1],grid->Ncolumns, grid->Nrows, grid->Nlevels);
}

void
data_simulate(datagrid *grid, int par_rank, int par_size, double T)
{
  int nsum, f, i, j, k;
  int *data = NULL;

#ifdef ZONAL_DATA
  for(k=0; k< grid->Nlevels; k++)
    {
    for(j=0; j< grid->Nrows; j++)
      {
      for(i=0; i < grid->Ncolumns; i++)
        {
        int idx = k*(grid->Ncolumns)*(grid->Nrows) + j*(grid->Ncolumns) + i;
        grid->data[0][idx] = cos((grid->rmesh_y[j] + T)*3.1415) * sin((grid->rmesh_x[i] + T)*3.1415);
/*
-90 < Y < 90
-pi/2 < arg < pi/2
 0 < cos() < 1
-1 < sin() < 1
*/
/*
        grid->data[1][idx] = cos((grid->rmesh_y[j]/180.0 + T)*3.1415) * sin((grid->rmesh_x[i]/360.0 + T)*3.1415); 
        grid->data[2][idx] = cos(grid->rmesh_y[j]*3.1415/360.0);
        grid->data[3][idx] = sin(grid->rmesh_x[i]*3.1415/360.0);
*/
        for(f=1; f < NFIELDS; f++)
          {
          grid->data[f][idx] = grid->rmesh_x[i] * grid->rmesh_x[i] - grid->rmesh_y[j] * grid->rmesh_y[j] - 
                                                                   grid->rmesh_x[i] * grid->rmesh_y[j] * grid->rmesh_z[k];;
          }
        }
      }
    }
#else
  for(k=0; k< grid->Nlevels+1; k++)
    {
    for(j=0; j< grid->Nrows+1; j++)
      {
      for(i=0; i < grid->Ncolumns+1; i++)
        {
        int idx = k*(grid->Ncolumns+1)*(grid->Nrows+1) + j*(grid->Ncolumns+1) + i;
        grid->data[0][idx] = grid->rmesh_x[i] * grid->rmesh_x[i] - grid->rmesh_y[j] * grid->rmesh_y[j] - 
                                                                   grid->rmesh_x[i] * grid->rmesh_y[j] * grid->rmesh_z[k];
/*
-90 < Y < 90
-pi/2 < arg < pi/2
 0 < cos() < 1
-1 < sin() < 1
*/
/*
        grid->data[1][idx] = cos((grid->rmesh_y[j]/180.0 + T)*3.1415) * sin((grid->rmesh_x[i]/360.0 + T)*3.1415); 
        grid->data[2][idx] = cos(grid->rmesh_y[j]*3.1415/360.0);
        grid->data[3][idx] = sin(grid->rmesh_x[i]*3.1415/360.0);
*/
        for(f=1; f < NFIELDS; f++)
          {
          grid->data[f][idx] = cos((grid->rmesh_y[j] + T)*3.1415) * sin((grid->rmesh_x[i] + T)*3.1415);
          }
        }
      }
    }
#endif

}

void
simulation_data_ctor(simulation_data *sim)
{
    sim->updateplots = 1;
    sim->cycle = 0;
    sim->time = 0.;
    sim->done = 0;
    sim->savingFiles = 0;
    sim->saveCounter = 0;
    sim->par_rank = 0;
    sim->par_size = 1;
    strcpy(sim->filename, "null");
    grid_data_ctor(sim->par_size, &sim->grid);
}

void
simulation_data_dtor(simulation_data *sim)
{
    grid_data_dtor(&(sim->grid));
}

const char *cmd_names[] = {"halt", "step", "run", "Update ON/OFF", "Save ON/OFF"};

void BOV_WriteData(simulation_data *sim, char *Filename, int tindex, const char *postfix)
{
  //int dimuids[3]={sim->global_dims[2], sim->global_dims[1], sim->global_dims[0]};
  int i, rank, rc, ucount;
  FILE *fp;
  gzFile *gp;
  // given a rank, we need the index in the cartesian grid by flipping the indices
  char lname[256];
  sprintf(lname, "%s.%05d.%04d.%s", Filename, sim->cart_dims[0]*sim->coords[1]+sim->coords[0], tindex, postfix);
  ucount = sim->grid.Ncolumns * sim->grid.Nrows * sim->grid.Nlevels;

  if(gzipped) {
    gp = gzopen(lname, "wb");
    gzwrite(gp, (const void *)sim->grid.data[0], sizeof(float) * ucount);
    gzclose(gp);
  }
  else {
    fp = fopen(lname, "wb");
    fwrite(sim->grid.data[0], sizeof(float), ucount, fp);
    fclose(fp);
  }


}

void MPIIO_WriteData(simulation_data *sim, char *Filename)
{
  int dimuids[3]={sim->global_dims[2], sim->global_dims[1], sim->global_dims[0]};
  int f, rc, ustart[3], ucount[3];
  MPI_Offset disp = 0;
  long offset;
  MPI_File      filehandle;
  MPI_Datatype  filetype;

  rc = MPI_File_open(sim->comm_cart, Filename, 
		MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &filehandle);

  ustart[2] = sim->grid.Ncolumns * sim->coords[0];
  ustart[1] = sim->grid.Nrows * sim->coords[1];
  ustart[0] = 0;
  ucount[2] = sim->grid.Ncolumns;
  ucount[1] = sim->grid.Nrows;
  ucount[0] = sim->grid.Nlevels;

 // Create the subarray representing the local block
  MPI_Type_create_subarray(3, dimuids, ucount, ustart,
					MPI_ORDER_C, MPI_FLOAT, &filetype);
  MPI_Type_commit(&filetype);

	
  for(f=0; f < NFIELDS; f++)
    {
    MPI_File_set_view(filehandle, disp, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);

    MPI_File_write_all(filehandle, sim->grid.data[f],
                       ucount[0]*ucount[1]*ucount[2],
                       MPI_FLOAT, MPI_STATUS_IGNORE);
    disp += sim->global_dims[2] * sim->global_dims[1] * sim->global_dims[0] * sizeof(float);
    }
  MPI_File_close(&filehandle);
  MPI_Type_free(&filetype);
}

void DumpData(simulation_data *sim, char *Filename, int tindex, const char *postfix)
{
  char lname[512];
#ifdef ADIOS
  if(sim->cycle == 1){
    //fprintf(stderr,"Calling adios_init()\n");
    //adios_init ("/users/jfavre/Projects/ADIOS/benchmark.xml");
    //adios_init ("/users/jfavre/Projects/ADIOS/benchmark.xml", sim->comm_cart);
    adios_init ("benchmark.xml", sim->comm_cart);
    }
  sprintf(lname,"%s.%04d.%s", Filename, tindex, postfix);
  ADIOS_WriteData(sim, lname);
  if(sim->cycle == NUMBER_OF_ITERATIONS){
    //fprintf(stderr, "Calling adios_finalize()\n");
    adios_finalize (sim->par_rank);
    }
#elif NETCDF4
  sprintf(lname,"%s.%04d.%s", Filename, tindex, postfix);
  NETCDF4_WriteData(sim, Filename);
#elif HDF5
  sprintf(lname,"%s.%04d.%s", Filename, tindex, postfix);
  HDF5_WriteData(sim, lname);
#elif BOV
  if(gzipped)
    BOV_WriteData(sim, Filename, tindex, "bof.gz");
  else
    BOV_WriteData(sim, Filename, tindex, "bof");
  if(sim->par_rank == 0)
    {
    sprintf(lname,"%s.%04d.bov", Filename, tindex);
    FILE *fp = fopen(lname,"w");
    fprintf(fp,"# BOV version: 1.0\n");
    fprintf(fp,"# file written by IO benchmark program\n");
    if(gzipped) 
      fprintf(fp,"DATA_FILE: %s.%%05d.%04d.bof.gz\n", Filename, tindex);
    else
      fprintf(fp,"DATA_FILE: %s.%%05d.%04d.bof\n", Filename, tindex);
    fprintf(fp,"DATA SIZE: %d %d %d\n", sim->global_dims[0], sim->global_dims[1], sim->global_dims[2]);
    fprintf(fp,"DATA_BRICKLETS: %d %d %d\n", sim->grid.Nrows, sim->grid.Ncolumns, sim->grid.Nlevels);
    fprintf(fp,"DATA FORMAT: FLOAT\n");
    fprintf(fp,"VARIABLE: node_data\n");
    //fprintf(fp,"VARIABLE PALETTE MIN: 0\n");
    //fprintf(fp,"VARIABLE PALETTE MAX: 14.7273\n");
    fprintf(fp,"BRICK ORIGIN: 0.0 0.0 0.0\n");
    fprintf(fp,"BRICK SIZE: %f %f %f\n", 1.0*sim->global_dims[0], 1.0*sim->global_dims[1], 1.0*sim->global_dims[2]);
    fprintf(fp,"BRICK X_AXIS: 1.000 0.000 0.000\n");
    fprintf(fp,"BRICK Y_AXIS: 0.000 1.000 0.000\n");
    fprintf(fp,"BRICK Z_AXIS: 0.000 0.000 1.000\n");
    fprintf(fp,"DATA_ENDIAN: LITTLE\n");
    fprintf(fp,"CENTERING: nodal\n");
    fprintf(fp,"BYTE_OFFSET: 0\n");
    fclose(fp);
    }
#else
  sprintf(lname,"%s.%04d.%s", Filename, tindex, postfix);
  MPIIO_WriteData(sim, lname);
#endif
}

void
simulate_one_timestep(simulation_data *sim)
{
  ++sim->cycle;
    if(sim->cycle == NUMBER_OF_ITERATIONS)
      {
      sim->done = 1;
      }
  sim->time += .01;

    /* Simulate the current round of grid. */

  if(calculate)
    data_simulate(&sim->grid, sim->par_rank, sim->par_size, sim->time);
#ifdef VISIT_INSITU
    if(sim->updateplots)
       {
       VisItTimeStepChanged();
       VisItUpdatePlots();
       }
    if(sim->savingFiles)
      {
      char filename[100];
      sprintf(filename, "/scratch/eiger/jfavre/updateplots1k%04d.jpg", sim->saveCounter);
        if(VisItSaveWindow(filename, 128, 128, VISIT_IMAGEFORMAT_JPEG) == VISIT_OKAY)
        {
            sim->saveCounter++;
            if(sim->par_rank == 0) printf("Saved %s\n", filename);
        }
        else if(sim->par_rank == 0)
            printf("The image could not be saved to %s\n", filename);
    }

#endif
}

void mainloop(simulation_data *sim)
{
  int blocking, visitstate, err = 0;
  char Filename[256];

#ifdef VISIT_INSITU
  if(sim->runMode == SIM_STOPPED)
    simulate_one_timestep(sim);
  do
    {
    blocking = (sim->runMode == SIM_RUNNING) ? 0 : 1;

// Get input from VisIt or timeout so the simulation can run.
    if(sim->par_rank == 0)
      visitstate = VisItDetectInput(blocking, -1);
    MPI_Bcast(&visitstate, 1, MPI_INT, 0, sim->comm_cart);

// Do different things depending on the output from VisItDetectInput
    switch(visitstate)
      {
      case 0:
// There was no input from VisIt, return control to sim
        simulate_one_timestep(sim);
      break;
      case 1:
// VisIt is trying to connect to sim.
        if(VisItAttemptToCompleteConnection() == VISIT_OKAY)
          {
          //fprintf(stderr, "VisIt connected\n");
          VisItSetCommandCallback(ControlCommandCallback, (void*)sim);
          VisItSetSlaveProcessCallback2(SlaveProcessCallback, (void*)sim);

          VisItSetGetMetaData(SimGetMetaData, (void*)sim);
          VisItSetGetMesh(SimGetMesh, (void*)sim);
          VisItSetGetVariable(SimGetVariable, (void*)sim);
          VisItSetGetDomainList(SimGetDomainList, (void*)sim);
          }
        else
          {
          char *err = VisItGetLastError();
          fprintf(stderr, "VisIt did not connect: %s\n", err);
          free(err);
          }
      break;
      case 2:
// VisIt wants to tell the engine something
        if(!ProcessVisItCommand(sim))
          {
// Disconnect on an error or closed connection.
           VisItDisconnect();
// Start running again if VisIt closes.
           sim->runMode = SIM_RUNNING;
          }
      break;
      default:
        fprintf(stderr, "Can't recover from error %d!\n", visitstate);
        err = 1;
      break;
      }
    } while(!sim->done && err == 0);
#else
  do
    {
    simulate_one_timestep(sim);
    char postfix[4];
#ifdef ADIOS
    strcpy(postfix, "bp");
    DumpData(sim, sim->filename, sim->saveCounter++, postfix);
#elif NETCDF4
#ifdef PARIO
    strcpy(postfix, "p.nc");
    DumpData(sim, sim->filename, sim->saveCounter++, postfix);
#else
    strcpy(postfix, "s.nc");
    DumpData(sim, sim->filename, sim->saveCounter++, postfix);
#endif
#elif HDF5
    strcpy(postfix, "h5");
    DumpData(sim, sim->filename, sim->saveCounter++, postfix);
#elif BOV
    strcpy(postfix, "bof");
    char name[512];
    char cmd[512];
    int l0 =  strlen(strrchr(sim->filename, '/'));
    int l =  1+strlen(sim->filename) - l0;
    strncpy(name, sim->filename, l);
    sprintf(&name[l], "%05d/", sim->saveCounter);
    sprintf(cmd, "mkdir %s", &name[0]);
    if(!sim->par_rank){
      system(cmd);
    }
    MPI_Barrier(sim->comm_cart);
    strncpy(&name[l+6], &sim->filename[l], l0);
    //fprintf(stderr, "creating a subdirectory %s %d\n",name, l0);
    //fprintf(stderr, "creating a subdirectory %s %d\n",sim->filename, sim->saveCounter);
    DumpData(sim, name, sim->saveCounter++, postfix);
#else
    strcpy(postfix, "bin");
    DumpData(sim, sim->filename, sim->saveCounter++, postfix);
#endif

    //if(sim->par_rank == 0) printf("Saved file %s\n", Filename);

    if(sim->cycle == NUMBER_OF_ITERATIONS)
      {
      sim->done = 1;
      }
    } while(!sim->done && err == 0);
#endif
}

int main(int argc, char **argv)
{
  char *env = NULL;
  simulation_data sim;
  simulation_data_ctor(&sim);
  sim.runMode = SIM_STOPPED;
  int r, c, ndims=2, dims[2]={8,8}, periods[2]={0,0};
  int xdivs=2, ydivs=2, zdivs=1;
  double before, after, *mem;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &sim.par_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &sim.par_size);
  if(!sim.par_rank)
    get_procmem(&before);

#ifdef VISIT_INSITU
    SimulationArguments(argc, argv);
#endif
#define OPTS 1
#ifdef OPTS
  if(!sim.par_rank)
    {
    while(1)
      {
       static struct option long_options[] =
             {
               {"xdivs",  required_argument, 0, 'x'},
               {"ydivs",  required_argument, 0, 'y'},
               {"blocksize",  required_argument, 0, 'b'},
               {"timesteps",  required_argument, 0, 't'},
               {"file",    required_argument, 0, 'f'},
               {0, 0, 0, 0}
             };
           /* getopt_long stores the option index here. */
           int option_index = 0;
     
           c = getopt_long (argc, argv, "x:y:b:t:f:",
                            long_options, &option_index);
     
           /* Detect the end of the options. */
           if (c == -1)
             break;
     
           switch (c)
             {
             case 0:
               /* If this option set a flag, do nothing else now. */
               if (long_options[option_index].flag != 0)
                 break;
               printf ("option %s", long_options[option_index].name);
               if (optarg)
                 printf (" with arg %s", optarg);
               printf ("\n");
               break;
     
             case 'x':
               //printf ("option -x with value `%s'\n", optarg);
               dims[0] = atoi(optarg);
               break;
     
             case 'y':
               //printf ("option -y with value `%s'\n", optarg);
               dims[1] = atoi(optarg);
               break;

             case 'b':
               //printf ("option -b with value `%s'\n", optarg);
               BlockSize = atoi(optarg);
               break;

             case 't':
               //printf ("option -t with value `%s'\n", optarg);
               NUMBER_OF_ITERATIONS = atoi(optarg);
               break;
     
             case 'f':
               //printf ("option -f with value `%s'\n", optarg);
               strcpy(sim.filename, optarg);
               break;
     
             case '?':
               /* getopt_long already printed an error message. */
               break;
     
             default:
               abort ();
             }
         }
    }
#endif
  MPI_Bcast(dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&BlockSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&NUMBER_OF_ITERATIONS, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // broadcast name, even if it was not set
  MPI_Bcast(sim.filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);

  if(sim.par_size != dims[0]*dims[1])
    {
    if(sim.par_rank == 0)
      fprintf(stderr,"The cartesian topology requested (%d nodes) does not allow to map all tasks. Resubmit\n", sim.par_size);
    MPI_Finalize();
    exit(1);
    }

  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 1, &sim.comm_cart);
  MPI_Cart_coords(sim.comm_cart, sim.par_rank, 2, sim.coords);
  //sim.comm_cart = MPI_COMM_WORLD;
  sim.cart_dims[0] = dims[0];
  sim.cart_dims[1] = dims[1];

 /* Adjust the partitioning */
  sim.grid.Ncolumns = BlockSize; //LONGITUDE_SIZE; /* assume they divide evenly */
  sim.global_dims[0] = BlockSize * dims[0];

  sim.grid.Nrows = BlockSize; //LATITUDE_SIZE; /* assume they divide evenly */
  sim.global_dims[1] = BlockSize * dims[1];

  sim.global_dims[2] = BlockSize; //LEVELS_SIZE;
  sim.grid.Nlevels = BlockSize;

  //fprintf(stderr,"p=[%d], coords[%dx%d], array[%d*%d*%d] global[%d*%d*%d]\n", sim.par_rank, sim.coords[0], sim.coords[1], sim.grid.Ncolumns, sim.grid.Nrows, sim.grid.Nlevels, sim.global_dims[0], sim.global_dims[1], sim.global_dims[2]);

#ifdef VISIT_INSITU
    FILE *fp;
    if(sim.par_rank == 0){
        env = VisItGetEnvironment();
    //fp = fopen("/scratch/rosa/jfavre/foo.txt", "w");
    //fprintf(fp, "%s\n", env);
    }
      

    /* Pass the environment to all other processors collectively. */
    VisItSetupEnvironment2(env);
    if(env != NULL)
        free(env);

    /* Install callback functions for global communication. */
    VisItSetBroadcastIntFunction2(visit_broadcast_int_callback, (void*)&sim);
    VisItSetBroadcastStringFunction2(visit_broadcast_string_callback, (void*)&sim);
    /* Tell whether the simulation is parallel. */
    VisItSetParallel(sim.par_size > 1);
    VisItSetParallelRank(sim.par_rank);
    if(sim.par_rank == 0)
      {
      //fprintf(fp, "\n--reached before callbacks\n", env);
      //fclose(fp);
      }


    if(sim.par_rank == 0)
    {
        VisItInitializeSocketAndDumpSimFile("Cosmo",
            "Code YYY coupled with VisIt's libsim",
            "/path/to/where/sim/was/started",
            NULL, NULL, "/users/jfavre/.visit/simulations/rosa.sim2");
    }
#endif

  grid_data_allocate(&sim);

  mpi_elapsed = MPI_Wtime();
  mainloop(&sim);
  mpi_elapsed = MPI_Wtime() - mpi_elapsed;

  grid_data_dtor(&sim.grid);

  if(!sim.par_rank){
    get_procmem(&after);
    fprintf(stderr, "Testing get_procmem... %f %f %f\n",
          before, after, after-before);
    }

  MPI_Barrier(sim.comm_cart);
 
  if(sim.par_rank == 0)
    {
    printf("written grids of %d,%d,%d\n",  BlockSize, BlockSize, BlockSize);
    printf("written %d iterations\n", NUMBER_OF_ITERATIONS);
    printf("MPI Elapsed time: %f sec\n", mpi_elapsed);
    printf("average  %f Gbytes/sec\n", (BlockSize*BlockSize*BlockSize*4.0*sim.par_size*NUMBER_OF_ITERATIONS) / (1024*1024*1024*mpi_elapsed));
    }

  MPI_Finalize();

  return 0;
}


#ifdef VISIT_INSITU

visit_handle SimGetMetaData(void *cbdata)
{
    visit_handle md = VISIT_INVALID_HANDLE;
    simulation_data *sim = (simulation_data *)cbdata;
    int f;
    /* Create metadata. */
    if(VisIt_SimulationMetaData_alloc(&md) == VISIT_OKAY)
    {
        int i;
        visit_handle m1 = VISIT_INVALID_HANDLE;
        visit_handle vmd = VISIT_INVALID_HANDLE;
        visit_handle cmd = VISIT_INVALID_HANDLE;

        /* Set the simulation state. */
        VisIt_SimulationMetaData_setMode(md, (sim->runMode == SIM_STOPPED) ?
                VISIT_SIMMODE_STOPPED : VISIT_SIMMODE_RUNNING);
        VisIt_SimulationMetaData_setCycleTime(md, sim->cycle, sim->time);

        /* Set the first mesh's properties.*/
        if(VisIt_MeshMetaData_alloc(&m1) == VISIT_OKAY)
        {
            /* Set the mesh's properties.*/
            VisIt_MeshMetaData_setName(m1, "mesh");
            VisIt_MeshMetaData_setMeshType(m1, VISIT_MESHTYPE_RECTILINEAR);
            VisIt_MeshMetaData_setTopologicalDimension(m1, 3);
            VisIt_MeshMetaData_setSpatialDimension(m1, 3);
            VisIt_MeshMetaData_setNumDomains(m1, sim->par_size);
            VisIt_MeshMetaData_setXUnits(m1, "");
            VisIt_MeshMetaData_setYUnits(m1, "");
            VisIt_MeshMetaData_setZUnits(m1, "");
            VisIt_MeshMetaData_setXLabel(m1, "");
            VisIt_MeshMetaData_setYLabel(m1, "");
            VisIt_MeshMetaData_setZLabel(m1, "");
            VisIt_SimulationMetaData_addMesh(md, m1);
        }

        for(f=0; f < NFIELDS; f++)
          {
          char name[8];
          sprintf(name,"data%03d", f);
          if(VisIt_VariableMetaData_alloc(&vmd) == VISIT_OKAY)
            {
            VisIt_VariableMetaData_setName(vmd, name);
            VisIt_VariableMetaData_setMeshName(vmd, "mesh");
            VisIt_VariableMetaData_setType(vmd, VISIT_VARTYPE_SCALAR);
            VisIt_VariableMetaData_setCentering(vmd, VISIT_VARCENTERING_ZONE);

            VisIt_SimulationMetaData_addVariable(md, vmd);
            }
          }

        /* Add some custom commands. */
        for(i = 0; i < sizeof(cmd_names)/sizeof(const char *); ++i)
        {
            visit_handle cmd = VISIT_INVALID_HANDLE;
            if(VisIt_CommandMetaData_alloc(&cmd) == VISIT_OKAY)
            {
                VisIt_CommandMetaData_setName(cmd, cmd_names[i]);
                VisIt_SimulationMetaData_addGenericCommand(md, cmd);
            }
        }
    }

    return md;
}

visit_handle SimGetMesh(int domain, const char *name, void *cbdata)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    simulation_data *sim = (simulation_data *)cbdata;

    if(strcmp(name, "mesh") == 0)
    {
        if(VisIt_RectilinearMesh_alloc(&h) != VISIT_ERROR)
        {
            int minRealIndex[3], maxRealIndex[3];
            minRealIndex[0] = 0;
            minRealIndex[1] = 0;
            minRealIndex[2] = 0;
            maxRealIndex[0] = sim->grid.Ncolumns;
            maxRealIndex[1] = sim->grid.Nrows;
            maxRealIndex[2] = sim->grid.Nlevels;
#ifdef SHOW_GHOST_ARRAY
            // only modify in X and Y directions, since Z does not have ghost
            minRealIndex[0]++;
            minRealIndex[1]++;
            maxRealIndex[0]--;
            maxRealIndex[1]--;
#endif
    //fprintf(stderr,"minRealIndex[%dx%dx%d], maxRealIndex[%dx%dx%d]\n", minRealIndex[0], minRealIndex[1],minRealIndex[2],maxRealIndex[0], maxRealIndex[1],maxRealIndex[2]);
            visit_handle hxc, hyc, hzc;
            VisIt_VariableData_alloc(&hxc);
            VisIt_VariableData_alloc(&hyc);
            VisIt_VariableData_alloc(&hzc);
            VisIt_VariableData_setDataF(hxc, VISIT_OWNER_SIM, 1, sim->grid.Ncolumns+1, sim->grid.rmesh_x);
            VisIt_VariableData_setDataF(hyc, VISIT_OWNER_SIM, 1, sim->grid.Nrows+1,    sim->grid.rmesh_y);
            VisIt_VariableData_setDataF(hzc, VISIT_OWNER_SIM, 1, sim->grid.Nlevels+1,  sim->grid.rmesh_z);
            VisIt_RectilinearMesh_setCoordsXYZ(h, hxc, hyc, hzc);
#ifdef SHOW_GHOST_ARRAY
            VisIt_RectilinearMesh_setRealIndices(h, minRealIndex, maxRealIndex);
#endif
        }
    }
    return h;
}

visit_handle SimGetVariable(int domain, const char *name, void *cbdata)
{
  visit_handle h = VISIT_INVALID_HANDLE;
  int f, ret = -1, nComponents = 1, nTuples;
  simulation_data *sim = (simulation_data *)cbdata;

  if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
    {
    nTuples = (sim->grid.Ncolumns) * (sim->grid.Nrows) *  (sim->grid.Nlevels);
    for (f = 0; f < NFIELDS; f++)
      {
      char simdname[8];
      sprintf(simdname, "data%03d", f);
      if(strcmp(name, simdname) == 0)
        {
        ret = VisIt_VariableData_setDataF(h, VISIT_OWNER_SIM, nComponents, nTuples, sim->grid.data[f]);
        break;
        }
      }
     if(ret == VISIT_ERROR)
        {
        VisIt_VariableData_free(h);
        h = VISIT_INVALID_HANDLE;
        }
    }
    return h;
}
visit_handle SimGetDomainList(const char *name, void *cbdata)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    if(VisIt_DomainList_alloc(&h) != VISIT_ERROR)
    {
        visit_handle hdl;
        int i, *iptr = NULL;
        simulation_data *sim = (simulation_data *)cbdata;

        iptr = (int *)malloc(sizeof(int));
        *iptr = sim->par_rank;

        if(VisIt_VariableData_alloc(&hdl) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataI(hdl, VISIT_OWNER_VISIT, 1, 1, iptr);
            VisIt_DomainList_setDomains(h, sim->par_size, hdl);
        }
    }
    return h;
}

/* Callback function for control commands, which are the buttons in the
 * GUI's Simulation window. This type of command is handled automatically
 * provided that you have registered a command callback such as this.
 */
void ControlCommandCallback(const char *cmd, const char *args, void *cbdata)
{
    simulation_data *sim = (simulation_data *)cbdata;

    if(strcmp(cmd, "halt") == 0)
      {
      sim->runMode = SIM_STOPPED;
      }
    else if(strcmp(cmd, "step") == 0)
      {
      simulate_one_timestep(sim);
      }
    else if(strcmp(cmd, "run") == 0)
      {
      sim->runMode = SIM_RUNNING;
      }
    else if(strcmp(cmd, "Update ON/OFF") == 0)
      {
      sim->updateplots = sim->updateplots ? 0 : 1;
      }
    else if(strcmp(cmd, "Save ON/OFF") == 0)
      {
      if (sim->savingFiles) 
        sim->savingFiles = 0;
      if (!sim->savingFiles) 
        sim->savingFiles = 1;
      }

}

static int visit_broadcast_int_callback(int *value, int sender, void *cbdata)
{
  simulation_data *sim = (simulation_data *)cbdata;
  return MPI_Bcast(value, 1, MPI_INT, sender, sim->comm_cart);
}

static int visit_broadcast_string_callback(char *str, int len, int sender, void *cbdata)
{
  simulation_data *sim = (simulation_data *)cbdata;
  return MPI_Bcast(str, len, MPI_CHAR, sender, sim->comm_cart);
}

#define VISIT_COMMAND_PROCESS 0
#define VISIT_COMMAND_SUCCESS 1
#define VISIT_COMMAND_FAILURE 2

/* Helper function for ProcessVisItCommand */
static void BroadcastSlaveCommand(int *command, simulation_data *sim)
{
  MPI_Bcast(command, 1, MPI_INT, 0, sim->comm_cart);
}

/* Callback involved in command communication. */
void SlaveProcessCallback(void *cbdata)
{
   simulation_data *sim = (simulation_data *)cbdata;
   int command = VISIT_COMMAND_PROCESS;
   BroadcastSlaveCommand(&command, sim);
}

/* Process commands from viewer on all processors. */
int ProcessVisItCommand(simulation_data *sim)
{
    int command;
    if (sim->par_rank==0)
    {
        int success = VisItProcessEngineCommand();

        if (success == VISIT_OKAY)
        {
            command = VISIT_COMMAND_SUCCESS;
            BroadcastSlaveCommand(&command, sim);
            return 1;
        }
        else
        {
            command = VISIT_COMMAND_FAILURE;
            BroadcastSlaveCommand(&command, sim);
            return 0;
        }
    }
    else
    {
        /* Note: only through the SlaveProcessCallback callback
         * above can the rank 0 process send a VISIT_COMMAND_PROCESS
         * instruction to the non-rank 0 processes. */
        while (1)
        {
            BroadcastSlaveCommand(&command, sim);
            switch (command)
            {
            case VISIT_COMMAND_PROCESS:
                VisItProcessEngineCommand();
                break;
            case VISIT_COMMAND_SUCCESS:
                return 1;
            case VISIT_COMMAND_FAILURE:
                return 0;
            }
        }
    }
}

#endif
