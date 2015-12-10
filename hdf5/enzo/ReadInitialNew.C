// sed "s-CSCS CSCS-64-" ReadInitialNew.templatejg > ReadInitialNew.C

#ifdef USEMPI
#include "mpi.h"
#endif

#define NPES1D 8

#include "hdf5.h"
#include <stdio.h>
#include <iostream>
#include <string>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
using namespace std;

main(int argc, char* argv[]) {

  int mype;
  int npes;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  if(mype == 0) printf("PARALLEL MPI EXECUTION WITH %d PROCESSORS\n", npes);
  string inputfilename;
  string objname;
  long nx=2400;
  long ndata=nx*nx*nx;

  inputfilename = argv[1];
  objname = argv[2];
  if(mype == 0)printf("%s\n",inputfilename.c_str());
  if(mype == 0)printf("%s\n",objname.c_str());

  hid_t file_id, dset_id, mem_type_id; 
  herr_t status;
  herr_t h5_status;

  hid_t plist_id, pobj_id;
  MPI_Info info = MPI_INFO_NULL;
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  //H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
  H5Pset_fapl_mpiposix(plist_id, MPI_COMM_WORLD, 0);
  file_id = H5Fopen(inputfilename.c_str(), H5F_ACC_RDONLY, plist_id);
  if(mype == 0)printf("File open\n"); 

  pobj_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(pobj_id, H5FD_MPIO_COLLECTIVE);
  dset_id = H5Dopen2(file_id, objname.c_str(), H5P_DEFAULT);
  if(mype == 0)printf("Object open\n"); 

// FROM HERE MIMICING ORIGINAL READFILE.C IN ENZO
// interface parameters
    int Rank;
    

//
    hid_t      file_dsp_id, mem_dsp_id;
    hsize_t    slab_stride[4], slab_count[4], slab_block[4];
    hsize_t    slab_offset[4];
    hsize_t    xfer_size;

/*
    long nxpes = nx/npes;
    long nxpes0=nxpes;
*/
    int Part=0;
    slab_offset[0] = 0;
    slab_stride[0] = 1;
    slab_count[0] = 1;
    slab_block[0] = 1;

/* The following setup is for planar domain decomposition
    long nxpes = nx/(npes-1);
    long nxpes0=nxpes;
    if(mype == npes-1)nxpes=nx%(npes-1);

    slab_offset[1] = mype*nxpes0;
    slab_stride[1] = 1;
    slab_count[1] = nxpes;
    slab_block[1] = 1;
    slab_offset[2] = 0;
    slab_stride[2] = 1;
    slab_count[2] = nx;
    slab_block[2] = 1;
    slab_offset[3] = 0;
    slab_stride[3] = 1;
    slab_count[3] = nx;
    slab_block[3] = 1;
*/

/* The following setup is for cubic domain decomposition */
 
    int npes1D=NPES1D;
    int ix = mype/(npes1D*npes1D);
    int iy = (mype%(npes1D*npes1D))/npes1D;
    int iz = mype-iy*npes1D-ix*npes1D*npes1D;
    MPI_Barrier(MPI_COMM_WORLD);
/*
    long nxpesx = nx/(npes1D-1);
    long nxpesy = nx/(npes1D-1);
    long nxpesz = nx/(npes1D-1);
*/
    long nxpesx = nx/npes1D;
    long nxpesy = nx/npes1D;
    long nxpesz = nx/npes1D;
    long nxpesx0=nxpesx;
    long nxpesy0=nxpesy;
    long nxpesz0=nxpesz;
/*
    if(ix == npes1D-1)nxpesx=nx%(npes-1);
    if(iy == npes1D-1)nxpesy=nx%(npes-1);
    if(iz == npes1D-1)nxpesz=nx%(npes-1);
*/
    //printf(">>> %d %d %d %d %ld %ld %ld\n", mype, ix, iy, iz, nxpesx, nxpesy,nxpesz);

    slab_offset[1] = ix*nxpesx0;
    slab_stride[1] = 1;
    slab_count[1] = nxpesx;
    slab_block[1] = 1;
    slab_offset[2] = iy*nxpesy0;
    slab_stride[2] = 1;
    slab_count[2] = nxpesy;
    slab_block[2] = 1;
    slab_offset[3] = iz*nxpesz0;
    slab_stride[3] = 1;
    slab_count[3] = nxpesz;
    slab_block[3] = 1;


    if(mype == 0)printf("End of setup\n"); 

    // Data in memory is considered 1D, stride 1, with zero offset

    xfer_size = nxpesx*nxpesy*nxpesz;
    if(mype == 0)printf("Buffer allocated\n"); 
    hsize_t    mem_stride[1], mem_count[1], mem_block[1];
    hsize_t    mem_offset[1];
    mem_stride[0] = 1;           // contiguous elements
    mem_count[0] = xfer_size;    // number of elements in field
    mem_offset[0] = 0;           // zero offset in buffer
    mem_block[0] = 1;            // single element blocks

    hsize_t     Slab_Dims[4];
    int         Slab_Rank=4;
    Slab_Dims[0]=1;
    Slab_Dims[1]=nx;
    Slab_Dims[2]=nx;
    Slab_Dims[3]=nx;
    double * tempbuffer = new double [xfer_size];
    //printf("DATA SIZE OF %d = %ld\n", mype,xfer_size);

    double t0 = MPI_Wtime();
    mem_dsp_id = H5Screate_simple(1, mem_count, NULL);
    h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, mem_offset, mem_stride, mem_count, mem_block);
    file_dsp_id = H5Screate_simple(Slab_Rank, Slab_Dims, NULL);
    h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, slab_block);
    if(mype == 0)
    printf("Start reading (PE %d)\n", mype); 
    h5_status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, mem_dsp_id, file_dsp_id, pobj_id, tempbuffer);
    if(mype == 0)
    printf("End reading (PE %d)\n", mype); 

    //if(mype == 0)
    printf("%d %f %f\n", mype, tempbuffer[0], tempbuffer[100]); 
    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    double ttot=t1-t0;
    if(mype == 0)
    printf("TIME TO LOAD DATA %f\n",ttot); 
    
    H5Sclose(mem_dsp_id);
    H5Sclose(file_dsp_id);
    H5Pclose(plist_id);
    H5Pclose(pobj_id);
    H5Dclose(dset_id);  
    H5Fclose(file_id);

    MPI_Finalize();

}
