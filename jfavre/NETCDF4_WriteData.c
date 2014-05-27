/*
// Test for parallel I/O benchmark
// Written by Jean M. Favre, Swiss National Supercomputer Center
// Last update: Wed Jun 8 2011
*/
void NETCDF4_WriteData(simulation_data *sim, char *Filename)
{
  int ret = NC_NOERR;
  int f, old_fill_mode;
  int ncid, dimuids[4], varid2[NFIELDS];
  size_t ustart[4], ucount[4];
  size_t chunksize[4] = {1, sim->global_dims[2], sim->global_dims[1], sim->global_dims[0]};
  MPI_Info info = MPI_INFO_NULL;
  char var_name2[NFIELDS][8];

  //fprintf(stderr,"opening %s, %s\n", Filename, "NC_NETCDF4|NC_COLLECTIVE");
#ifdef PARIO
  if(ret = nc_set_chunk_cache( sizeof(float) * sim->global_dims[2]*sim->global_dims[1]*sim->global_dims[0] * 1, 23, 0.0))
  //if(ret = nc_set_chunk_cache( sizeof(float) * LEVELS_SIZE*LATITUDE_SIZE*LONGITUDE_SIZE * 1, 1, 0.0))
    fprintf(stderr, "Error: %s\n", nc_strerror(ret));

  if(ret = nc_create_par(Filename, NC_NETCDF4|NC_MPIIO, sim->comm_cart, info, &ncid))
    fprintf(stderr, "Error: %s\n", nc_strerror(ret));

  nc_def_dim(ncid, "time", NC_UNLIMITED, &dimuids[0]);
  nc_def_dim(ncid, "levels", sim->global_dims[2], &dimuids[1]);
  nc_def_dim(ncid, "lat", sim->global_dims[1], &dimuids[2]);
  nc_def_dim(ncid, "long", sim->global_dims[0], &dimuids[3]);

  /* Create NFIELDS variables. */
  for (f = 0; f < NFIELDS; f++)
    {
    sprintf(var_name2[f], "data-%02d", f);
    if(ret = nc_def_var(ncid, var_name2[f], NC_FLOAT, 4, &dimuids[0], &varid2[f]))
      fprintf(stderr, "Error: %s\n", nc_strerror(ret));
    nc_def_var_fill(ncid, varid2[f], NC_NOFILL,  NULL);
    nc_var_par_access(ncid, varid2[f], NC_COLLECTIVE); //NC_INDEPENDENT);
    nc_def_var_chunking(ncid, varid2[f], NC_CHUNKED, chunksize);
    }
  nc_set_fill(ncid, NC_NOFILL, &old_fill_mode);
  nc_enddef(ncid);
#else
  if(sim->par_rank == 0)
    {
    if(ret = nc_set_chunk_cache( sizeof(float) * sim->global_dims[2]*sim->global_dims[1]*sim->global_dims[0]*1, 23, 0.0))
      fprintf(stderr, "Error: %s\n", nc_strerror(ret));
    if (ret = nc_create(Filename, NC_NETCDF4|NC_CLOBBER, &ncid))
      fprintf(stderr, "Error: %s\n", nc_strerror(ret));

    nc_def_dim(ncid, "time", NC_UNLIMITED, &dimuids[0]);
    nc_def_dim(ncid, "levels", sim->global_dims[2], &dimuids[1]);
    nc_def_dim(ncid, "lat", sim->global_dims[1], &dimuids[2]);
    nc_def_dim(ncid, "long", sim->global_dims[0], &dimuids[3]);

  /* Create NFIELDS variables. */
    for (f = 0; f < NFIELDS; f++)
      {
      sprintf(var_name2[f], "data-%02d", f);
      if(ret = nc_def_var(ncid, var_name2[f], NC_FLOAT, 4, &dimuids[0], &varid2[f]))
        fprintf(stderr, "Error: %s\n", nc_strerror(ret));
      nc_def_var_fill(ncid, varid2[f], NC_NOFILL,  NULL);
      nc_def_var_chunking(ncid, varid2[f], NC_CHUNKED, chunksize);
      }
    nc_set_fill(ncid, NC_NOFILL, &old_fill_mode);
    nc_enddef(ncid);
  } // only rank 0 creates and initializes file
#endif

/********************* end of defs ******************/

#ifdef PARIO
  ustart[0] = 0;
  ustart[1] = 0;
  ustart[2] = sim->grid.Nrows * sim->coords[1];
  ustart[3] = sim->grid.Ncolumns * sim->coords[0];

  ucount[0] = 1;
  ucount[1] = sim->grid.Nlevels;
  ucount[2] = sim->grid.Nrows;
  ucount[3] = sim->grid.Ncolumns;

  for(f=0; f < NFIELDS; f++)
    {
    if (ret = nc_put_vara_float(ncid, varid2[f], &ustart[0], &ucount[0], sim->grid.data[f]))
      {
      fprintf(stderr, "Error: %s\n", nc_strerror(ret));
      }
    }
 
#else
  int proc;

  MPI_Status status;

  if(sim->par_rank == 0)
    {
    ustart[0] = 0;
    ustart[1] = 0;
    ustart[2] = sim->grid.Nrows * sim->coords[1];
    ustart[3] = sim->grid.Ncolumns * sim->coords[0];

    ucount[0] = 1;
    ucount[1] = sim->grid.Nlevels;
    ucount[2] = sim->grid.Nrows;
    ucount[3] = sim->grid.Ncolumns;
  // first rank 0 writes its own data
    for(f=0; f < NFIELDS; f++)
      {
      if (ret = nc_put_vara_float(ncid, varid2[f], &ustart[0], &ucount[0], sim->grid.data[f]))
        {
        fprintf(stderr, "Error: %s\n", nc_strerror(ret));
        }
      }
  // then rank 0 allocates a buffer to receive from his neighbors, then writes their data
  int size = (sim->grid.Ncolumns) * (sim->grid.Nrows) * (sim->grid.Nlevels);
  float *buffer = (float *)malloc(sizeof(float)*size);
  for(proc=1; proc < sim->par_size; proc++)
    {
    int neighbors_coords[2];
    ustart[0] = 0;
    ustart[1] = 0;
    // next 3 lines are valid if using an MPI topology. Otherwise (for VisIt, use my own. TBD)
    MPI_Cart_coords(sim->comm_cart, proc, 2, neighbors_coords);
    ustart[2] = sim->grid.Nrows * neighbors_coords[1];
    ustart[3] = sim->grid.Ncolumns * neighbors_coords[0];
    //ustart[2] = sim->grid.Nrows * sim->coords[1];
    //ustart[3] = sim->grid.Ncolumns * sim->coords[0];

    for(f=0; f < NFIELDS; f++)
      {
      MPI_Recv(buffer, size, MPI_FLOAT, proc, f, sim->comm_cart, &status);
      //fprintf(stderr,"received from proc %d\n", proc);
      if ((ret = nc_put_vara_float(ncid, varid2[f], &ustart[0], &ucount[0], buffer)))
        {
        fprintf(stderr, "Error: %s\n", nc_strerror(ret));
        }
      }
    }
  free(buffer);
   }
 else // the case of all other MPI ranks different than 0
   {
    int size = (sim->grid.Ncolumns) * (sim->grid.Nrows) * (sim->grid.Nlevels);
    for(f=0; f < NFIELDS; f++)
      {
      MPI_Send(sim->grid.data[f], size, MPI_FLOAT, 0, f, sim->comm_cart);
      }
   }
#endif
  nc_close(ncid);
}

