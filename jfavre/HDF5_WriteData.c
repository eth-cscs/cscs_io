/*
// Test for parallel I/O benchmark
// Written by Jean M. Favre, Swiss National Supercomputer Center
// Last update: Wed Jun 8 2011
*/
void HDF5_WriteData(simulation_data *sim, char *Filename)
{
  int f;
  hid_t fapl_id, file_id, dset_id;
  hid_t filespace, memspace;
  hsize_t ucount[3], ustart[3];
  hsize_t dimuids[3]={sim->global_dims[2], sim->global_dims[1], sim->global_dims[0]};
  MPI_Info info = MPI_INFO_NULL;
  herr_t status;
  hid_t	plist_id;
#define  N 6 // for bounding box
// Set up file access property list with parallel I/O access

  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, sim->comm_cart, info);

// Create a new file collectively and release property list identifier.

  file_id = H5Fcreate(Filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

// add an attribute containing 6 floats, it is written collectively, i.e. all MPI tasks must do it
   double attr_data[N] = {0.1, 0.2, 0.3, 1.1, 2.2, 3.3}; // a bounding box with lower left and upper right corners
   hsize_t dims = N;
   hid_t dataspace_id = H5Screate_simple(1, &dims, NULL);
   hid_t attribute_id = H5Acreate(file_id, "bbox", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, attr_data);
   status = H5Aclose(attribute_id);
   status = H5Sclose(dataspace_id);

// Create the dataspace for the dataset.
  filespace = H5Screate_simple(3, dimuids, NULL); 

// Each process defines dataset in memory and writes it to the hyperslab
  ustart[0] = 0;
  ustart[1] = sim->grid.Nrows * sim->coords[1];
  ustart[2] = sim->grid.Ncolumns * sim->coords[0];

  ucount[0] = sim->grid.Nlevels;
  ucount[1] = sim->grid.Nrows;
  ucount[2] = sim->grid.Ncolumns;

  memspace = H5Screate_simple(3, ucount, NULL);

// Create the dataset with default properties and close filespace.
  char var_name[16];
  for(f=0; f < NFIELDS; f++)
    {
    sprintf(var_name, "data-%02d", f);
    dset_id = H5Dcreate(file_id, var_name, H5T_NATIVE_FLOAT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

// Select hyperslab in the file.
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, ustart, NULL, ucount, NULL);

// Create property list for collective dataset write.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    if(status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
		      plist_id, sim->grid.data[f]) < 0)
      {
      fprintf(stderr, "Error writing dataset: %s\n", var_name);
      }

    H5Dclose(dset_id);
    }

  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);
  H5Fclose(file_id);
}
