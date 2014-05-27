adios_groupsize = 4 \
                + 4 \
                + 4 \
                + 4 \
                + 4 \
                + 4 \
                + 4 \
                + 4 \
                + 4 \
                + 4 * (sim->grid.Nlevels) * (sim->grid.Nrows) * (sim->grid.Ncolumns) \
                + 4 * (sim->grid.Nlevels) * (sim->grid.Nrows) * (sim->grid.Ncolumns) \
                + 4 * (sim->grid.Nlevels) * (sim->grid.Nrows) * (sim->grid.Ncolumns) \
                + 4 * (sim->grid.Nlevels) * (sim->grid.Nrows) * (sim->grid.Ncolumns) \
                + 4 * (sim->grid.Nlevels) * (sim->grid.Nrows) * (sim->grid.Ncolumns) \
                + 4 * (sim->grid.Nlevels) * (sim->grid.Nrows) * (sim->grid.Ncolumns) \
                + 4 * (sim->grid.Nlevels) * (sim->grid.Nrows) * (sim->grid.Ncolumns) \
                + 4 * (sim->grid.Nlevels) * (sim->grid.Nrows) * (sim->grid.Ncolumns) \
                + 4 * (sim->grid.Nlevels) * (sim->grid.Nrows) * (sim->grid.Ncolumns) \
                + 4 * (sim->grid.Nlevels) * (sim->grid.Nrows) * (sim->grid.Ncolumns);
adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
adios_write (adios_handle, "nx_global", &sim->global_dims[0]);
adios_write (adios_handle, "ny_global", &sim->global_dims[1]);
adios_write (adios_handle, "nz_global", &sim->global_dims[2]);
adios_write (adios_handle, "offs_x", &offs_x);
adios_write (adios_handle, "offs_y", &offs_y);
adios_write (adios_handle, "offs_z", &offs_z);
adios_write (adios_handle, "Ncolumns", &sim->grid.Ncolumns);
adios_write (adios_handle, "Nrows", &sim->grid.Nrows);
adios_write (adios_handle, "Nlevels", &sim->grid.Nlevels);
adios_write (adios_handle, "/data/data-00", sim->grid.data[0]);
adios_write (adios_handle, "/data/data-01", sim->grid.data[1]);
adios_write (adios_handle, "/data/data-02", sim->grid.data[2]);
adios_write (adios_handle, "/data/data-03", sim->grid.data[3]);
adios_write (adios_handle, "/data/data-04", sim->grid.data[4]);
adios_write (adios_handle, "/data/data-05", sim->grid.data[5]);
adios_write (adios_handle, "/data/data-06", sim->grid.data[6]);
adios_write (adios_handle, "/data/data-07", sim->grid.data[7]);
adios_write (adios_handle, "/data/data-08", sim->grid.data[8]);
adios_write (adios_handle, "/data/data-09", sim->grid.data[9]);
