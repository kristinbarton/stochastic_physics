import os
import numpy as np
import netCDF4 as nc

def generate_ca_test_files(res=96, ncells=4, nca=1, nca_g=1):
    os.makedirs("INPUT", exist_ok=True)
    os.makedirs("RESTART", exist_ok=True)
    
    # Grid specifications
    nx_global = res
    ny_global = res
    # Finer resolution calculated matching the Fortran domain setup
    nxc = (res // ncells) * ncells 
    nyc = (res // ncells) * ncells

    for tile in range(1, 7):
#        # 1. Create the base input condition file required for the run
#        cond_filename = f"INPUT/C{res}_ca_condition.tile{tile}.nc"
#        with nc.Dataset(cond_filename, "w", format="NETCDF4") as rootgrp:
#            # Dimensions
#            rootgrp.createDimension("grid_xt", nx_global)
#            rootgrp.createDimension("grid_yt", ny_global)
#            rootgrp.createDimension("time", None)
#            
#            # Variables
#            xt = rootgrp.createVariable("grid_xt", "f4", ("grid_xt",))
#            yt = rootgrp.createVariable("grid_yt", "f4", ("grid_yt",))
#            cond = rootgrp.createVariable("ca_condition", "f4", ("grid_xt", "grid_yt", "time"))
#            
#            # Populate baseline arrays
#            xt[:] = np.arange(1, nx_global + 1)
#            yt[:] = np.arange(1, ny_global + 1)
#            # Fill with 1.0 so CA calculations have valid initialization state weights
#            cond[:, :, 0] = np.ones((nx_global, ny_global), dtype=np.float32)
#        
#        print(f"Generated condition file: {cond_filename}")

        # 2. Create a dummy baseline mid_run restart file
        restart_filename = f"INPUT/ca_data.tile{tile}.nc"
        with nc.Dataset(restart_filename, "w", format="NETCDF4") as rootgrp:
            rootgrp.createDimension("nx_cell", nxc)
            rootgrp.createDimension("ny_cell", nyc)
            rootgrp.createDimension("nca", nca)
            rootgrp.createDimension("nca_g", nca_g)
            
            # Replicates arrays written out via write_ca_restart()
            iini = rootgrp.createVariable("iini", "i4", ("nx_cell", "ny_cell", "nca"))
            ilives = rootgrp.createVariable("ilives_in", "i4", ("nx_cell", "ny_cell", "nca"))
            
            iini[:, :, :] = np.zeros((nxc, nyc, nca), dtype=np.int32)
            ilives[:, :, :] = np.zeros((nxc, nyc, nca), dtype=np.int32)
            
        print(f"Generated dummy restart placeholder: {restart_filename}")

if __name__ == "__main__":
    # Matches the baseline RES=96 value from run_unit_tests_ca.sh
    generate_ca_test_files(res=96)
