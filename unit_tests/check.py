import numpy as np
import netCDF4 as nc

def validate_gaussian_grid(lons, lats):
    print("--- Gaussian Grid Validation Report ---\n")
    
    n_lat = len(lats)
    n_lon = len(lons)
    n = n_lat // 2
    print(f"Detected Shape: N={n} grid ({n_lat} latitudes, {n_lon} longitudes)")
    
    # 1. Pole Check
    has_poles = np.any(np.isclose(np.abs(lats), 90.0))
    print(f"[1] No Exact Pole Points:  {'FAIL' if has_poles else 'PASS'}")
    
    # 2. Symmetry Check
    is_symmetric = np.allclose(lats, -lats[::-1], atol=1e-10)
    print(f"[2] Equator Symmetry:      {'PASS' if is_symmetric else 'FAIL'}")
    
    # 3. Longitude Spacing
    lon_diffs = np.diff(lons)
    expected_spacing = 360.0 / n_lon
    is_uniform_lon = np.allclose(lon_diffs, expected_spacing, atol=1e-10)
    print(f"[3] Uniform Longitudes:    {'PASS' if is_uniform_lon else 'FAIL'} (Expected spacing: {expected_spacing:.5f}°)")
    
    # 4. The Mathematical Check (Legendre Roots)
    mu_computed = np.sin(np.radians(lats))
    roots_exact, _ = np.polynomial.legendre.leggauss(n_lat)
    max_error = np.max(np.abs(mu_computed - roots_exact))
    
    tolerance = 1e-12
    is_gaussian = max_error < tolerance
    
    print(f"[4] Legendre Root Match:   {'PASS' if is_gaussian else 'FAIL'}")
    print(f"    -> Maximum deviation from ideal root: {max_error:.4e}")
    
    print("\n--- Final Conclusion ---")
    if is_gaussian and is_symmetric and is_uniform_lon and not has_poles:
        print("✅ VALID: This is a mathematically perfect Regular Gaussian grid.")
    else:
        print("❌ INVALID: This grid failed one or more Gaussian geometry tests.")


# --- Read Data from NetCDF ---
file_path = "run_ggca/ggca_lats.nc"

try:
    # Open the dataset in read mode
    dataset = nc.Dataset(file_path, 'r')
    
    # Extract the variables as numpy arrays
    # Using [:] ensures we pull the raw data values out of the NetCDF variable object
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    
    # Close the file safely
    dataset.close()
    
    print(f"Successfully loaded data from {file_path}\n")
    
    # Run the validation
    validate_gaussian_grid(lon, lat)
    
except FileNotFoundError:
    print(f"Error: Could not find the file '{file_path}'. Please check the path.")
except KeyError as e:
    print(f"Error: Missing expected variable in NetCDF file - {e}")
