# Create animation of CA on GG and FV3 grids
import argparse
import os
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from PIL import Image

def main(args):
    # Open gaussian grid dataset
    ggfile = f"{args.dir}/ca_restart.nc"

    print(f"Gaussian grid: {ggfile}")    
    ggds = xr.open_dataset(ggfile)
    ggfield = ggds['field'].isel(nca=args.nca)

    # Get min/max values for colorbar
    vmin = ggfield.min().values
    vmax = ggfield.max().values

    # Gather  FV3 tiles 
    print(f"FV3 grid: {args.dir}/ca_out.tile*.nc")
    fv3_dss = []
    fv3_lls = []
    for i in range(1,7):
        fv3file = f"{args.dir}/ca_out.tile{i}.nc"
        fv3grid = f"{args.gridpre}{i}.nc"
        fv3_dss.append(xr.open_dataset(fv3file))
        fv3_lls.append(xr.open_dataset(fv3grid))

    # Prep three different projections
    proj_np  = ccrs.Orthographic(central_longitude=-90.0, central_latitude=90.0) # North pole centered
    proj_m90 = ccrs.Orthographic(central_longitude=-90.0, central_latitude=0.0 ) # Equatorial, US-centered
    proj_p90 = ccrs.Orthographic(central_longitude=90.0,  central_latitude=0.0 ) # Equatorial, Asia-centered
    projs = [proj_np, proj_m90, proj_p90]

    pngfiles = []

    for t in range (ggds.sizes['time']):
        print(f"Generating frame {t+1}/{ggds.sizes['time']}", end='\r')
        fig = plt.figure(figsize=(16,10))

        time_slice = ggfield.isel(time=t)

        # This cyclic point prevents blank seam line in global plots (?)
        data_cyc, lons_cyc = add_cyclic_point(time_slice.values, coord=ggds.lon.values)
        lats = ggds.lat.values

        fv3_axs = []
        for col in range(3):
            for row in range(2):
                ax = fig.add_subplot(2, 3, row*3 + col + 1, projection=projs[col])
                ax.set_global()
                ax.coastlines(linewidth=1.0,color='white')
                if (col==0):
                    row_label = "Gaussian Grid" if row==0 else "FV3 Grid"
                    ax.text(-0.15, 0.5, row_label, transform=ax.transAxes, rotation=90, va='center', ha='center', fontsize=14)
                if (row==0): # Gaussian Grid
                    im = ax.pcolormesh(lons_cyc, lats, data_cyc, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap='viridis', shading='auto')
                else: # FV3 Grid 
                    fv3_axs.append(ax)

        # Plot all FV3 tiles
        for tile in range(6):
            data_fv3 = fv3_dss[tile][f"ca{args.nca}"].isel(time=t).values
            lons_fv3 = fv3_lls[tile].geolon.values
            lats_fv3 = fv3_lls[tile].geolat.values

            for ax in fv3_axs:
                im = ax.pcolormesh(lons_fv3, lats_fv3, data_fv3, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap='viridis', shading='auto')


        fig.suptitle(f"Time Step: {t}", fontsize=16, y=0.95)

        # Save temp PNGs to stitch into GIF later
        fname = f"tempframe_{t:04d}.png"
        plt.savefig(fname, bbox_inches='tight', dpi=100)
        if (t==1):
            plt.show()
        plt.close(fig)
        pngfiles.append(fname)

    fps = 5
    frame_duration = int(1000 / fps) # milliseconds
    frames = [Image.open(f) for f in pngfiles]
    frames[0].save(args.output, format='GIF', append_images=frames[1:], save_all=True, duration=frame_duration, loop=0) 

    for frame in frames:
        frame.close()
    for filename in pngfiles:
        os.remove(filename)

    print(f"\nDone! GIF saved to: {args.output}")

if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser(description="Generate CA output animation")
    parser.add_argument("--dir", default="../run_ggca/RESTART", help="Path to output files")
    parser.add_argument("--nca", default=1,type=int, help="Which CA number to plot")
    parser.add_argument("--gridpre", default='../run_ggca/INPUT/C96.mx025_oro_data.tile', help="/path/tile prefix of location containing grid lat/lon data")
    parser.add_argument("-o", "--output", default="output.gif", help="Output GIF file location")
    args = parser.parse_args()

    main(args)
