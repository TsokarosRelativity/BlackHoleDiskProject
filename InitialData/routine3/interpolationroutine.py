import pandas as pd
import numpy as np
from scipy.spatial import KDTree
import argparse
import time
import concurrent.futures
from tqdm import tqdm


parser = argparse.ArgumentParser(
    prog="Interpolation Routine",
    description="Interpolating the 3D data",
    epilog="Routine takes in 3D initial data and interpolates the initial data to match grids for Evolution Routine",
)
parser.add_argument(
    "data_folder",
    type=str,
    help="Folder in which the Initial Data is being stored and calculated",
)
args = parser.parse_args()

print("======================AHHHHHHHHHHHH=======================")
s1 = time.time()
output_dir = "grid_data/"
gridfile = "grids_bh_disk_patrik"
print("LOADING IN DATA... processing")
df = pd.read_hdf(args.data_folder + "3D_data/all_data.h5", key="df")
df = df.map(lambda x: np.float64(0) if abs(x) < 1e-10 else x)

treemap = KDTree(df[["x", "y", "z"]].values)
print("FINISHED PREPROCESSING DATA :)")

cols = list(df.columns)

def process_line(line):
    parts = line.split()[1:]  # Ignore the first column
    new_grid = list(map(np.float64, parts))
    output_file = f"CTS_bin-proc{parts[-1]}.d"     
    return new_grid, output_file

def gridmaker(new_grid, outputfile, df, kval):
    ts = time.time()
    try:
        x_min, y_min, z_min, dx, dy, dz, Nx, Ny, Nz, _ = new_grid
        x_arr = np.arange(x_min, x_min + dx * Nx, dx)
        y_arr = np.arange(y_min, y_min + dy * Ny, dy)
        z_arr = np.arange(z_min, z_min + dz * Nz, dz)
        grid = np.array(np.meshgrid(x_arr, y_arr, z_arr)).T.reshape(-1, 3)
        dMat, iMat = treemap.query(grid, k=kval, workers=-1)
        interpolated_data = np.array([
                interpolate_grid(p, d, i, df) 
                for p, d, i in zip(grid, dMat, iMat)
                ])

        ntot = Nx * Ny * Nz
        with open(args.data_folder + output_dir + outputfile, "wb") as f:
            f.write(np.array([Nx, Ny, Nz, ntot], dtype=np.float64).tobytes())
            f.write(interpolated_data.astype(np.float64).tobytes())
        print(f"Finished {outputfile} in {time.time() - ts} sec")
    except Exception as e:
        print(f"Error in {outputfile}: {e}")
        print(f"time elapse: {time.time() - ts} sec")
    return 



def interpolate_grid(p, dist, ind, df):
    if np.any(dist == 0):
        t = df.iloc[ind[np.argmin(dist)]]
        return np.hstack([p, t.values[3:]])
    w = 1/dist
    w = w.reshape((-1, 1))
    r = np.sum(w * df.iloc[ind].values[:, 3:], axis=0) / np.sum(w)
    return np.hstack([p, r])

linedata = []
with open(gridfile, "r") as f:
    lines = f.readlines()
    linedata = [process_line(l) for l in lines]

with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = []
    for ld in tqdm(linedata, desc="Processing grids"):
        future = executor.submit(gridmaker, ld[0], ld[1], df, 8)
        futures.append(future)
    
    concurrent.futures.wait(futures)
#    with concurrent.futures.ThreadPoolExecutor() as executor:
#        # Submit all tasks to the executor
#        futures = [executor.submit(gridmaker, ld[0], ld[1], df, 8) for ld in linedata ]
#        
#        # Wait for all tasks to complete
#        concurrent.futures.wait(futures)



print(f"Finished all grids in {time.time() - s1} sec")



