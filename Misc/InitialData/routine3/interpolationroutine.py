import pandas as pd
import numpy as np
from scipy.spatial import KDTree
import argparse
import time
from joblib import Parallel, delayed
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
df = df.mask(abs(df) < 1e-14, np.float64(0))
df = df.drop_duplicates(subset=["x", "y", "z"], keep="last")
treemap = KDTree(df[["x", "y", "z"]].to_numpy())
print("FINISHED PREPROCESSING DATA :)")

cols = list(df.columns)

def process_line(line):
    parts = line.split()[1:]  # Ignore the first column
    new_grid = list(map(np.float64, parts))
    output_file = f"CTS_bin-proc{parts[-1]}.d"     
    return new_grid, output_file

def gridmaker(new_grid, outputfile, df):
    ts = time.time()
    try:
        x_min, y_min, z_min, dx, dy, dz, Nx, Ny, Nz, _ = new_grid
        x_arr = np.arange(x_min, x_min + dx * Nx, dx)
        y_arr = np.arange(y_min, y_min + dy * Ny, dy)
        z_arr = np.arange(z_min, z_min + dz * Nz, dz)
        grid = np.array(np.meshgrid(x_arr, y_arr, z_arr)).T.reshape(-1, 3)
        inds = [locatingpoints(p, treemap) for p in grid]
        interpolated_data = np.array([
            interpolate_grid(p, i, df)
            for p, i in zip(grid, inds)
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

def locatingpoints(p, tree):
    dis, ind = tree.query(p, k=8, workers=-1)
    if dis[0] == 0:
        return np.array([ind[0]])
    mxd = np.max(dis)
    idxs = tree.query_ball_point(p, r=mxd, workers=-1)
    return idxs


def interpolate_grid(p, ind, df):
    if len(ind) == 1:
        return df.iloc[ind].to_numpy()
    data = df.iloc[ind].to_numpy()
    dis = np.linalg.norm(data[:, :3] - p, ord=2, axis=1)
    wet = 1/dis
    wet = wet.reshape((-1, 1))
    r = np.sum(wet*data[:, 3:], axis=0) / np.sum(wet) 
    return np.hstack([p, r])

linedata = []
with open(gridfile, "r") as f:
    lines = f.readlines()
    linedata = [process_line(l) for l in lines]
    print("LINE DATA LOADED AHHHHHHHHHHHHHHHH")


r = Parallel(n_jobs=-1)(delayed(gridmaker)(l[0], l[1], df) for l in linedata)
#
#with open(gridfile, "r") as f:
#    lines = f.readlines()
#    r = Parallel(n_jobs=-1)(delayed(grid_writer)(l, df) for l in lines)
#

#with concurrent.futures.ThreadPoolExecutor() as executor:
#    futures = []
#    for ld in tqdm(linedata, desc="Processing grids"):
#        future = executor.submit(gridmaker, ld[0], ld[1], df)
#        futures.append(future)
#    
#    concurrent.futures.wait(futures)



print(f"Finished all grids in {time.time() - s1} sec")

