import pandas as pd
import numpy as np
from scipy.spatial import KDTree
from scipy.interpolate import RBFInterpolator, BarycentricInterpolator
from numba import njit, prange
from concurrent.futures import ThreadPoolExecutor
import argparse
import sys

parser = argparse.ArgumentParser(
    prog="Interpolation Three",
    description="interpolation script thats lowkey cheating but whatevs",
    epilog="I couldnt quite interpolate stuff quickly so im doing a bunch of them all in batches",
)
parser.add_argument(
    "gridfile",
    type=str,
)

args = parser.parse_args()
print("Loading in data...")
df_3d = pd.read_hdf("/data/sjammi6/thesisproject/data/3D_data/all_data_routine1.h5", key="df")
df_3d = df_3d.drop_duplicates()
df_3d = df_3d.sort_values(by=["r", "theta", "phi"], ignore_index=True)

# Precompute unique arrays
radarr = np.sort(df_3d['r'].unique())
thetaarr = np.sort(df_3d['theta'].unique())
phiarr = np.sort(df_3d['phi'].unique())

# Convert DataFrame to NumPy array
df_3d_np = df_3d.to_numpy()
print("Data loaded successfully!")

print("Creating KD-tree...")
points = df_3d[["r", "theta", "phi"]].values
# Create KD-tree
tree = KDTree(points)

print("KD-tree created successfully!")


def cartesiantospherical(arr):
    r = np.sqrt(arr[:, 0]**2 + arr[:, 1]**2 + arr[:, 2]**2)
    theta = np.where(r != 0.0, np.arccos(np.clip(arr[:, 2] / r, -1.0, 1.0)), 0.0)
    phi = np.arctan2(arr[:, 1], arr[:, 0])
    new_grid = np.column_stack((r, theta, phi))
    return new_grid

@njit(parallel=True)
def cartesiantospherical_parallel(arr):
    n = arr.shape[0]
    new_grid = np.zeros((n, 3), dtype=np.float64)
    for i in prange(n):  # Parallel loop
        x, y, z = arr[i]
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arccos(z / r) if r != 0 else 0.0
        phi = np.arctan2(y, x)
        new_grid[i, 0] = 0 if np.abs(r) < 1e-14 else r
        new_grid[i, 1] = 0 if np.abs(theta) < 1e-14 else theta
        new_grid[i, 2] = 0 if np.abs(phi) < 1e-14 else phi
    return new_grid

def process_line(line):
    try:
        xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, MPI_ID = map(np.float64, line.split()[1:]) 
        xarr = np.linspace(xmin, xmin + nx*dx, int(nx), dtype=np.float64)
        yarr = np.linspace(ymin, ymin + ny*dy, int(ny), dtype=np.float64)
        zarr = np.linspace(zmin, zmin + nz*dz, int(nz), dtype=np.float64)    
        cartgrid = np.array(np.meshgrid(xarr, yarr, zarr)).T.reshape(-1, 3)
        sphericalgrid = cartesiantospherical(cartgrid)
        output_file = f"CTS_bin-proc{int(MPI_ID)}.d"
        return [cartgrid, sphericalgrid, nx, ny, nz, output_file]
    except:
        print("this line could not be processed:")
        print(line)
        sys.exit()
        return

def interpolate_point_parallel(point):
    interpdata = np.zeros(27)
    interpdata[0:3] = point
    rid = np.searchsorted(radarr, point[0])
    thetad = np.searchsorted(thetaarr, point[1])
    phid = np.searchsorted(phiarr, point[2])
    rmatch = np.isclose(radarr[rid], rid, atol=1e-14)
    thetamatch = np.isclose(thetaarr[thetad], thetad, atol=1e-14)
    phimatch = np.isclose(phiarr[phid], phid, atol=1e-14)

    #case 1
    if rmatch and thetamatch and phimatch:
        _, i = tree.query(point)
        return df_3d_np[i]
    #case 2
    if rmatch and thetamatch and not phimatch:
        ps = np.zeros(3)
        if phid < 2:
            ps = np.array([0, 1, 2])
        if phid > len(phiarr) - 3:
            ps = np.array([len(phiarr) - 1, len(phiarr) - 2, len(phiarr) - 3])
        else:
            ps = np.array([phid - 1, phid, phid + 1])

        coords = np.array(np.meshgrid(radarr[rid], thetaarr[thetad], phiarr[ps])).T.reshape(-1, 3)
        _, ind = tree.query(coords, workers=-1)
        interpdata[3:] = np.array([BarycentricInterpolator(df_3d_np[ind, 2], df_3d_np[ind, i], interpdata[2]) for i in range(3, 26)], dtype=np.float64)
        return interpdata 
    #case 3
    if rmatch and not thetamatch and phimatch:
        ts = np.zeros(3)
        if thetad < 2:
            ts = np.array([0, 1, 2])
        if thetad > len(thetaarr) - 3:
            ts = np.array([len(thetaarr) - 1, len(thetaarr) - 2, len(thetaarr) - 3])
        else:
            ts = np.array([phid - 1, phid, phid + 1])
        coords = np.array(np.meshgrid(radarr[rid], thetaarr[ts], phiarr[phid])).T.reshape(-1, 3)
        _, ind = tree.query(coords, workers=-1)
        interpdata[3:] = np.array([BarycentricInterpolator(df_3d_np[ind, 1], df_3d_np[ind, i], interpdata[1]) for i in range(3, 26)], dtype=np.float64)
        return interpdata 
    #case 4
    if not rmatch and thetamatch and phimatch:
        rs = np.zeros(3)
        if rid < 2:
            rs = np.array([0, 1, 2])
        if rid > len(radarr) - 3:
            rs = np.array([len(radarr) - 1, len(radarr) - 2, len(radarr) - 3])
        else:
            rs = np.array([rid - 1, rid, rid + 1])

        coords = np.array(np.meshgrid(radarr[rs], thetaarr[thetad], phiarr[phid])).T.reshape(-1, 3)
        _, ind = tree.query(coords, workers=-1)

        interpdata[3:] = np.array([BarycentricInterpolator(df_3d_np[ind, 0], df_3d_np[ind, i], interpdata[0]) for i in range(4, 27)], dtype=np.float64)
        return interpdata 
    #case 5
    if rmatch and not thetamatch and not phimatch:
        ps = np.zeros(3)
        if phid < 2:
            ps = np.array([0, 1, 2])
        if phid > len(phiarr) - 3:
            ps = np.array([len(phiarr) - 1, len(phiarr) - 2, len(phiarr) - 3])
        else:
            ps = np.array([phid - 1, phid, phid + 1])
        ts = np.zeros(3)
        if thetad < 2:
            ts = np.array([0, 1, 2])
        if thetad > len(thetaarr) - 3:
            ts = np.array([len(thetaarr) - 1, len(thetaarr) - 2, len(thetaarr) - 3])
        else:
            ts = np.array([phid - 1, phid, phid + 1])
        coords = np.array(np.meshgrid(radarr[rid], thetaarr[ts], phiarr[ps])).T.reshape(-1, 3)
        _, ind = tree.query(coords, workers=-1)

        rbfinterp = RBFInterpolator(df_3d_np[ind, :3], df_3d_np[ind, 3:], kernel="linear")
        
        interpdata[3:] = rbfinterp(np.array([point]))[0]
        return interpdata
    #case 6
    if not rmatch and thetamatch and not phimatch:
        rs = np.zeros(3)
        if rid < 2:
            rs = np.array([0, 1, 2])
        if rid > len(radarr) - 3:
            rs = np.array([len(radarr) - 1, len(radarr) - 2, len(radarr) - 3])
        else:
            rs = np.array([rid - 1, rid, rid + 1])
        ps = np.zeros(3)
        if phid < 2:
            ps = np.array([0, 1, 2])
        if phid > len(phiarr) - 3:
            ps = np.array([len(phiarr) - 1, len(phiarr) - 2, len(phiarr) - 3])
        else:
            ps = np.array([phid - 1, phid, phid + 1])
        coords = np.array(np.meshgrid(radarr[rs], thetaarr[thetad], phiarr[ps])).T.reshape(-1, 3)
        _, ind = tree.query(coords, workers=-1)
        rbfinterp = RBFInterpolator(df_3d_np[ind, :3], df_3d_np[ind, 3:], kernel="linear")
        
        interpdata[3:] = rbfinterp(np.array([point]))[0]
        return interpdata
    #case 7
    if not rmatch and not thetamatch and phimatch:
        ts = np.zeros(3)
        if thetad < 2:
            ts = np.array([0, 1, 2])
        if thetad > len(thetaarr) - 3:
            ts = np.array([len(thetaarr) - 1, len(thetaarr) - 2, len(thetaarr) - 3])
        else:
            ts = np.array([phid - 1, phid, phid + 1])
        rs = np.zeros(3)
        if rid < 2:
            rs = np.array([0, 1, 2])
        if rid > len(radarr) - 3:
            rs = np.array([len(radarr) - 1, len(radarr) - 2, len(radarr) - 3])
        else:
            rs = np.array([rid - 1, rid, rid + 1])
        coords = np.array(np.meshgrid(radarr[rs], thetaarr[thetad], phiarr[phid])).T.reshape(-1, 3)
        _, ind = tree.query(coords, workers=-1)
        rbfinterp = RBFInterpolator(df_3d_np[ind, :3], df_3d_np[ind, 3:], kernel="linear")
        
        interpdata[3:] = rbfinterp(np.array([point]))[0]
        return interpdata
    #case 8
    else:
        ps = np.zeros(3)
        if phid < 2:
            ps = np.array([0, 1, 2])
        if phid > len(phiarr) - 3:
            ps = np.array([len(phiarr) - 1, len(phiarr) - 2, len(phiarr) - 3])
        else:
            ps = np.array([phid - 1, phid, phid + 1])
        ts = np.zeros(3)
        if thetad < 2:
            ts = np.array([0, 1, 2])
        if thetad > len(thetaarr) - 3:
            ts = np.array([len(thetaarr) - 1, len(thetaarr) - 2, len(thetaarr) - 3])
        else:
            ts = np.array([phid - 1, phid, phid + 1])
        rs = np.zeros(3)
        if rid < 2:
            rs = np.array([0, 1, 2])
        if rid > len(radarr) - 3:
            rs = np.array([len(radarr) - 1, len(radarr) - 2, len(radarr) - 3])
        else:
            rs = np.array([rid - 1, rid, rid + 1])
        coords = np.array(np.meshgrid(radarr[rs], thetaarr[ts], phiarr[ps])).T.reshape(-1, 3)
        _, ind = tree.query(coords, workers=-1)
        rbfinterp = RBFInterpolator(df_3d_np[ind, :3], df_3d_np[ind, 3:], kernel="linear")
        
        interpdata[3:] = rbfinterp(np.array([point]))[0]
        return interpdata


    #return [cartgrid, sphericalgrid, nx, ny, nz, output_file]
def routine(gridparam):
    interpolated_data = np.array([
                interpolate_point_parallel(d[:3])
                for d in gridparam[1]
            ], dtype=np.float64)
    cartgrid = gridparam[0]
    interpolated_data[:, :3] = cartgrid
    outputfile = gridparam[5]
    nx = gridparam[2]
    ny = gridparam[3]
    nz = gridparam[4]
    ntot = nx*ny*nz

    with open(outputfile, "wb") as f:
        f.write(np.array([nx], dtype=np.int32).tobytes())
        f.write(np.array([ny], dtype=np.int32).tobytes())
        f.write(np.array([nz], dtype=np.int32).tobytes())
        f.write(np.array([ntot], dtype=np.int32).tobytes())
        f.write(interpolated_data.astype(np.float64).tobytes())
    return



# preprocess data
grid_data = []
lines = []
print("Processing Grid Parameters...")

with open(args.gridfile, "r") as f:
    lines = [l for l in f]

with ThreadPoolExecutor() as executor:
    grid_data = list(executor.map(process_line, lines))

print("Grid Parameters loaded!")

for g in grid_data:
    routine(g)

print("finished")
