import pandas as pd 
import numpy as np
from scipy.interpolate import RBFInterpolator
from scipy.spatial import KDTree
from scipy.interpolate import BarycentricInterpolator
from tqdm import tqdm

print("Loading in data...")
df_3d = pd.read_hdf("/data/sjammi6/thesisproject/data/3D_data/all_data_routine2.h5", key="df")
df_3d = df_3d.drop_duplicates()
df_3d = df_3d.sort_values(by=["r", "theta", "phi"], ignore_index=True)
rad_arr = np.sort(df_3d.r.unique())
theta_arr = np.sort(df_3d.theta.unique())
phi_arr = np.sort(df_3d.phi.unique())
print("Data loaded successfully!")
print("Creating KD-tree...")
# Create a KD-tree for efficient nearest-neighbor lookups
points = df_3d[["r", "theta", "phi"]].values
tree = KDTree(points)
print("KD-tree created successfully!")

def get_nearest_point_index(r, theta, phi, tree):
    """Find the index of the nearest point within a given tolerance."""
    _, index = tree.query([r, theta, phi])
    return index

def interpolationpoint(radarr, thetaarr, phiarr, point, tree, df):
    rp, thetap, phip = point
    ppoint = np.array([rp, thetap, phip]).reshape(1, -1)
    rid = np.searchsorted(radarr, rp)
    thetad = np.searchsorted(thetaarr, thetap)
    phid = np.searchsorted(phiarr, phip)
    rmatch = np.isclose(radarr[rid], rp, atol=1e-14)
    thetamatch = np.isclose(thetaarr[thetad], thetap, atol=1e-14)
    phimatch = np.isclose(phiarr[phid], phip, atol=1e-14)
    
    # case matching
    if rmatch and thetamatch and phimatch:
        id = get_nearest_point_index(ppoint[0], ppoint[1], ppoint[2], tree)
        return df.iloc[id].to_numpy()
    if rmatch and thetamatch and not phimatch:
        ps = []
        # edge cases:
        if phid < 2:
            ps = [0, 1, 2, 3, 4]
        elif phid > len(phiarr) - 3:
            ps = [len(phiarr) - 5, len(phiarr) - 4, len(phiarr) - 3, len(phiarr) - 2, len(phiarr) - 1]
        else:
            ps = [phid - 2, phid - 1, phid, phid + 1, phid + 2]
        
        cords = []
        for p in ps:
            cords.append([radarr[rid], thetaarr[thetad], phiarr[p]])
        id = [get_nearest_point_index(c[0], c[1], c[2], tree) for c in cords]
        data = df.iloc[id].sort_values(by=["r", "theta", "phi"]).drop_duplicates().to_numpy()
        
        
        vals = []
        for i in range(3, 27):
            interp = BarycentricInterpolator(data[:, 2], data[:, i])
            vals.append(interp(phip))
        vals = np.array(vals).reshape(1, -1)
        #interp = RBFInterpolator(data[:, :3], data[:, 3:], kernel="linear")
        #vals = interp(np.array([[rp, thetap, phip]]).reshape(1, -1))
        return vals

    if rmatch and not thetamatch and phimatch:
        ts = []
        # edge cases:
        if thetad < 2:
            ts = [0, 1, 2, 3, 4]
        elif thetad > len(thetaarr) - 3:
            ts = [len(thetaarr) - 5, len(thetaarr) - 4, len(thetaarr) - 3, len(thetaarr) - 2, len(thetaarr) - 1]
        else:
            ts = [thetad - 2, thetad - 1, thetad, thetad + 1, thetad + 2]
        
        cords = []
        for t in ts:
            cords.append([radarr[rid], thetaarr[t], phiarr[phid]])
        id = [get_nearest_point_index(c[0], c[1], c[2], tree) for c in cords]
        data = df.iloc[id].sort_values(by=["r", "theta", "phi"]).drop_duplicates().to_numpy()
        
        
        vals = []
        for i in range(3, 27):
            interp = BarycentricInterpolator(data[:, 1], data[:, i])
            vals.append(interp(thetap))
        vals = np.array(vals).reshape(1, -1)
        #interp = RBFInterpolator(data[:, :3], data[:, 3:], kernel="linear")
        #vals = interp(np.array([[rp, thetap, phip]]).reshape(1, -1))
        return vals
    
    if not rmatch and thetamatch and phimatch:
        rs = []
        # edge cases:
        if rid < 2:
            rs = [0, 1, 2, 3, 4]
        elif rid > len(radarr) - 3:
            rs = [len(radarr) - 5, len(radarr) - 4, len(radarr) - 3, len(radarr) - 2, len(radarr) - 1]
        else:
            rs = [rid - 2, rid - 1, rid, rid + 1, rid + 2]
        
        cords = []
        for r in rs:
            if r == 0:
                cords.append([np.float64(0), np.float64(0), np.float64(0)])
            else:
                cords.append([radarr[r], thetaarr[thetad], phiarr[phid]])
        id = [get_nearest_point_index(c[0], c[1], c[2], tree) for c in cords]
        data = df.iloc[id].sort_values(by=["r", "theta", "phi"]).drop_duplicates().to_numpy()
        
        
        vals = []
        for i in range(3, 27):
            interp = BarycentricInterpolator(data[:, 0], data[:, i])
            vals.append(interp(rp))
        vals = np.array(vals).reshape(1, -1)
        #interp = RBFInterpolator(data[:, :3], data[:, 3:], kernel="linear")
        #vals = interp(np.array([[rp, thetap, phip]]).reshape(1, -1))
        return vals
        
    if rmatch and not thetamatch and not phimatch:
        ts = []
        ps = []
        # edge cases:
        if thetad < 2:
            ts = [0, 1, 2, 3, 4]
        elif thetad > len(thetaarr) - 3:
            ts = [len(thetaarr) - 5, len(thetaarr) - 4, len(thetaarr) - 3, len(thetaarr) - 2, len(thetaarr) - 1]
        else:
            ts = [thetad - 2, thetad - 1, thetad, thetad + 1, thetad + 2]
        if phid < 2:
            ps = [0, 1, 2, 3, 4]
        elif phid > len(phiarr) - 3:
            ps = [len(phiarr) - 5, len(phiarr) - 4, len(phiarr) - 3, len(phiarr) - 2, len(phiarr) - 1]
        else:
            ps = [phid - 2, phid - 1, phid, phid + 1, phid + 2]

        cords = []
        for t in ts:
            for p in ps:
                cords.append([radarr[rid], thetaarr[t], phiarr[p]])
        id = [get_nearest_point_index(c[0], c[1], c[2], tree) for c in cords]
        data = df.iloc[id].sort_values(by=["r", "theta", "phi"]).drop_duplicates().to_numpy()
        
        
        interp = RBFInterpolator(data[:, :3], data[:, 3:], kernel="linear")
        vals = interp(np.array([[rp, thetap, phip]]).reshape(1, -1))
        return vals
        
    if not rmatch and thetamatch and not phimatch:
        rs = []
        ps = []
        # edge cases:
        if rid < 2:
            rs = [0, 1, 2, 3, 4]
        elif rid > len(radarr) - 3:
            rs = [len(radarr) - 5, len(radarr) - 4, len(radarr) - 3, len(radarr) - 2, len(radarr) - 1]
        else:
            rs = [rid - 2, rid - 1, rid, rid + 1, rid + 2]
        if phid < 2:
            ps = [0, 1, 2, 3, 4]
        elif phid > len(phiarr) - 3:
            ps = [len(phiarr) - 5, len(phiarr) - 4, len(phiarr) - 3, len(phiarr) - 2, len(phiarr) - 1]
        else:
            ps = [phid - 2, phid - 1, phid, phid + 1, phid + 2]

        cords = []
        for r in rs:
            for p in ps:
                if r == 0:
                    cords.append([np.float64(0), np.float64(0), np.float64(0)])
                else:
                    cords.append([radarr[r], thetaarr[thetad], phiarr[p]])
        id = [get_nearest_point_index(c[0], c[1], c[2], tree) for c in cords]
        data = df.iloc[id].sort_values(by=["r", "theta", "phi"]).drop_duplicates().to_numpy()
        
        
        interp = RBFInterpolator(data[:, :3], data[:, 3:], kernel="linear")
        vals = interp(np.array([[rp, thetap, phip]]).reshape(1, -1))
        return vals
        
    if not rmatch and not thetamatch and phimatch:
        rs = []
        ts = []
        # edge cases:
        if rid < 2:
            rs = [0, 1, 2, 3, 4]
        elif rid > len(radarr) - 3:
            rs = [len(radarr) - 5, len(radarr) - 4, len(radarr) - 3, len(radarr) - 2, len(radarr) - 1]
        else:
            rs = [rid - 2, rid - 1, rid, rid + 1, rid + 2]
        if thetad < 2:
            ts = [0, 1, 2, 3, 4]
        elif thetad > len(thetaarr) - 3:
            ts = [len(thetaarr) - 5, len(thetaarr) - 4, len(thetaarr) - 3, len(thetaarr) - 2, len(thetaarr) - 1]
        else:
            ts = [thetad - 2, thetad - 1, thetad, thetad + 1, thetad + 2]

        cords = []
        for r in rs:
            for t in ts:
                if r == 0:
                    cords.append([np.float64(0), np.float64(0), np.float64(0)])
                else:
                    cords.append([radarr[r], thetaarr[t], phiarr[phid]])
        id = [get_nearest_point_index(c[0], c[1], c[2], tree) for c in cords]
        data = df.iloc[id].sort_values(by=["r", "theta", "phi"]).drop_duplicates().to_numpy()
        
        
        interp = RBFInterpolator(data[:, :3], data[:, 3:], kernel="linear")
        vals = interp(np.array([[rp, thetap, phip]]).reshape(1, -1))
        return vals
        
    else:
        # if no rmatch, no thetamatch, no phimatch
        rs = []
        ps = []
        ts = []
        # edge cases:
        if rid < 2:
            rs = [0, 1, 2, 3, 4]
        elif rid > len(radarr) - 3:
            rs = [len(radarr) - 5, len(radarr) - 4, len(radarr) - 3, len(radarr) - 2, len(radarr) - 1]
        else:
            rs = [rid - 2, rid - 1, rid, rid + 1, rid + 2]
        
        if phid < 2:
            ps = [0, 1, 2, 3, 4]
        elif phid > len(phiarr) - 3:
            ps = [len(phiarr) - 5, len(phiarr) - 4, len(phiarr) - 3, len(phiarr) - 2, len(phiarr) - 1]
        else:
            ps = [phid - 2, phid - 1, phid, phid + 1, phid + 2]

        if thetad < 2:
            ts = [0, 1, 2, 3, 4]
        elif thetad > len(thetaarr) - 3:
            ts = [len(thetaarr) - 5, len(thetaarr) - 4, len(thetaarr) - 3, len(thetaarr) - 2, len(thetaarr) - 1]
        else:
            ts = [thetad - 2, thetad - 1, thetad, thetad + 1, thetad + 2]



        cords = []
        for r in rs:
            for p in ps:
                for t in ts:
                    if r == 0:
                        cords.append([np.float64(0), np.float64(0), np.float64(0)])
                    else:
                        cords.append([radarr[r], thetaarr[t], phiarr[p]])
        id = [get_nearest_point_index(c[0], c[1], c[2], tree) for c in cords]
        data = df.iloc[id].sort_values(by=["r", "theta", "phi"]).drop_duplicates().to_numpy()
       
        
        interp = RBFInterpolator(data[:, :3], data[:, 3:], kernel="thin_plate_spline")
        vals = interp(np.array([[rp, thetap, phip]]).reshape(1, -1))
        return vals


def process_line(line):
    parts = line.split()[1:]  # Ignore the first column
    new_grid = list(map(np.float64, parts))
    output_file = f"CTS_bin-proc{parts[-1]}.d"     
    return new_grid, output_file

def gridmaker(gridparams):
    xmin, ymin, zmin, dx, dy, dz, nx, ny, nz, _ = gridparams
    xarr = np.arange(xmin, xmin + nx*dx, dx)
    yarr = np.arange(ymin, ymin + ny*dy, dy)
    zarr = np.arange(zmin, zmin + nz*dz, dz)
    gridcoords = np.array(np.meshgrid(xarr, yarr, zarr)).T.reshape(-1, 3)
    return gridcoords, nx, ny, nz

def cartesiantospherical(arr):
    r = np.sqrt(arr[:, 0]**2 + arr[:, 1]**2 + arr[:, 2]**2)
    theta = np.where(r != 0.0, np.arccos(np.clip(arr[:, 2] / r, -1.0, 1.0)), 0.0)
    phi = np.arctan2(arr[:, 1], arr[:, 0])
    new_grid = np.column_stack((r, theta, phi))
    new_grid = np.where(abs(new_grid) < 1e-14, np.float64(0), new_grid)
    return new_grid

def routine(line):
    gridparams, outputfile = process_line(line)
    cartgrid, Nx, Ny, Nz = gridmaker(gridparams)
    sphercalgrid = cartesiantospherical(cartgrid)
    interpolated_data = [
                interpolationpoint(rad_arr, theta_arr, phi_arr, point, tree, df_3d)
                for point in sphercalgrid
            ]

    interpolated_data = np.array(interpolated_data).reshape(-1, 24)
    interpolated_data = np.column_stack((cartgrid, interpolated_data))

    ntot = Nx * Ny * Nz
    
    with open(outputfile, "wb") as f:
        f.write(np.array([Nx, Ny, Nz, ntot], dtype=np.float64).tobytes())
        f.write(interpolated_data.astype(np.float64).tobytes())
    return

with open("/data/sjammi6/thesisproject/grids_bh_disk_patrik", "r") as f:
    for line in tqdm(f):
        routine(line)
print("FINISHED INTERPOLATION :)")

