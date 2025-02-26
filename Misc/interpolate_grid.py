import numpy as np
import pandas as pd
from interpolate import *
import time



data_folder = "September_20_2024_16-53/" 
print('STARTING WOOHOO')
output_dir = "processed_grids/" ##end with /
#output_file = 'test_data.3d'

print('LOADING IN DATA')
s1 = time.time()
df = pd.read_hdf(data_folder + '3D_data/all_data_updated_jacobian.h5', key='df')
print(f'loaded data in {time.time() - s1} sec')

print('PREPROCESSING DATA')
s2 = time.time()
rad_arr = np.sort(df.r.unique())
theta_arr = np.sort(df.theta.unique())
phi_arr = np.sort(df.phi.unique())
idx_point_map = {(r, theta, phi): index for index, (r, theta, phi) in enumerate(zip(df['r'], df['theta'], df['phi']))}
print(f'processed data in {time.time() - s2} sec')
def interpolate_grid(new_grid, df):
    ### generate new grid
    x_min, y_min, z_min, dx, dy, dz, Nx, Ny, Nz, MPI_ID = new_grid
    x_arr = np.arange(x_min, x_min + dx*Nx, dx)
    y_arr = np.arange(y_min, y_min + dy*Ny, dy)
    z_arr = np.arange(z_min, z_min + dz*Nz, dz)
    s2 = time.time()
    data = []
    total_points = Nx*Ny*Nz
    processed_points = 0
    for z in z_arr:
        for y in y_arr:
            for x in x_arr:
                processed_points += 1
                # progress = processed_points / total_points * 100
                # print(f'Progress: {progress:.2f}%', end='\r')
                p = np.array([x,y,z])
                # new_line = [x, y, z]
                new_line = interpolate_point(p, df, rad_arr, theta_arr, phi_arr, idx_point_map)
                data.append(new_line)
    return np.concatenate(data)

def process_line(line):
    parts = line.split()[1:]  # Ignore the first column
    new_grid = list(map(float, parts))  
    output_file = f"CTS_bin-proc{parts[-1]}.d"  # Use the last value as the name of the output file
    return new_grid, output_file

def interpolate_grids(input_file, df):
    with open(input_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            start = time.time()
            new_grid, output_file = process_line(line)
            interpolated_grid = interpolate_grid(new_grid, df)
            # np.savetxt(output_dir + output_file, interpolated_grid)
            nx, ny, nz = new_grid[6:9]
            ntot = nx*ny*nz
            with open(data_folder + output_dir + output_file, 'wb') as f:
                f.write(np.array([nx, ny, nz, ntot], dtype=np.int32).tobytes())
                f.write(interpolated_grid.astype(np.float64).tobytes())
            print(f"File '{output_file}' created in {time.time() - start} seconds")


s3 = time.time()

gridfile = "grids_bh_disk_patrik"
interpolate_grids(gridfile, df)

### For test purposes
# input_file = 'patryk_grids/test.txt'
# interpolate_grids(input_file, df)
