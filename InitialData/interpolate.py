import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import pandas as pd

#### BH Horizon Radius r_s = 0.5 * np.sqrt(kerrm**2 - kerra**2)
kerrm = 1
kerra = 0.8
r_s = 0.5 * np.sqrt(kerrm**2 - kerra**2)

def locate_point(r_arr,theta_arr,phi_arr, p): #r_arr, theta_arr, phi_arr define the ID grid, p defines a point in the new grid
    # p should be a tuple (x,y,z) in cartesian coords
    x,y,z = p
    r_0 = np.sqrt(x**2 + y**2 + z**2)
    
    theta_0 = np.arctan2(np.sqrt(x**2+y**2),z)
    phi_0 = np.arctan2(y,x)
    phi_0 = phi_0 % (2 * np.pi) #ensures phi takes values between [0,2pi]
    theta_0 = theta_0 % (np.pi) #ensures theta takes values between [0, phi]
    rs, ts, ps = None, None, None
    for i in range(len(r_arr)-1):
        if r_0>=r_arr[i] and r_0<=r_arr[i+1]:
            rs = (r_arr[i], r_arr[i+1])
            break
    
    for i in range(len(theta_arr)-1):
        if theta_0>=theta_arr[i] and theta_0<=theta_arr[i+1]:
            ts = (theta_arr[i], theta_arr[i+1])
            break
    
    for i in range(len(phi_arr)-1):
        if phi_0>=phi_arr[i] and phi_0<=phi_arr[i+1]:
            ps = (phi_arr[i], phi_arr[i+1])
            break

    if rs is None or ts is None or ps is None:
        print(f'p={p} is not in our grid domain')
        return 
    
    else:
        combinations = []
        
        if rs[0] == 0:
            for t in ts:
                for p in ps:
                    combinations.append([rs[1],t,p])
            combinations.append([0,0,0])
        else:
            for r in rs:
                for t in ts:
                    for p in ps:
                        ### save points in spherical
                        combinations.append([r,t,p])

        return combinations
    

########## For testing purposes
# # Assuming the file is in the same directory as your script
# file_path = 'extrisic_curvature_vars.dat'

# # Load data from the file
# data = np.loadtxt(file_path)

# # Extract the values from the 3rd and 4th columns
# r_arr = data[:, 2]  # 3rd column, Python indexing starts from 0
# t_arr = data[:, 3]  # 4th column

# # Remove all duplicates and sort in ascending order
# r_arr = np.sort(np.unique(r_arr))
# theta_arr = np.sort(np.unique(t_arr))

# num_phi = 500
# phi_arr = np.linspace(0, 2*np.pi, num_phi+1)
########## 

def plot_points(r_arr,theta_arr,phi_arr, p):

    # p = (5,3,10)
    # Define the list of 8 points

    points = locate_point(r_arr,theta_arr,phi_arr, p)

    print(points)
    # Extract x, y, z coordinates for both points and p
    x_points, y_points, z_points = zip(*points)
    x_p, y_p, z_p = p

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the 8 points in blue
    ax.scatter(x_points, y_points, z_points, c='blue', label='Points')

    # Plot the single point p in red
    ax.scatter(x_p, y_p, z_p, c='red', label='Point p')

    #Change view and limits of plot
    ax.view_init(elev=80, azim=0)  

    # ax.set_xlim(0, 5)  # Set limits for the X-axis
    # ax.set_ylim(-5, 5)  # Set limits for the Y-axis
    # ax.set_zlim(-5, 5)  # Set limits for the Z-axis

    # Add labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Add a legend
    ax.legend()

    plt.savefig('points_plot.png')

    # Show the plot
    # plt.show()

def calc_simple_average(points, scalar_values): #each point should be a numpy array with (x,y,z), scalar values is a numpy array
    return np.mean(scalar_values)


def calc_weighted_average(points, scalar_values, p): #points should be a numpy array
    # Compute distances from each point to p
    distances = np.linalg.norm(points - np.array(p), axis=1)
    
    # Compute weights (inverse distances)
    weights = 1.0 / distances
    
    # Compute weighted sum of scalar values
    weighted_sum = np.sum(weights * scalar_values)
    
    # Normalize by dividing by the sum of weights
    normalized_weighted_average = weighted_sum / np.sum(weights)
    
    return normalized_weighted_average


def load_3d_data(data_3d_path):
    data = np.loadtxt(data_3d_path)
    return data


def spherical2cart(r,t,p):
    x = r*np.sin(t)*np.cos(p)
    y = r*np.sin(t)*np.sin(p)
    z = r*np.cos(t)
    return (x,y,z)

# data = load_3d_data('3D_data/ext_cur.3d')
def interpolate_point(p, df, rad_arr, theta_arr, phi_arr, idx_point_map):
    num_scalars = len(df.columns) - 3
    surrounding_pts = locate_point(rad_arr,theta_arr,phi_arr, p)

    if not surrounding_pts:
        print(f'could not find point in grid, grid has these properties:')
        print(f'r_min, r_max, num_r = {rad_arr[0], rad_arr[-1], len(rad_arr)}, t_min, t_max, num_t = {theta_arr[0], theta_arr[-1], len(theta_arr)}, p_min, p_max, num_p = {phi_arr[0], phi_arr[-1], len(phi_arr)}')
        return ValueError
    surrounding_pts = np.array(surrounding_pts)

    points = []
    scalars = []

    for point in surrounding_pts:
        r, theta, phi = point
        # Find the closest matching row in df based on spherical coordinates
        row_idx = idx_point_map[(r, theta, phi)]
        closest_row = df.iloc[row_idx]
        # Check if any rows are found
        if not closest_row.empty:
            point_array = closest_row
            x,y,z = spherical2cart(point_array.iloc[0], point_array.iloc[1], point_array.iloc[2])
            points.append([x,y,z])
            scalars.append(point_array[3:])
    points = np.array(points)
    scalars = np.array(scalars)
    interpolated_vals = []
    for i in range(num_scalars):
        interpolated_vals.append(calc_weighted_average(points, scalars[:,i], p))
    return np.array(interpolated_vals)

### For testing purposes
# s1 = time.time()
# print('starting preprocessing')
# df = pd.read_hdf('3D_data/all_data.h5', key='df')
# rad_arr = np.sort(df.r.unique())
# theta_arr = np.sort(df.theta.unique())
# phi_arr = np.sort(df.phi.unique())
# idx_point_map = {(r, theta, phi): index for index, (r, theta, phi) in enumerate(zip(df['r'], df['theta'], df['phi']))}

# print(f'preprocessing data took {time.time() - s1} sec')
# p = (-0.099, -0.2875, -0.025)
# p1 = (-0.249, 0.1625, -0.05)
# print(interpolate_point(p, df, rad_arr, theta_arr, phi_arr, idx_point_map))
# print(interpolate_point(p1, df, rad_arr, theta_arr, phi_arr, idx_point_map))



