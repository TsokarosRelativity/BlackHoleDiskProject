import numpy as np
import pandas as pd


#### BH Horizon Radius r_s = 0.5 * np.sqrt(kerrm**2 - kerra**2)
kerrm = 1
kerra = 0.8
r_s = 0.5 * np.sqrt(kerrm**2 - kerra**2)
def locate_point(r_arr, theta_arr, phi_arr, p):
    x, y, z = p
    r_0 = np.sqrt(x**2 + y**2 + z**2)
    theta_0 = np.arccos(z / r_0) if r_0 != 0 else np.float64(0)
    phi_0 = np.arctan2(y, x) % (2 * np.pi)

    # Find nearest indices
    r_idx_1 = np.searchsorted(r_arr, r_0) - 1
    r_idx_2 = np.searchsorted(r_arr, r_0, side="right")
    theta_idx_1 = np.searchsorted(theta_arr, theta_0) - 1
    theta_idx_2 = np.searchsorted(theta_arr, theta_0, side="right")
    phi_idx_1 = np.searchsorted(phi_arr, phi_0) - 1
    phi_idx_2 = np.searchsorted(phi_arr, phi_0, side="right")

    theta_edge_case = False
    if theta_idx_1 == -1:
        theta_idx_1 = 1
        theta_edge_case = True
    if theta_idx_2 == len(theta_arr):
        theta_idx_2 = theta_idx_1
        theta_edge_case = True
    if phi_idx_1 == -1:
        phi_idx_1 = len(phi_arr) - 1
    if phi_idx_2 == len(phi_arr):
        phi_idx_2 = 0
    if phi_idx_1 == len(phi_arr):
        phi_idx_1 = len(phi_arr) - 1
    # Handle edge cases
    # Get surrounding points
    rs = [r_arr[r_idx_1], r_arr[r_idx_2]]
    ts = [theta_arr[theta_idx_1], theta_arr[theta_idx_2]]
    ps = [phi_arr[phi_idx_1], phi_arr[phi_idx_2]]
    if theta_edge_case:
        phi_idx_3 = (phi_idx_1 + len(phi_arr) // 2) % len(phi_arr)
        phi_idx_4 = (phi_idx_2 + len(phi_arr) // 2) % len(phi_arr)
        ps = [phi_arr[phi_idx_3], phi_arr[phi_idx_4]]
    # Generate all combinations
    combinations = np.array(np.meshgrid(rs, ts, ps)).T.reshape(-1, 3)
    z = combinations[:, 0] == np.float64(0)
    combinations[z] = np.array([np.float64(0), np.float64(0), np.float64(0)])
    combinations = np.unique(combinations, axis=0)
    return combinations
def calc_weighted_average(points, scalar_vals, p):
    d = np.sum(np.abs((points - p)), axis=1)
    weights = 1 / d
    weights = weights.reshape((-1, 1))
    return np.sum(weights * scalar_vals, axis=0) / np.sum(weights)


def spherical2cart(r, t, p):
    x = r * np.sin(t) * np.cos(p)
    y = r * np.sin(t) * np.sin(p)
    z = r * np.cos(t)
    return (x, y, z)


def interpolate_point(p, df, rad_arr, theta_arr, phi_arr, idx_point_map):
    surroundingpoints = locate_point(rad_arr, theta_arr, phi_arr, p)  # numpy array
    carpointdata = []
    scalardata = []
    for point in surroundingpoints:
        tmpr, tmpt, tmpp = point
        rowIndex = idx_point_map[(tmpr, tmpt, tmpp)]
        rowdata = df.iloc[rowIndex]
        tmpx, tmpy, tmpz = spherical2cart(tmpr, tmpt, tmpp)
        carpointdata.append([tmpx, tmpy, tmpz])
        scalardata.append(rowdata[3:])
    carpointdata = np.array(carpointdata)
    scalardata = np.array(scalardata)
    interpolateddata = calc_weighted_average(carpointdata, scalardata, p)
    return np.hstack((p, interpolateddata))
