import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from numba import njit


#### BH Horizon Radius r_s = 0.5 * np.sqrt(kerrm**2 - kerra**2)
kerrm = 1
kerra = 0.8
r_s = 0.5 * np.sqrt(kerrm**2 - kerra**2)


def locate_point(
    r_arr, theta_arr, phi_arr, p
):  # r_arr, theta_arr, phi_arr define the ID grid, p defines a point in the new grid
    # p should be a tuple (x,y,z) in cartesian coords
    phi_arr_extended = np.zeros(len(phi_arr) + 2, dtype=np.float64)
    phi_arr_extended[0] = phi_arr[-2]
    phi_arr_extended[1:-1] = phi_arr
    phi_arr_extended[-1] = phi_arr[1]    
    x, y, z = p
    r_0 = np.sqrt(x**2 + y**2 + z**2)
    theta_0 = np.arccos(z / r_0)
    skip = False
    origin = False
    phi_0 = np.arctan2(y, x)
    phi_0 = phi_0 % (2 * np.pi)  # ensures phi takes values between [0,2pi]
    ts = np.zeros(2, dtype=np.float64)
    ps = np.zeros(4, dtype=np.float64)
    rs = np.zeros(2, dtype=np.float64)
    ep = np.array([phi_arr[-2], phi_arr[1]])

    for i in range(1, len(phi_arr_extended) - 1):
        if phi_0 == phi_arr_extended[i]:
            ps[0] = phi_arr_extended[i - 1]
            ps[1] = phi_arr_extended[i + 1]
            break
        if phi_arr_extended[i] < phi_0 < phi_arr_extended[i + 1]:
            ps[0] = phi_arr_extended[i]
            ps[1] = phi_arr_extended[i + 1]
            break

    if theta_0 == theta_arr[0]:
        ts[0] = theta_arr[1]
        skip = True
    if theta_0 == theta_arr[-1]:
        ts[0] = theta_arr[-2]
        skip = True
    if skip:
        ps[2] = ep[0] - ps[0]
        ps[3] = ep[1] - ps[1]

    if not skip:
        for i in range(len(theta_arr)):
            if theta_0 == theta_arr[i]:
                ts[0] = theta_arr[i - 1]
                ts[1] = theta_arr[i + 1] 
                break
            if theta_arr[i] < theta_0 < theta_arr[i + 1]:
                ts[0] = theta_arr[i]
                ts[1] = theta_arr[i + 1]
                break

    for i in range(len(r_arr)):
        if r_0 == r_arr[i]:
            if i == 0:
                rs[0] = r_arr[i]
                rs[1] = r_arr[i + 1]
                origin = True
                break
            else:
                rs[0] = r_arr[i - 1]
                rs[1] = r_arr[i + 1]
                break
        if r_arr[i] < r_0 < r_arr[i + 1]:
            rs[0] = r_arr[i]
            rs[1] = r_arr[i + 1]
            break
    tc = 8
    if origin:
        tc = 5
    combinations = np.zeros((tc, 3), dtype=np.float64)

    if origin:
        combinations[0] = np.array([np.float64(0), np.float64(0), np.float64(0)])
        if skip:
            combinations[1] = np.array([rs[1], ts[0], ps[0]])
            combinations[2] = np.array([rs[1], ts[0], ps[1]])
            combinations[3] = np.array([rs[1], ts[0], ps[2]])
            combinations[4] = np.array([rs[1], ts[0], ps[3]])
        else:
            combinations[1] = np.array([rs[1], ts[0], ps[0]])
            combinations[2] = np.array([rs[1], ts[0], ps[1]])
            combinations[3] = np.array([rs[1], ts[1], ps[0]])
            combinations[4] = np.array([rs[1], ts[1], ps[1]])
    else:
        if skip:
            combinations[0] = np.array([rs[0], ts[0], ps[0]])
            combinations[1] = np.array([rs[0], ts[0], ps[1]])
            combinations[2] = np.array([rs[0], ts[0], ps[2]])
            combinations[3] = np.array([rs[0], ts[0], ps[3]])
            combinations[4] = np.array([rs[1], ts[0], ps[0]])
            combinations[5] = np.array([rs[1], ts[0], ps[1]])
            combinations[6] = np.array([rs[1], ts[0], ps[2]])
            combinations[7] = np.array([rs[1], ts[0], ps[3]])
        else:
            combinations[0] = np.array([rs[0], ts[0], ps[0]])
            combinations[1] = np.array([rs[0], ts[0], ps[1]])
            combinations[2] = np.array([rs[0], ts[1], ps[0]])
            combinations[3] = np.array([rs[0], ts[1], ps[1]])
            combinations[4] = np.array([rs[1], ts[0], ps[0]])
            combinations[5] = np.array([rs[1], ts[0], ps[1]])
            combinations[6] = np.array([rs[1], ts[1], ps[0]])
            combinations[7] = np.array([rs[1], ts[1], ps[1]])

    return combinations 


def calc_weighted_average(points, scalar_vals, p):
    d = np.sum(np.abs((points - p)) ** (1/2), axis=1)
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
