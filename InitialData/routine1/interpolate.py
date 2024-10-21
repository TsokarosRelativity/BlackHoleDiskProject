import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import pandas as pd

#### BH Horizon Radius r_s = 0.5 * np.sqrt(kerrm**2 - kerra**2)
kerrm = 1
kerra = 0.8
r_s = 0.5 * np.sqrt(kerrm**2 - kerra**2)


def locate_point(
    r_arr, theta_arr, phi_arr, p
):  # r_arr, theta_arr, phi_arr define the ID grid, p defines a point in the new grid
    # p should be a tuple (x,y,z) in cartesian coords
    x, y, z = p
    r_0 = np.sqrt(x**2 + y**2 + z**2)
    theta_0 = np.arccos(z / r_0)
    phi_0 = np.arctan2(y, x)
    phi_0 = phi_0 % (2 * np.pi)  # ensures phi takes values between [0,2pi]
    rs, ts, ps = None, None, None

    for i in range(len(r_arr) - 1):
        if r_arr[i] <= r_0 <= r_arr[i + 1]:
            rs = (r_arr[i], r_arr[i + 1])
            break
    for i in range(len(theta_arr) - 1):
        if theta_arr[i] <= theta_0 <= theta_arr[i + 1]:
            ts = (theta_arr[i], theta_arr[i + 1])
            break
    for i in range(len(phi_arr) - 1):
        if phi_arr[i] <= phi_0 <= phi_arr[i + 1]:
            ps = (phi_arr[i], phi_arr[i + 1])
            break

    if rs is None or ts is None or ps is None:
        return None

    combinations = []
    if rs[0] == np.float64(0):
        combinations.append([np.float64(0), np.float64(0), np.float64(0)])
        for t in ts:
            for p in ps:
                combinations.append([rs[1], t, p])
    else:
        for r in rs:
            for t in ts:
                for p in ps:
                    combinations.append([r, t, p])

    return np.array(combinations)


def calc_weighted_average(points, scalar_vals, p):
    d = np.sum((points - p) ** 2, axis=1)
    if any(d <= 1e-4):
        return scalar_vals[np.argmin(d)]
    weights = 1 / ((d) ** 2)
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
