import numpy as np
import pandas as pd
from numba import njit, prange


#### BH Horizon Radius r_s = 0.5 * np.sqrt(kerrm**2 - kerra**2)
kerrm = 1
kerra = 0.8
r_s = 0.5 * np.sqrt(kerrm**2 - kerra**2)


@njit(parallel=True)
def locate_point(
    r_arr, theta_arr, phi_arr, p
):  # r_arr, theta_arr, phi_arr define the ID grid, p defines a point in the new grid
    # p should be a tuple (x,y,z) in cartesian coords
    phi_arr = np.append(phi_arr[-2], np.append(phi_arr, phi_arr[1]))
    x, y, z = p
    r_0 = np.sqrt(x**2 + y**2 + z**2)
    theta_0 = np.arccos(z / r_0)
    skip = False
    phi_0 = np.arctan2(y, x)
    phi_0 = phi_0 % (2 * np.pi)  # ensures phi takes values between [0,2pi]
    rs, ts, ps = [], [], []

    for i in prange(1, len(phi_arr) - 1):
        if phi_0 == phi_arr[i]:
            ps = [phi_arr[i - 1], phi_arr[i + 1]]
            break
        if phi_arr[i] < phi_0 < phi_arr[i + 1]:
            ps = [phi_arr[i], phi_arr[i + 1]]
            break
    ps = np.array(ps)

    if theta_0 == theta_arr[0]:
        ts = [theta_arr[1]]
        skip = True
        ps = np.append(ps, [phi_arr[-2] - ps[0], phi_arr[-2] - ps[1]])
    if theta_0 == theta_arr[-1]:
        ts = [theta_arr[-2]]
        ps = np.append(ps, [phi_arr[-2] - ps[0], phi_arr[-2] - ps[1]])
        skip = True
    if not skip:
        for i in prange(len(theta_arr)):
            if theta_0 == theta_arr[i]:
                ts = (theta_arr[i - 1], theta_arr[i + 1])
                break
            if theta_arr[i] < theta_0 < theta_arr[i + 1]:
                ts = (theta_arr[i], theta_arr[i + 1])
                break
    ts = np.array(ts)
    for i in prange(len(r_arr)):
        if r_0 == r_arr[i]:
            if i == 0:
                rs = [r_arr[i], r_arr[1]]
                break
            else:
                rs = (r_arr[i - 1], r_arr[i + 1])
                break
        if r_arr[i] < r_0 < r_arr[i + 1]:
            rs = (r_arr[i], r_arr[i + 1])
            break
    rs = np.array(rs)
    combinations = []
    if rs[0] == r_arr[0]:
        combinations.append([np.float64(0), np.float64(0), np.float64(0)])
        for tmpt in ts:
            for tmpp in ps:
                combinations.append([rs[1], tmpt, tmpp])
    else:
        for tmpr in rs:
            for tmpt in ts:
                for tmpp in ps:
                    combinations.append([tmpr, tmpt, tmpp])
    return np.array(combinations)


@njit(parallel=True)
def calc_weighted_average(points, scalar_vals, p):
    d = np.sum((points - p) ** (1 / 2), axis=1)
    if any(d <= 1e-7):
        return scalar_vals[np.argmin(d)]
    weights = 1 / d
    weights = weights.reshape((-1, 1))
    return np.sum(weights * scalar_vals, axis=0) / np.sum(weights)


@njit(parallel=True)
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
