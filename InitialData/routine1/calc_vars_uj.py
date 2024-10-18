import numpy as np
import time
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    prog="Calc Vars: Updated Jacobian",
    description="2D -> 3D Initial Data: Calculate Variables Routine",
    epilog="Routine takes in 2D initial data and converts 2D spherical initial data into 3D cartesian initial data",
)
parser.add_argument(
    "folder",
    type=str,
    help="args.folder in which the Initial Data is being stored and calculated",
)
parser.add_argument(
    "--v", action="store_true", help="showing details while running routine"
)
args = parser.parse_args()


start_time = time.time()
if args.v:
    print("2D data -> 3D data: Calculate Variables Routine")


######### INPUTS ###########################
file_path = "2D_data/all_vars.dat"
rs = np.float64(0.3)
n_phi = 100

K = 5.3206252168330817e-002  # this is K from the polytropic equation P = K*(rho**Gamma)... K can be found in x.dat as k
gamma = 4 / 3  # can also be found in x.dat
G = 6.674 * 10 ** (-11)  # gravitational constant in SI units
c = 3 * 10 ** (8)  # speed of light in SI units
epsilon_0 = 8.854 * 10 ** (-12)  # permitivity of free space in SI units
unit_conv = (
    G ** (1 / 2) * c ** (-1) * epsilon_0 ** (-1 / 2) * 10 ** (4)
)  # conversion from geometrized units for magnetic B field to Gauss (CGS units)

######################################################################################################################
if args.v:
    print("         args.folder: {}".format(args.folder))
    print("         file_path: {}".format(file_path))
    print("         Number of Phi Values: {}".format(n_phi))
    print("         Horizon: {}".format(rs))
    print("         K (constant of proportionality from polytropic EOS): {}".format(K))
    print("         gamma: {}".format(gamma))
    print("         G (gravitational constant): {}".format(G))
    print("         c (speed of light): {}".format(c))
    print("         epsilon_0 (permitivity of free space): {}".format(epsilon_0))


#### READING IN DATA
s1 = time.time()
column_names = [
    "x",
    "z",
    "r",
    "theta",
    "alpha",
    "psi",
    "q",
    "betak",
    "betat",
    "dbetatdr",
    "dbetatdt",
    "He",
    "Hf",
    "Omega",
    "b2",
    "rho",
]
df = pd.read_csv(
    args.folder + file_path, delimiter=r"\s+", header=None, names=column_names
)
s2 = time.time()
if args.v:
    print(f"read in the data in {s2-s1} secs")

tolerance = 0.00001
# df = df[~(np.abs(df.r - rs) < tolerance)]

# Convert all negative values of b2 to 0
df["b2"] = df["b2"].apply(lambda x: max(x, 0))


# isnt this section of code reflecting the coordinates across the z axis, giving -x values not - z values

#### TURNING 2D DATA INTO 3D DATA
# Reflect over x-y plane (-z values)
df_temp = df.copy()
df_temp["theta"] = np.pi - df_temp["theta"]
df = pd.concat([df, df_temp]).reset_index(drop=True)
df = df.drop_duplicates()
s3 = time.time()

if args.v:
    print(f"reflected to -z in {s3-s2} secs")

    print("-----------------------------------------------")
# Calculate evenly spaced phi values
df = pd.DataFrame(np.repeat(df.values, n_phi, axis=0), columns=df.columns)
num_rows = len(df)
phi_values = np.linspace(0, 2 * np.pi, num_rows % n_phi or n_phi)
df["phi"] = np.tile(phi_values, num_rows // n_phi + 1)[:num_rows]
s4 = time.time()
if args.v:

    print(f"added phi values in {s4-s3} secs")

    print("-----------------------------------------------")
#### PERFORMING CALCULATIONS
s5 = time.time()
df["x"] = df.r * np.sin(df.theta) * np.cos(df.phi)
df["y"] = df.r * np.sin(df.theta) * np.sin(df.phi)
df["z"] = df.r * np.cos(df.theta)
df["beta"] = df.betat + df.betak
df["P"] = K * (df.rho**gamma)

# Assemble Jacobian (3x3)
df["J_00"] = df.x / df.r
df["J_01"] = df.y / df.r
df["J_02"] = df.z / df.r

df["J_10"] = np.cos(df.phi) * df.z / df.r**2
df["J_11"] = np.sin(df.phi) * df.z / df.r**2
df["J_12"] = -np.sqrt(df.x**2 + df.y**2) / df.r**2

df["J_20"] = -np.sin(df.phi) / (df.r * np.sin(df.theta))
df["J_21"] = np.cos(df.phi) / (df.r * np.sin(df.theta))
df["J_22"] = 0

df["zero"] = np.zeros(len(df), dtype=np.float64)

# Only nonzero components of beta
df["beta_phi"] = (df.psi**4) * (df.r**2) * (np.sin(df.theta)) * df.beta
df["beta__phi"] = df.beta

df["beta__x"] = df.J_02 * df.beta__phi
df["beta__y"] = df.J_12 * df.beta__phi
df["beta__z"] = df.J_22 * df.beta__phi

# Defining K (extrinsic curvature) in different regions
df["K_rp"] = np.zeros(len(df), dtype=np.float64)
df["K_tp"] = np.zeros(len(df), dtype=np.float64)

cond1 = np.abs(df["r"] - rs) < 0.00001  # on horizon
cond2 = df.z >= 0  # positive z

df.loc[cond1, "K_rp"] = (df.He * np.sin(df.theta) ** 2) / (df.psi**2 * df.r**2)
df.loc[cond1, "K_tp"] = (df.Hf * np.sin(df.theta)) / (df.psi**2 * df.r)
df.loc[~cond1, "K_rp"] = (df.He * np.sin(df.theta) ** 2) / (df.psi**2 * df.r**2) + (
    2 * df.alpha
) ** (-1) * df.psi**4 * df.r**2 * np.sin(df.theta) ** 2 * df.dbetatdr
df.loc[(~cond1) & (cond2), "K_tp"] = (df.Hf * np.sin(df.theta)) / (
    df.psi**2 * df.r
) + (2 * df.alpha) ** (-1) * df.psi**4 * df.r**2 * np.sin(
    df.theta
) ** 2 * df.dbetatdt
df.loc[(~cond1) & (~cond2), "K_tp"] = (df.Hf * np.sin(df.theta)) / (
    df.psi**2 * df.r
) - (2 * df.alpha) ** (-1) * df.psi**4 * df.r**2 * np.sin(
    df.theta
) ** 2 * df.dbetatdt

# Simplifying K expressions
df["K_rp"] = df["K_rp"].apply(lambda x: x if not np.isnan(x) else 0)
df["K_tp"] = df["K_tp"].apply(lambda x: x if not np.isnan(x) else 0)

# Assemble conformal 3-metric in spherical coords (covariant form)
df["gamma_00"] = df.psi**4 * np.exp(2 * df.q)  # gamma_rr
df["gamma_11"] = df.psi**4 * np.exp(2 * df.q) * df.r**2  # gamma_thetatheta
df["gamma_22"] = df.psi**4 * df.r**2 * np.sin(df.theta) ** 2  # gamma_phiphi

# Convert gamma to covariant form
df["gamma__00"] = 1 / df["gamma_00"]
df["gamma__11"] = 1 / df["gamma_11"]
df["gamma__22"] = 1 / df["gamma_22"]

# Calc Fluid Velocity (These are upper indices indicated by two lower '__')
df["u__t"] = np.zeros(len(df), dtype=np.float64)
df.loc[(df.rho == 0), "u__t"] = 0.0
df.loc[~(df.rho == 0), "u__t"] = (
    df.alpha**2 - df.gamma_22 * (df.Omega + df.beta) ** 2
) ** (-1 / 2)
df["u__phi"] = df.Omega * df.u__t
df["u__x"] = df.J_02 * df.u__phi
df["u__y"] = df.J_12 * df.u__phi
df["u__z"] = np.zeros(len(df), dtype=np.float64)

# Calc Mag Field Components (These are lower indices)
df["b2_phi"] = df.u__t**2 * df.gamma_22 * df.alpha**2 * df.b2
df["b__phi"] = df.gamma__22 * df.b2_phi ** (
    1 / 2
)  # contravariant form of the phi component of b (the only nonzero component)
df["b_t"] = -df.Omega * df.b2_phi ** (1 / 2)

# Convert b__phi to b_phi using gamma_22
df["b_phi"] = df.b__phi * df.gamma_22  # Convert upper index to lower index

# Convert Units
G = 6.674 * 10 ** (-11)  # gravitational constant in SI units
c = 3 * 10 ** (8)  # speed of light in SI units
epsilon_0 = 8.854 * 10 ** (-12)  # permitivity of free space in SI units
unit_conv = (
    G ** (1 / 2) * c ** (-1) * epsilon_0 ** (-1 / 2) * 10 ** (4)
)  # conversion from geometrized units for magnetic B field to Gauss (CGS units)

df["B__x"] = (
    (4 * np.pi) ** (1 / 2) * df.J_02 * df.b__phi * unit_conv
)  # contravariant mag field in CGI (gauss) units
df["B__y"] = (4 * np.pi) ** (1 / 2) * df.J_12 * df.b__phi * unit_conv
df["B__z"] = np.zeros(len(df), dtype=np.float64)
s6 = time.time()

if args.v:
    print(f"calculations done in {s6-s5} secs")

    print("-----------------------------------------------")
s7 = time.time()
# Explicitly calculating gamma_ab in Cartesian coordinates
df["gxx"] = (df.psi**4 * np.exp(2 * df.q)) * (
    ((df.x / df.r) ** 2) + ((df.z * np.cos(df.phi) / df.r) ** 2)
) + np.sin(df.phi) ** 2 * df.psi**4
df["gxy"] = (
    df.psi**4
    * np.exp(2 * df.q)
    * (df.x * df.y / df.r**2 + df.z**2 * np.cos(df.phi) * np.sin(df.phi) / df.r**2)
    - np.sin(df.phi) * np.cos(df.phi) * df.psi**4
)
df["gxz"] = (
    df.psi**4
    * np.exp(2 * df.q)
    * (
        df.x * df.z / df.r**2
        - df.z * np.cos(df.phi) * np.sqrt(df.x**2 + df.y**2) / df.r**2
    )
)
df["gyy"] = (
    df.psi**4
    * np.exp(2 * df.q)
    * (df.y**2 / df.r**2 + df.z**2 * np.sin(df.phi) ** 2 / df.r**2)
    + np.cos(df.phi) ** 2 * df.psi**4
)
df["gyz"] = (
    df.psi**4
    * np.exp(2 * df.q)
    * (
        df.y * df.z / df.r**2
        - df.z * np.sin(df.phi) * np.sqrt(df.x**2 + df.y**2) / df.r**2
    )
)
df["gzz"] = (
    df.psi**4 * np.exp(2 * df.q) * (df.x**2 + df.y**2) / df.r**2
    + df.psi**4 * np.exp(2 * df.q) * df.z**2 / df.r**2
)


df.loc[cond1, "Kxx"] = (
    -2
    * np.sin(df.phi)
    * (df.Hf * df.z * np.cos(df.phi) + df.He * df.x * np.sin(df.theta))
    / (df.psi**2 * df.r**4)
)
df.loc[cond1, "Kxy"] = (
    df.Hf * df.z * np.cos(2 * df.phi)
    + df.He * (df.x * np.cos(df.phi) - df.y * np.sin(df.phi)) * np.sin(df.theta)
) / (df.psi**2 * df.r**4)
df.loc[cond1, "Kxz"] = (
    np.sin(df.phi)
    * (df.Hf * np.sqrt(df.x**2 + df.y**2) - df.He * df.z * np.sin(df.theta))
    / (df.psi**2 * df.r**4)
)
df.loc[cond1, "Kyy"] = (
    2
    * np.cos(df.phi)
    * (df.Hf * df.z * np.sin(df.phi) + df.He * df.y * np.sin(df.theta))
    / (df.psi**2 * df.r**4)
)
df.loc[cond1, "Kyz"] = (
    np.cos(df.phi)
    * (-df.Hf * np.sqrt(df.x**2 + df.y**2) + df.He * df.z * np.sin(df.theta))
    / (df.psi**2 * df.r**4)
)
df.loc[cond1, "Kzz"] = np.zeros(cond1.sum(), dtype=np.float64)

# For the positive z region
df.loc[~cond1 & cond2, "Kxx"] = -(
    (2 * df.alpha * df.He + df.dbetatdr * df.psi**6 * df.r**4)
    * df.x
    * np.sin(df.phi)
    * np.sin(df.theta)
    + df.z
    * np.cos(df.phi)
    * (
        df.dbetatdt * df.psi**6 * df.r**3 * np.sin(df.phi) * np.sin(df.theta)
        + 2 * df.alpha * df.Hf * np.sin(df.phi)
    )
) / (df.alpha * df.psi**2 * df.r**4)
df.loc[~cond1 & cond2, "Kxy"] = (
    (2 * df.alpha * df.He + df.dbetatdr * df.psi**6 * df.r**4)
    * (df.x * np.cos(df.phi) - df.y * np.sin(df.phi))
    * np.sin(df.theta)
    + df.z
    * np.cos(2 * df.phi)
    * (2 * df.alpha * df.Hf + df.dbetatdt * df.psi**6 * df.r**3 * np.sin(df.theta))
) / (2 * df.alpha * df.psi**2 * df.r**4)
df.loc[~cond1 & cond2, "Kxz"] = (
    (
        -2 * df.alpha * df.He * df.z * np.sin(df.phi) * np.sin(df.theta)
        + df.psi**6
        * df.r**3
        * (df.dbetatdt * np.sqrt(df.x**2 + df.y**2) - df.dbetatdr * df.r * df.z)
        * np.sin(df.phi)
        * np.sin(df.theta)
        + 2 * df.alpha * df.Hf * np.sqrt(df.x**2 + df.y**2)
    )
    * np.sin(df.phi)
) / (2 * df.alpha * df.psi**2 * df.r**4)
df.loc[~cond1 & cond2, "Kyy"] = (
    np.cos(df.phi)
    * (
        (2 * df.alpha * df.He + df.dbetatdr * df.psi**6 * df.r**4)
        * df.y
        * np.sin(df.phi)
        * np.sin(df.theta)
        + df.z
        * (
            df.dbetatdt * df.psi**6 * df.r**3 * np.sin(df.phi) * np.sin(df.theta)
            + 2 * df.alpha * df.Hf * np.sin(df.phi)
        )
    )
) / (df.alpha * df.psi**2 * df.r**4)
df.loc[~cond1 & cond2, "Kyz"] = -(
    np.cos(df.phi)
    * (
        -2 * df.alpha * df.He * df.z * np.sin(df.theta)
        + df.psi**6
        * df.r**3
        * (df.dbetatdt * np.sqrt(df.x**2 + df.y**2) - df.dbetatdr * df.r * df.z)
        * np.sin(df.theta)
        + 2 * df.alpha * df.Hf * np.sqrt(df.x**2 + df.y**2)
    )
) / (2 * df.alpha * df.psi**2 * df.r**4)
df.loc[~cond1 & cond2, "Kzz"] = np.zeros((~cond1 & cond2).sum(), dtype=np.float64)

# For the negative z region
df.loc[~cond1 & ~cond2, "Kxx"] = -(
    (2 * df.alpha * df.He + df.dbetatdr * df.psi**6 * df.r**4)
    * df.x
    * np.sin(df.phi)
    * np.sin(df.theta)
    + df.z
    * np.cos(df.phi)
    * (
        -df.dbetatdt * df.psi**6 * df.r**3 * np.sin(df.phi) * np.sin(df.theta)
        + 2 * df.alpha * df.Hf * np.sin(df.phi)
    )
) / (df.alpha * df.psi**2 * df.r**4)
df.loc[~cond1 & ~cond2, "Kxy"] = (
    (2 * df.alpha * df.He + df.dbetatdr * df.psi**6 * df.r**4)
    * (df.x * np.cos(df.phi) - df.y * np.sin(df.phi))
    * np.sin(df.theta)
    + df.z
    * np.cos(2 * df.phi)
    * (2 * df.alpha * df.Hf - df.dbetatdt * df.psi**6 * df.r**3 * np.sin(df.theta))
) / (2 * df.alpha * df.psi**2 * df.r**4)
df.loc[~cond1 & ~cond2, "Kxz"] = -(
    np.sin(df.phi)
    * (
        -2 * df.alpha * df.Hf * np.sqrt(df.x**2 + df.y**2)
        + 2 * df.alpha * df.He * df.z * np.sin(df.theta)
        + df.psi**6
        * df.r**3
        * (
            df.dbetatdt * np.sqrt(df.x**2 + df.y**2) * np.sin(df.theta)
            + df.dbetatdr * df.r * df.z * np.sin(df.theta)
        )
    )
) / (2 * df.alpha * df.psi**2 * df.r**4)
df.loc[~cond1 & ~cond2, "Kyy"] = (
    np.cos(df.phi)
    * (
        (2 * df.alpha * df.He + df.dbetatdr * df.psi**6 * df.r**4)
        * df.y
        * np.sin(df.phi)
        * np.sin(df.theta)
        + df.z
        * (
            -df.dbetatdt * df.psi**6 * df.r**3 * np.sin(df.phi) * np.sin(df.theta)
            + 2 * df.alpha * df.Hf * np.sin(df.phi)
        )
    )
) / (df.alpha * df.psi**2 * df.r**4)
df.loc[~cond1 & ~cond2, "Kyz"] = -(
    np.cos(df.phi)
    * (
        -2 * df.alpha * df.He * df.z * np.sin(df.theta)
        + df.psi**6
        * df.r**3
        * (df.dbetatdt * np.sqrt(df.x**2 + df.y**2) - df.dbetatdr * df.r * df.z)
        * np.sin(df.theta)
        + 2 * df.alpha * df.Hf * np.sqrt(df.x**2 + df.y**2)
    )
) / (2 * df.alpha * df.psi**2 * df.r**4)
df.loc[~cond1 & ~cond2, "Kzz"] = np.zeros((~cond1 & ~cond2).sum(), dtype=np.float64)

s8 = time.time()
if args.v:
    print(f"K and gamma in Cartesian coordinates calculated in {s8-s7} secs")

    # savecols = df[['r', 'theta', 'phi', 'alpha', 'beta__x', 'beta__y', 'beta__z', 'psi', 'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz', 'Kxx', 'Kxy', 'Kxz', 'Kyy', 'Kyz', 'Kzz', 'rho', 'u__x', 'u__y', 'u__z']]
    # ret_df = dfall[['r', 'theta', 'phi', 'alpha', 'beta__x', 'beta__y', 'beta__z', 'psi', 'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz', 'Kxx', 'Kxy', 'Kxz', 'Kyy', 'Kyz', 'Kzz', 'rho', 'u__x', 'u__y', 'u__z', 'A_x', 'A_y', 'A_z']]
    # df.to_hdf(args.folder + '3D_data/databeforevectpotential.h5', key='df', mode='w')

    print("-----------------------------------------------")
    print(f"saving done in {time.time() - s8} secs ")
    # print(dict(df.dtypes))
    print("-----------------------------------------------")
    print("performing vector potential calculations")
s12 = time.time()
# Vector Potential Calculation in lower indices form
# print(df)
# print(df.describe())

df["F_vp"] = np.where(
    (~(df["b2"] < 0)) & (~(df["psi"] == 0)),
    np.sign(df["b_phi"])
    * np.abs(np.sqrt(df["b2"]))
    * df["r"]
    * (np.exp(2 * df["q"]))
    * ((df["psi"]) ** 4),
    np.float64(0),
)
# df['F_vp'] = df['F_vp'].replace([np.nan, np.inf, -np.inf], 0)

pivot = df.pivot_table(
    index="theta", columns="r", values="F_vp", aggfunc="first"
).astype("float64")
pivot.columns = pivot.columns.astype("float64")
A_phi = np.cumsum(np.diff(pivot, prepend=0, axis=1) * np.cumsum(pivot, axis=1), axis=1)
A_phi = np.where(np.isfinite(A_phi), A_phi, np.float64(0))
A_phi[np.isnan(A_phi)] = np.float64(0)  # Note: use '=' instead of '=='

# Convert A_phi back to a DataFrame
A_phi_df = pd.DataFrame(
    A_phi, index=pivot.index, columns=pivot.columns.astype("float64")
)
# print(A_phi_df)
result = A_phi_df.reset_index().melt(id_vars="theta", var_name="r", value_name="A_phi")
result["r"] = result["r"].astype("float64")
# print(result)
# result.to_hdf(args.folder + '3D_data/vectorpotentialresult.h5', key='df', mode='w')
df = pd.merge(df, result, on=["theta", "r"], how="left")
# Theta 0 -> all r values
# ...
# Theta N
df["A_x"] = df["J_20"] * df["A_phi"]
df["A_y"] = df["J_21"] * df["A_phi"]
df["A_z"] = np.float64(0)
df["A_x"] = df["A_x"].replace([np.nan, np.inf, -np.inf], np.float64(0))
df["A_y"] = df["A_y"].replace([np.nan, np.inf, -np.inf], np.float64(0))
# print(df)
# print(df.describe())
if args.v:
    print("-----------------------------------------------")
    # df.to_hdf(args.folder + '3D_data/dfaftervectorpotential.h5', key='df', mode='w')
    print("calculations completed in {} seconds".format(time.time() - s12))

    print("-----------------------------------------------")


# populating below the horizon
s9 = time.time()
df = df.drop_duplicates()
horizondf = df.copy()

horizondf = horizondf[~(np.abs(horizondf.r - rs) < 0.00001)]
horizondf["r"] = horizondf["r"].apply(lambda x: (rs**2) / x)
horizondf["alpha"] = -horizondf["alpha"]
horizondf["u__z"] = -horizondf["u__z"]
horizondf["u__x"] = -horizondf["u__x"]
horizondf["u__y"] = -horizondf["u__y"]

dfall = pd.concat([df, horizondf], axis=0, ignore_index=True)

dfall = dfall.drop_duplicates()
s10 = time.time()
if args.v:
    print(f"data under the horizon has been populated in {s10 - s9} seconds")
# Save the necessary values
# ret_df = dfall[['r', 'theta', 'phi', 'alpha', 'beta__x', 'beta__y', 'beta__z', 'psi', 'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz', 'Kxx', 'Kxy', 'Kxz', 'Kyy', 'Kyz', 'Kzz', 'rho', 'u__x', 'u__y', 'u__z']]

ret_df.to_hdf(args.folder + "3D_data/3d_datadump_routine1.h5", key="df", mode="w")

ret_df = dfall[
    [
        "r",
        "theta",
        "phi",
        "alpha",
        "beta__x",
        "beta__y",
        "beta__z",
        "psi",
        "gxx",
        "gxy",
        "gxz",
        "gyy",
        "gyz",
        "gzz",
        "Kxx",
        "Kxy",
        "Kxz",
        "Kyy",
        "Kyz",
        "Kzz",
        "rho",
        "u__x",
        "u__y",
        "u__z",
        "A_x",
        "A_y",
        "A_z",
    ]
]

origin = {
    "r": [np.float64(0)],
    "theta": [np.float64(0)],
    "phi": [np.float64(0)],
    "alpha": [np.float64(1)],
    "beta__x": [np.float64(0)],
    "beta__y": [np.float64(0)],
    "beta__z": [np.float64(0)],
    "psi": [np.float64(1)],
    "gxx": [np.float64(1)],
    "gxy": [np.float64(0)],
    "gxz": [np.float64(0)],
    "gyy": [np.float64(1)],
    "gyz": [np.float64(0)],
    "gzz": [np.float64(1)],
    "Kxx": [np.float64(0)],
    "Kxy": [np.float64(0)],
    "Kxz": [np.float64(0)],
    "Kyy": [np.float64(0)],
    "Kyz": [np.float64(0)],
    "Kzz": [np.float64(0)],
    "rho": [np.float64(0)],
    "u__x": [np.float64(0)],
    "u__y": [np.float64(0)],
    "u__z": [np.float64(0)],
    "A_x": [np.float64(0)],
    "A_y": [np.float64(0)],
    "A_z": [np.float64(0)],
}

origin_df = pd.DataFrame(origin)
ret_df = pd.concat([ret_df, origin_df], axis=0, ignore_index=True)

ret_df = ret_df.drop_duplicates()
# dfall = pd.concat([df, horizondf], axis=0, ignore_index=True)
ret_df.to_hdf(args.folder + "3D_data/all_data_routine1.h5", key="df", mode="w")


if args.v:
    print(f"saving done in {time.time() - s8} secs ")

    print(
        f"Run took: {time.time() - start_time} seconds 3D_data/all_data_updated_jacobian.h5"
    )
