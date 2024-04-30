import numpy as np
import time
import pandas as pd
start_time = time.time()

######### INPUTS
file_path = '2D_data/all_vars.dat'
rs = np.float64(0.3)
n_phi = 100

K = 5.3206252168330817E-002 #this is K from the polytropic equation P = K*(rho**Gamma)... K can be found in x.dat as k
gamma = 4/3 #can also be found in x.dat

s1 = time.time()
column_names = ['x','z','r', 'theta', 'alpha', 'psi', 'q', 'betak', 'betat', 'dbetatdr', 'dbetatdt', 'He', 'Hf', 'Omega', 'b2', 'rho']
df = pd.read_csv(file_path, delimiter='\s+', header=None, names=column_names)
s2 = time.time()
print(f'read in the data in {s2-s1} secs')

# Exclude rows where the values in column 'r' are within the tolerance of rs
tolerance = 0.00001 
df = df[~(np.abs(df['r'] - rs) < tolerance)]

s3 = time.time()

# Calculate evenly spaced phi values
df = pd.DataFrame(np.repeat(df.values, n_phi, axis=0), columns=df.columns)
num_rows = len(df)
phi_values = np.linspace(0, 2*np.pi, num_rows % n_phi or n_phi)

# Add phi values as a new column
df['phi'] = np.tile(phi_values, num_rows // n_phi + 1)[:num_rows]
s4 = time.time()
print(f'added phi values in {s4-s3} secs')

df['beta'] = df.betat + df.betak
df['K_rp'] = (df.He*np.sin(df.theta)**2) / (df.psi**2 * df.r**2) +  (2*df.alpha)**(-1) * df.psi**4 * df.r**2 * np.sin(df.theta)**2 * df.dbetatdr
df['K_tp'] = (df.Hf*np.sin(df.theta)) / (df.psi**2 * df.r) + (2*df.alpha)**(-1) * df.psi**4 * df.r**2 * np.sin(df.theta)**2 * df.dbetatdt
df['P'] = K*(df.rho**gamma)

# Add theta values for -z values
df_temp = df.copy()
df_temp['theta'] = np.pi - df_temp['theta']
df_temp['K_tp'] = -df_temp['K_tp']
df = pd.concat([df,df_temp]).reset_index(drop=True)
s5= time.time()
print(f'reflected to -z in {s5-s4} secs')

# Assemble Jacobian (4x4)

df['J_00'] = 1
df['J_11'] = np.sin(df.theta)*np.cos(df.phi)
df['J_12'] = df.r * np.cos(df.theta) * np.cos(df.phi)
df['J_13'] = -df.r * np.sin(df.theta) * np.sin(df.phi)

df['J_21'] = np.sin(df.theta) * np.sin(df.phi)
df['J_22'] = df.r * np.cos(df.theta) * np.sin(df.phi)
df['J_23'] = df.r * np.sin(df.theta) * np.cos(df.phi)

df['J_31'] = np.cos(df.theta)
df['J_32'] = -df.r*np.sin(df.theta)
df['J_33'] = 0

df['zero'] = 0

# Assemble metric in spherical coords (covariant form)
df['g_11'] = df.psi**4 * np.exp(2*df.q)                                 #g_rr
df['g_22'] = df.g_11*df.r**2                                            #g_thetatheta
df['g_33'] = df.psi**4 * df.r**2 * np.sin(df.theta)**2                  #g_phiphi
df['g_03'] = df.g_33 * df.beta                                          #shift
df['g_00'] = -df.alpha**2 + df.g_33 * df.beta**2                        #lapse

# Assemble metric in spherical coords (contravariant form)
df['g__00'] = -1/(df.alpha**2)
df['g__03'] = df.beta/(df.alpha**2)
df['g__11'] = 1/df.g_11
df['g__22'] = 1/df.g_22
df['g__33'] = 1/df.g_33 - (df.beta/df.alpha)**2


# Calc Fluid Velocity (These are upper indices indicated by two lower '__')
df['u__t'] = (df.alpha**2 - df.g_33 * (df.Omega + df.beta)**2)**(-1/2) 
df['u__phi'] = df.Omega * df.u__t
df['u__x'] = df.J_13 * df.u__phi
df['u__y'] = df.J_23 * df.u__phi
df['u__z'] = 0

# Calc Mag Field Components (These are lower indices)
df['b2_phi'] = df.u__t**2 * df.g_33 * df.alpha**2 * df.b2
df['b__phi'] = df.g__33 * df.b2_phi**(1/2)   #contravariant form of the phi component of b (the only nonzero component)
df['b_t'] = -df.Omega * df.b2_phi**(1/2)

G = 6.674 * 10**(-11)           #gravitational constant in SI units
c = 3 * 10**(8)                 #speed of light in SI units
epsilon_0 = 8.854 * 10**(-12)   #permitivity of free space in SI units
unit_conv = G**(1/2) * c**(-1) * epsilon_0**(-1/2) * 10**(4) #conversion from geometrized units for magnetic B field to Gauss (CGS units)

df['B__x'] = (4*np.pi)**(1/2) * df.J_13 * df.b__phi * unit_conv #contravariant mag field in CGI (gauss) units
df['B__y'] = (4*np.pi)**(1/2) * df.J_23 * df.b__phi * unit_conv
df['B__z'] = 0


s6 = time.time()
print(f'calculations done in {s6-s5} secs')

#### Perform matrix multiplication to do a coordinate transformation from spherical to cartesian (using contravariant form)
g4 = df[['g__00', 'zero', 'zero', 'g__03', 'zero', 'g__11', 'zero', 'zero', 'zero', 'zero', 'g__22', 'zero', 'g__03', 'zero', 'zero', 'g__33']].values.reshape(-1, 4, 4)
J4 = df[['J_00', 'zero', 'zero', 'zero', 'zero', 'J_11', 'J_12', 'J_13', 'zero', 'J_21', 'J_22', 'J_23', 'zero', 'J_31', 'J_32', 'J_33']].values.reshape(-1, 4, 4)
J3 = J4[:,1:,1:]
K = df[['zero', 'zero', 'K_rp', 'zero', 'zero', 'K_tp', 'K_rp', 'K_tp', 'zero']].values.reshape(-1, 3, 3)

g_cart = np.matmul(np.transpose(J4, axes=(0,2,1)), np.matmul(g4, J4))
K_cart_covar = np.matmul(np.transpose(J3, axes=(0,2,1)), np.matmul(K, J3))

## convert from covariant to contravariant
g_cart3 = g_cart[:,1:, 1:]
K_cart = np.matmul(g_cart3, np.matmul(K_cart_covar,g_cart3))


g_cart = g_cart.reshape(-1, 16)
indices = [0, 1, 2, 3, 5, 6, 7, 10, 11, 15]
df[['gtt', 'gtx', 'gty', 'gtz', 'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']] = g_cart[:,indices]

K_cart = K_cart.reshape(-1,9)
indices = [0, 1, 2, 4, 5, 8]
df[['Kxx', 'Kxy', 'Kxz', 'Kyy', 'Kyz', 'Kzz']] = K_cart[:, indices]

df['x'] = df.r * np.sin(df.theta) * np.cos(df.phi)
df['y'] = df.r * np.sin(df.theta) * np.sin(df.phi)
df['z'] = df.r * np.cos(df.theta)

s7 = time.time()
print(f'matrix multiplication done in {s7-s6} secs')

#### Save only the necessary values
# ret_df = df[['r', 'theta', 'phi', 'alpha', 'psi', 'gtx', 'gty', 'gtz', 'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz', 'Kxx', 'Kxy', 'Kxz', 'Kyy', 'Kyz', 'Kzz', 'rho', 'u__x', 'u__y', 'u__z', 'B__x', 'B__y', 'B__z']]
# ret_df.to_hdf('3D_data/all_data.h5', key='df', mode='w')

rho_df = df[['x', 'y', 'z','rho']]
rho_df.to_hdf('3D_data/rho_data.h5', key='df', mode='w')
print(f"Run took: {time.time() - start_time} seconds 3D_data/all_data.h5")
