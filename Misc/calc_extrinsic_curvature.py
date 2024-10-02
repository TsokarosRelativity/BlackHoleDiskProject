import numpy as np
import time
start_time = time.time()

output_dir = "output_files/" ##end with /
output_file = 'test3D.ini'

# Open the input file
with open('extrisic_curvature_vars.dat', 'r') as f:
    lines = f.readlines()

# Create a new list to store the lines for the new file
new_lines = []

n_phi = 40

def calc_J(r, theta, phi):
    # Initialize the Jacobian matrix
    J = np.zeros((3, 3))

    # Calculate the Jacobian matrix elements
    J[0, 0] = np.cos(theta)*np.sin(phi)
    J[0, 1] = -r * np.sin(theta) * np.sin(phi)
    J[0, 2] = r * np.cos(theta) * np.cos(phi)

    J[1, 0] = np.sin(theta) * np.sin(phi)
    J[1, 1] = r * np.cos(theta) * np.sin(phi)
    J[1, 2] = r * np.sin(theta) * np.cos(phi)

    J[2, 0] = np.cos(theta)
    J[2, 1] = 0.0
    J[2, 2] = -r*np.sin(phi)

    return J
# Iterate through each line
for line in lines:
    # Split the line into columns
    columns = line.split()
    if len(columns) == 0:
        continue

    x,z,r,t,He,Hf,psi,alpha,dbetadr,dbetadt = [float(x) for x in columns]
    # rho = float(columns[2])
    # h = float(columns[3])
    if alpha == 0:
        continue
    coef = (1/(2*alpha)) * psi**4 * r**2 * np.sin(t)**2
    Kpsir = (He*np.sin(t)**2)/(psi**2 * r**2) + coef*dbetadr
    Kpsit = (Hf*np.sin(t))/(psi**2 * r) + coef*dbetadt

    #Assemble extrisic curvature K in cartesian coords
    K = np.zeros((3, 3))
    K[0,2] = Kpsir
    K[2,0] = Kpsir
    K[1,2] = Kpsit
    K[2,1] = Kpsit

    for i in range(n_phi):
        phi = 2*np.pi*(i/n_phi)

        x1 = r*np.cos(t)*np.cos(phi)
        y1 = r*np.sin(t)*np.sin(phi)
        z1 = r*np.cos(t)
        
        # Assemble Jacobian Matrix
        J = calc_J(r,t,phi)

        # Convert to cartesian
        K_cart = np.matmul(np.matmul(J.T,K),J).ravel() #.ravel() flattens to a 1D array

        Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz = K_cart #unpack
        

        new_line = f"{x1} {y1} {z1} {Kxx} {Kxy} {Kxz} {Kyx} {Kyy} {Kyz} {Kzx} {Kzy} {Kzz}\n"
        new_lines.append(new_line)


with open(output_dir + output_file, 'w') as f:
    f.writelines(new_lines)

print(f"New file '{output_file}' has been created.")
print(f"Run took: {time.time() - start_time} seconds")

