# BH Disk Initial Data Project
*Jonah Doppelt, Shreyas Jammi*

## Table of Contents
[TOC]

## How to run/use our code
The general pipeline for using this Initial Data generator and reader is fairly simple. First I will provide the general steps, then give a specific example. 

* Initialize input parameters and run sho100.f90 to generate 2D data in x-z plane
    * Specifics about the input parameters can be found in the [user input](#User-Input) section. 
    * Use `./run.sh sho100.90` to produce `sho100` executable
    * Run `./sho100` which create 2D data 
    * sho100 also generates a file called `x.dat` which contains diagnostic information about the system (ADM mass, disk total mass, etc)
* Process 2D variables into desired tensorial quantities with `calc_vars.py` to generate 3D data files
    * This operation specifies how many $\varphi$ values to include in the 3D grid (default is 100)
    * This step reflects the data to $-z$ values, imposes axisymmetry by including $\varphi$ values, and does any coordinate transformations or raising of indices as necessary
    * Produces a single h5 data file
    * Run `python calc_vars.py`
* Perform Interpolation onto new grids using `interpolate_grids.py`
    * Each new grid should be in the form: 
    `[x_min, y_min, z_min, dx, dy, dz, Nx, Ny, Nz, MPI_ID]`
    * Input should be a list of grids of the above form
    * Will produce a single h5 file per new grid containing all variable data at each point
    * Interpolation technique can be altered in `interpolate.py` and is default to be a distance-weighted average of the nearest 8 points in the old grid surrounding the new point

## New Code Full Documentation
### run.sh
This is simply a wrapper function which will create a folders for data organization and call the other scripts. It will first generate the 2D data from sho100, then run calculations using calc_vars. It takes in the sho100 fortran file as an argument as mentioned in the previous section. 

### calc_vars.py
This function is responsible for calculating all of the following necessary variables for the evolution of the system:

```['r', 'theta', 'phi', 'gtt', 'psi', 'gtx', 'gty', 'gtz', 'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz', 'Kxx', 'Kxy', 'Kxz', 'Kyy', 'Kyz', 'Kzz', 'P', 'rho', 'u__x', 'u__y', 'u__z', 'B__x', 'B__y', 'B__z']```

We choose to leave each point in our grid labelled in spherical coordinates because in the interpolation step it will be useful for us to leverage the regular structure of our grid wrt r, theta, and phi. 

##### Remove data points on the horizon
We first remove any data points in our grid whose radii are less than or equal to the radius of the BH horizon ($r_s$). This is because during the calculation of certain quantities we will run into divide by zero errors on the horizon. In particular, the formula for extrinsic curvature has a term which includes $1/(r-r_s)$. Though this poses a numerical issue, this does not result in a physical issue because it is multiplied by a term which is proportional to $(r-r_s)^3$. Nevertheless, it is useful to remove data points on the boundary to avoid numerical difficulties, and a fine enough grid will provide us enough information close to the BH horizon. 

##### Enforce axisymmetry (about z-axis) by generating $\varphi$ values 
It is important to first generate the new $\varphi$ grid points and populate with the 2D data with the tensorial components in spherical coordinates. We could not convert to cartesian coordinates and then attempt to rotate them about the z-axis as the tensor components would vary for different $\varphi$ values. 

##### Define extrinsic curvature and then reflect to -z coordinates
We would also like to populate the lower hemisphere with data, however there is a certain tensorial quantity which will indeed change after reflection. The non-zero component of the extrinsic curvature which has $\theta$ dependence will (of course) be impacted when we reflect about the x-y plane. It takes the following form:
$$
K_{\theta \varphi} = K_{\varphi \theta} = \frac{H_F \sin{\theta}}{\psi^2r} + \frac{1}{2\alpha} \psi^4 r^2 \sin^2{\theta} \partial_{\theta} \beta_T 
$$
When we apply our reflection, we are simply duplicating each grid point (along with its data) and making a new point such that $\theta \rightarrow \pi - \theta$. This will make $K_{\theta \varphi} \rightarrow -K_{\theta \varphi}$ (the sine function and partial derivate will aquire a negative sign). Therefore when we populate our negative hemisphere we must be careful to include this change. 

In this step we also calculated the fluid pressure $P = K \rho^\Gamma$ using the values of K and $\gamma$ as output form the sho100 x.dat file. 

##### Assemble the Jacobian and Metric
Next we assemble the 4x4 jacobian as well as the 4x4 metric (in its covariant and contravariant forms). The jacobian will be used to transform our tensors from spherical to cartesian coordinates. The metric will be used to raise and lower indices on other tensors. Notes on the exact mathematics (including formulas) of both are in the [mathematical review](#Mathematical-Review) section. 

##### Calculate Fluid Velocity and Magnetic Field components
Fluid velocity only has $u^\varphi$ component:
$$
u^\varphi = \Omega*(\alpha^2 - \psi^4r^2\sin^2\theta ( \Omega + \beta)^2)^{-1/2}
$$
We can use the Jacobian to convert from spherical to cartesian. A simple conversion 



### interpolation.py



## Notes on Patryk's Code (sho100.f95)
#### User Input
Begins around line 157
```fortran=157
! Black hole parameters: $m$ and $a$

kerrm = 1.0d0
kerra = 0.8d0*kerrm

! diskparams initial parameters: These are initial guesses for the values of the angular velocity
! $\Omega$ at the inner ($\Omega_1$) and outer $\Omega_2$ equatorial radii of the disk. The actual
! values are obtained later using the Newton-Raphson scheme.

omega1 = 0.03d0
omega2 = 0.004d0
...
```
Here is a list of the input parameters along with the values we're using and their meanings:
* BH Params
    * kerrm = 1.0: Mass of BH
    * kerra = 0.8: Angular Momentum of BH
* Torus Params
    * rr1 = 12: Inner radius of disk
    * rr2 = 30: Outer radius of disk
    * Omega1 = 0.03: Initial guess for angular momentum at inner radius
    * Omega2 = 0.004: Initial guess for angular momentum outer radius
    * rho0 = 8.8e-5: Maximum disk density
* Magnetic Field Params
    * nn = 1.0 (unclear)
    * c2 = 0.01 (unclear)
* Grid Params
    * nr = 800: Number of radii to be used in the grid
    * nt = 200: Number of theta to be used in the grid
* General Params/Constants
    * niter = 1000: Iterations in main loop to achieve convergence
    * restart = 0: Flag whether or not to use restart file
    * gam = 4/3: Polytropic constant
    * arotlaw = 0: A parameter defined in the Keplerian rotation law
    * pi = acos(-1.0)
    * tcut = 0.25*pi
    * adjustjdiskcut = 2




#### Defining the Grid
Begins around line 353
```fortran=353
! Constructing the grid

pi = acos(-1.0d0)

rin = rs
dr  = rin/50.0d0
fac = 1.01d0

r(1) = rin

do i = 2, nrr
   r(i) = r(1) + (fac**(i-1) - 1.0d0)*dr/(fac - 1.0d0)
end do
...
```

##### Radial Grid:
The radial grid is defined by having the inner most radius be equal to the BH radius $r_s = \frac{1}{2}\sqrt{m-a^2}$ where $m,a$ are the BH mass and angular momentum respectively. Then we define $dr = r_s/50$ and $fac = 1.01$. Then the radial array is generated using the formula:
$$
r(i) = r_s + \frac{fac^{i-1}-1}{fac -1}*dr
$$
for $i \in (2,n_r+2)$

This will generate a very fine grid near the BH but will become exponentially more spaced


##### Theta Grid:
The angular grid has two options:
1) Evenly spaced theta values between 0 and $\pi/2$ based on the number of theta values $nt$
2) Differential spacing such that there is a higher density of points closer to the equatorial plane (this is the default)

Option 1 is self explanatory and by default not used. 
Option 2 first defines the value $d \mu = \frac{1}{n_t}$ and then generates a grid of $\mu$ and $\theta$ values using the formula:
$$
\mu(j) = 1 + 0.5d\mu - (j-1)d\mu
$$
$$
\theta(j) = \arccos{\mu(j)}
$$

for $i \in (2,n_t+2)$

#### Variables Calculated/Output from Code
Deallocation of space begins on line 1668
```fortran=1668
! Deallocating storage space

!deallocate(wsp44 )
!deallocate(wsp45 )
...
...
deallocate(phi      )
deallocate(capitalb )
deallocate(betat    )
deallocate(q        )
deallocate(betak    )
deallocate(psi      )
```

This is convenient because in fortran you must declare every variable you plan on using (even intermediary) and then deallocate that space. Through their code, they calculate and store many very useful quantities. All of the available quantities are listed in this section of the code. We can then use these available quantities readily and output them directly without having to recalculate them. 

#### Example: Writing Variables to a File
```fortran=1634
! creating the file
open(114, file="file_name.dat", form="formatted")

! Loop over i and j (dummy indices for r and theta)
do i = nra-1, nrb+1
    do j = nta-1, ntb+1
        ! Write values
        write(114,*) r(i)*sin(t(j)), r(i)*cos(t(j)), r(i), t(j), &
                    rho(i,j), psi(i,j), alpha(i,j), &
                    dbetatdrm(i,j), dbetatdtm(i,j)
    end do
    ! Write a blank line after each row
    write(114,*)
end do

! Close file
close(114)
```
This is an example of how to simply write scalar variables over the grid; in this case I am storing $x,z,r,\theta,\rho,\psi,\alpha, \partial_r \beta_T, \partial_{\theta} \beta_T$ on each row of the file. Note that the order of the do loops determines the output structure; in this case loop over $\theta$'s with fixed r value, then change r. 

## Mathematical Review
#### Coordinate transformations
We begin with the jacobian matrix from spherical to cartesian coordinates in terms of spherical coordinates (including time):

$$ J = \frac{d(t,x,y,z)}{d(t, r, \theta, \phi)} = 
\begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & \sin(\theta) \cos(\varphi) &  r\cos(\theta)\cos(\varphi) & -r\sin(\theta)\sin(\varphi)  \\
0 & \sin(\theta)\sin(\varphi) & r\cos(\theta) \sin(\varphi) & r\sin(\theta)\cos(\varphi) \\
0 & \cos(\theta) & -r\sin(\theta) & 0
\end{pmatrix}
$$
The full 4x4 Jacobian will be useful for converting our 4x4 metric $g_{\mu \nu}$ from spherical to cartesian coordinates. For tensors without the time component (like the extrinsic curvature $K_{\mu \nu}$), only the 3x3 spatial components of the Jacobian are necessary.


For a given tensor quantity ($A_{\mu \nu}$) defined in spherical coordinates, the transformation to cartesian coordinates will be:

$$\tilde{A}_{\mu' \nu'} = \begin{pmatrix}
A_{xx} & A_{xy} & ... \\
A_{yx} & ... &  \\
... &  & A_{zz}
\end{pmatrix} = J^{T} A_{\mu \nu} J
$$

#### Raising and Lowering Indices for Tensors (Metric Information)
Tensors can either be in their covariant (lower index) or contravariant (upper index) form. If we want to change which form they are in, we use the metric tensor $g_{\mu \nu} = (g^{\mu \nu})^{-1}$ to raise or lower the indices. In our BH-disk system, our metric takes the following forms in spherical coordinates:
In 3+1:
$$
g = (-\alpha^2 + \beta_i \beta^i)dt^2 + 2\beta_i dt dx^i + \gamma_{ij}dx^i dx^j
$$
with nonzero components being
\begin{align}
(\beta^r, \beta^{\theta}, \beta^{\varphi}) &= (0, 0, \beta) \\
(\beta_r, \beta_{\theta}, \beta_{\varphi}) &= (0, 0, \psi^4 r^2 \beta \sin^2{\theta} ) \\
(\gamma_{rr}, \gamma_{\theta \theta}, \gamma_{\varphi \varphi}) &= (\psi^4 e^{2q}, \psi^4 e^{2q}r^2, \psi^4 r^2 \sin^2{\theta})
\end{align}
In matrix representation this gives:
\begin{equation}
g_{\mu \nu} = 
\begin{pmatrix}
-\alpha^2 + \psi^4 r^2 \beta^2 \sin^2{\theta} & 0 & 0 & \psi^4 r^2 \beta \sin^2{\theta}\\
0 & \psi^4 e^{2q} & 0 & 0 \\
0 & 0 & \psi^4 e^{2q}r^2 & 0 \\
\psi^4 r^2 \beta \sin^2{\theta} & 0 & 0 & \psi^4 r^2 \sin^2{\theta}
\end{pmatrix}
\end{equation}
\begin{equation}
g^{\mu \nu} =
\begin{pmatrix}
1/\alpha^2 & 0 & 0 & \beta/\alpha^2\\
0 & 1/(\psi^4 e^{2q}) & 0 & 0 \\
0 & 0 & 1/(\psi^4 e^{2q}r^2) & 0 \\
\beta/\alpha^2 & 0 & 0 & 1/(\psi^4 r^2 \sin^2{\theta}) - \beta^2/\alpha^2
\end{pmatrix}
\end{equation}

We would like to transform any tensorial value into its contravariant form before outputting it into out ID. If we have a tensor in its covariant form $A_{\mu \nu}$, we can raise its indices by performing:
\begin{equation}
A^{\tau \lambda} = g^{\tau\mu} A_{\mu \nu} g^{\lambda \nu}
\end{equation}

Note, there are times when we will both want to transform coordinates AND raise indices. We can do this by either changing our coordinates on both $A_{\mu \nu}$ and our metric, and then raising indices. Or raising indices in spherical coordinates and then performing a coordinate transformation. Here is a proof that these operations are equivalent (latin indices indicate cartesian coordinates while greek indices indicate spherical):
\begin{align}
A^{ij} &= g^{ik} A_{kl} g^{lj} \\
&= (J^Tg^{\tau\mu}J) (J^T A_{\mu \nu} J)(J^Tg^{\lambda \nu}J) \\
&= (J^Tg^{\tau\mu})(J J^T) A_{\mu \nu} (J J^T)(g^{\lambda \nu}J) \\
&= (J^Tg^{\tau\mu})A_{\mu \nu}(g^{\lambda \nu}J) \\
&= J^T(g^{\tau\mu}A_{\mu \nu}g^{\lambda \nu})J \\
&= J^T A^{\tau \lambda} J
\end{align}
Where in steps 3-4 we leveraged the fact that the Jacobian is a real unitary matrix and thus $J^T = J^{-1}$.

#### Application to Extrinsic Curvature $K_{\mu \nu}$
The only non-vanishing terms of $K_{\mu \nu}$ in spherical coordinates is are:
$$
K_{r \varphi} = K_{\varphi r} = \frac{H_E \sin^2{\theta}}{\psi^2r^2} + \frac{1}{2\alpha} \psi^4 r^2 \sin^2{\theta} \partial_r \beta_T 
$$
$$
K_{\theta \varphi} = K_{\varphi \theta} = \frac{H_F \sin{\theta}}{\psi^2r} + \frac{1}{2\alpha} \psi^4 r^2 \sin^2{\theta} \partial_{\theta} \beta_T 
$$

Where $H_E$ and $H_F$ are defined in the paper and calculated directly in Patryk's code for us. In fact, all variables needed to calcuate these components of the extrinsic curvature are provided calculated directly in the code.

We can then initialize our $K_{\mu \nu}$ tensor in spherical coordinates at every point in our 2D x-z plane as: 

$$K_{\mu \nu} = \begin{pmatrix}
0 & 0 & K_{r \varphi} \\
0 & 0 & K_{\theta \varphi} \\
K_{\varphi r} & K_{\varphi \theta } & 0
\end{pmatrix}
$$

We make sure to populate the $-z$ values by copying every value in our original grid but asserting that $\theta_{new} = \pi - \theta$. We can then populate a 3-D grid with axisymmetry by defining a parameter $n_{\varphi}$ which is the number of $\varphi$ values and then copying every value in our 2D grid and assigning it a serie of $\varphi$ values. This works because the tensors don't depend on $\varphi$ when defined in spherical coordinates. 

Now we perform the coordinate transformation using the Jacobian as outlined previously:

$$\tilde{K}_{\mu' \nu'} = \begin{pmatrix}
K_{xx} & K_{xy} & ... \\
K_{yx} & ... &  \\
... &  & K_{zz}
\end{pmatrix} = J^{T} K_{\mu \nu} J
$$

Our python script uses numpy's matrix multiplying function to perform this operation on each point in our now fully populated 3D grid and then writes the 6 unique cartesian components to a final file. (I am not sure how to do this operation in fortran but I'm sure it would be more efficient). 

## Cocal Initial Data (ID) code structure
Here is an example of the cocal initial data output structure for the extrinsic curvature
```fortran=1
subroutine IO_output_Kij_3D_WL
    use phys_constant, only : long
    ... 
    !imports and variable declarations
    ...
  !
  ! --- Metric potentials.
    oo3 = 1.0d0/3.0d0
  
    open(13,file='rnsgra_Kij_3D.las',status='unknown')
    write(13,'(5i5)') nrg, ntg, npg
    do ipg = 0, npg
      do itg = 0, ntg
        do irg = 0, nrg
          axx = tfkij_grid(irg,itg,ipg,1,1)
          axy = tfkij_grid(irg,itg,ipg,1,2)
          ... 
          !mathematical operations
          ...
          kyz = ayz + oo3*gyz*tk        
          kzz = azz + oo3*gzz*tk        
  
          write(13,'(1p,6e23.15)') kxx,kxy,kxz,kyy,kyz,kzz
        end do
      end do
    end do
    close(13)
  !
  end subroutine IO_output_Kij_3D_WL
```

The tfkij_grid has some structure such that each point is indexed monotonically. The order of output variables is fixed $\varphi$, $\theta$ and loop through indices in $r$. Then increase $\theta$, then increase $\varphi$. 

## 2D to 3D initial data volume generator

the initial data as created by the sho100 initial data is generated across the positive x-z axis, so the following describes the structure of how we convert this two dimensional data into three dimensional cartesian data that can be fed into the IllinoisGRMHD evolution codebase. 

From the previous section on Patryck's code notes, we understand that the 2D initial data is outputed in the following format:

```dat=
x_1 z_1 [data]
x_2 z_2 [data]
x_3 z_3 [data]
x_4 z_4 [data]
...
```

Knowing this, we need to write a python script to read the dat file. This function will read the dat file and load the numpy data array, aptly named **data_loader** because it loads the data. 

```python=
import numpy as np

def data_loader(filename: str):
    '''
    opens the file and loads the data into a numpy array

    Args:
        filename (str) : the name of the file to be opened
    Returns:
        data (np.array) : a numpy array with the data from the file
    '''
    data = []
    with open(filename, "r") as f:
        for line in f:
            line_tmp = line.strip().split(' ')
            data.append([float(x) for x in line_tmp])
    
    return np.array(data)
```




Working with numpy arrays is prefered because it is very effective with easily minipulating large datasets. From here, we want to flip the data across the z-axis. This can be done by generating a copy of the data, negating the z-axis coordinate, and adding the copy of the data to the original data, which can be seen in the following function:

```python=
import numpy as np

def grid_flipper(data: np.ndarray):
    '''
    flips the data about the z-axis

    Args:
        data (np.array) : the data to be flipped
    Returns:
        newdata (np.array) : the flipped data
    '''
    tmp_data = data.copy()
    tmp_data[:, 1] = -tmp_data[:, 1]
    newdata = np.vstack((tmp_data, data))
    return newdata
```

From here, we want to rotate the x-z data across the y-axis. In order to conceptualize this, it might be easier to correlate this to cylindrical data with varying magnitudes of equidistant phi values. Generally, we want to create an even distribution of phi values to cover the entire grid with proper resolution. The trade-off that we find here is that the more resolution that we have (which would be phi values in this case), the more computationally intensive the process is for generating data. To start, we will create 40 equidistant phi values between 0 to $2 \pi$. From here, we will take each point from the original data and rotate the data by the aforementioned phi value until we have completed all of the phi values and mapped all of the data properly. The following function is how this might be done:

```python=
import numpy as np

def cart_to_spherical(data: np.ndarray, phivals: int):
    '''
    converts the x-z data to x-y-z data by rotating it across the y-axis

    Args:
        data (np.array) : the data to be converted
        phivals (int) : the angle to rotate the data by
    Returns:
        converteddata (np.array) : the converted data
    '''
    phivalues = np.radians(np.arange(0, 360, phivals))
    converteddata = np.empty((0, data.shape[1]+1))
    #converteddata = np.empty((720, 1+array.shape[0]))
    for phi in phivalues:
        tmp = data.copy()
        xval  = np.cos(phi)*tmp[:, 0]
        yval = np.sin(phi)*tmp[:, 0]
        tmp = tmp[:, 1:]
        tmp = np.insert(tmp, 0, xval, axis=1)
        tmp = np.insert(tmp, 1, yval, axis=1)
        converteddata = np.vstack((converteddata, tmp))

    return converteddata
```



## Summary of our Initial Data Generator
There are areas for considerable optimization in this process, but here is the current pipeline of our data generation:

1. Input initial parameters into the beginning of sho100
2. Determine which output variables are necessary for calculating the desired grid quantities and write them to a file 
3. Use a python script on the output variables in order to do simple algebraic and matrix operations on each grid point, and output a 3D grid containing the desired scalar variables. 

Areas for optimization:
1. Python is likely slower than fortran to do these operations, so we should be able to do all the algebraic and matrix operations in the initial sho100 file and output the direct values we want without having to use a python script. 


