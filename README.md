# BH Disk Initial Data Project
Evolution of a Black Hole Disk with a Toroidally Magnetized Field Lines
*Jonah Doppelt, Shreyas Jammi*

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

$$K_{\theta \varphi} = K_{\varphi \theta} = \frac{H_F \sin{\theta}}{\psi^2r} + \frac{1}{2\alpha} \psi^4 r^2 \sin^2{\theta} \partial_{\theta} \beta_T $$

When we apply our reflection, we are simply duplicating each grid point (along with its data) and making a new point such that $\theta \rightarrow \pi - \theta$. This will make $K_{\theta \varphi} \rightarrow -K_{\theta \varphi}$ (the sine function and partial derivate will aquire a negative sign). Therefore when we populate our negative hemisphere we must be careful to include this change. 

In this step we also calculated the fluid pressure $P = K \rho^\Gamma$ using the values of K and $\gamma$ as output form the sho100 x.dat file. 

##### Assemble the Jacobian and Metric
Next we assemble the 4x4 jacobian as well as the 4x4 metric (in its covariant and contravariant forms). The jacobian will be used to transform our tensors from spherical to cartesian coordinates. The metric will be used to raise and lower indices on other tensors. Notes on the exact mathematics (including formulas) of both are in the mathematical review section. 

##### Calculate Fluid Velocity and Magnetic Field components
Fluid velocity only has $u^\varphi$ component:

$$u^\varphi = \Omega*(\alpha^2 - \psi^4r^2\sin^2\theta ( \Omega + \beta)^2)^{-1/2}$$

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

$$r(i) = r_s + \frac{fac^{i-1}-1}{fac -1}*dr$$

for $i \in (2,n_r+2)$

This will generate a very fine grid near the BH but will become exponentially more spaced


##### Theta Grid:
The angular grid has two options:
1) Evenly spaced theta values between 0 and $\pi/2$ based on the number of theta values $nt$
2) Differential spacing such that there is a higher density of points closer to the equatorial plane (this is the default)

Option 1 is self explanatory and by default not used. 
Option 2 first defines the value $d \mu = \frac{1}{n_t}$ and then generates a grid of $\mu$ and $\theta$ values using the formula:
$$\mu(j) = 1 + 0.5d\mu - (j-1)d\mu$$
$$\theta(j) = \arccos{\mu(j)}$$

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
