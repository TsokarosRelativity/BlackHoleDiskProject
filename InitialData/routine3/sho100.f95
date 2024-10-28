program torusy

!$ use omp_lib

implicit none

interface
subroutine solvebetak(betak, integrand, r)
  implicit none
  double precision, dimension(:,:) :: betak, integrand
  double precision, dimension(:)   :: r
end subroutine solvebetak

subroutine integration2d (integrand, outcome, r, t, storage)
implicit none
double precision, dimension(:,:) :: integrand
double precision :: outcome
double precision, dimension(:)   :: r, t, storage
end subroutine integration2d

subroutine integration1d (integrand, outcome, t)
implicit none
double precision, dimension(:) :: integrand
double precision :: outcome
double precision, dimension(:) :: t
end subroutine integration1d

subroutine exact(q, sq, a, b, c, r, t, boundary, kmax, lmax, wsp, IPIV, istart, &
   equ, axis, equator, horizon, betahor)
implicit none

double precision, dimension(:,:) :: q, sq, wsp
double precision, dimension(:,:) :: a, b, c
double precision, dimension(:)   :: r, t, boundary
integer :: kmax, lmax
integer, dimension (:) :: IPIV
logical :: istart
character :: equ
double precision, dimension(:)   :: axis, equator, horizon
double precision, dimension(:,:) :: betahor
end subroutine exact

subroutine derivr(f, dfdr, r)
implicit none
double precision, dimension(:,:) :: f, dfdr
double precision, dimension(:)   :: r
end subroutine derivr

subroutine derivt(f, dfdt, t)
implicit none
double precision, dimension(:,:) :: f, dfdt
double precision, dimension(:)   :: t
end subroutine derivt


subroutine diskparams(c1, w, r1, r2, omega1, omega2, beta1, beta2, alpha1, alpha2, psi1, psi2, a)
implicit none
double precision :: c1, w, r1, r2, omega1, omega2, beta1, beta2, alpha1, alpha2, psi1, psi2, a
end subroutine diskparams

subroutine diskomega(omega, r, t, alpha, beta, psi, w, a)
implicit none
double precision :: omega, r, t, alpha, beta, psi, w, a
end subroutine diskomega

subroutine magbernoulli(h, alpha, psi, r, t, gam, kk, nn, c2)
implicit none
double precision :: h, alpha, psi, r, t, gam, kk, nn, c2
end subroutine magbernoulli

subroutine derivr1(f, dfdr, r)
implicit none
double precision, dimension(:) :: f, dfdr
double precision, dimension(:) :: r
end subroutine derivr1

end interface


double precision, dimension(:,:), allocatable :: phi, capitalb, &
   betat, q, betak, psi, sphi, scapitalb, sbetat, sq, alpha, &
   capitalhe, capitalhf, b, capitala2, storage, integrand, stor1, stor2, stor3, stor4
double precision, dimension(:,:), allocatable :: a44, b44, c44, &
   a45, b45, c45, a46, b46, c46, a47, b47, c47
double precision, dimension(:,:), allocatable :: h, rho, rhocapitalh, &
   capitalj, capitalp
double precision, dimension(:,:), allocatable :: b2
double precision :: rr1, rr2, r1, r2, r3, alfa1, alfa2, psi1, psi2, jj, cc, &
   beta1, beta2, kk, gam, rho0, h0
integer :: i1, i2, i3
double precision, dimension(:),   allocatable :: r, t, mu, storage2, storage3
double precision, dimension(:,:),   allocatable :: boundary44, boundary45, &
  boundary46, boundary47

double precision, dimension(:), allocatable :: axis44, equator44, horizon44
double precision, dimension(:), allocatable :: axis45, equator45, horizon45
double precision, dimension(:), allocatable :: axis46, equator46, horizon46
double precision, dimension(:), allocatable :: axis47, equator47, horizon47

double precision :: kerra, kerrm, rk, sigmak, deltak
double precision :: rs, pi, rin, rout, fac
double precision :: dr1, dr2, dr3, dt1, dt2, dt3, dmu
double precision :: dbdr, dbdt, dphidr, dphidt
double precision :: dbetatdr, dbetatdt, temp
double precision :: capitalm1, capitalmast2, capitaljast, capitalb1, capitalj1, q1
integer :: nr, nt, i, j, iter, niter, restart
double precision :: capitalmast
double precision :: ah, kappa, mirr, mbh, omegah, mh, jh, madm, mt, mt2, mh2
double precision :: mc, ce, cp

double precision :: dqdr, dqdr2, dqdt, dqdt2
double precision, dimension(:,:), allocatable :: test, test2, test3
double precision :: rr, tt, dr, dt
integer :: ntt, nrr, nta, ntb, nra, nrb

double precision, dimension(:,:), allocatable :: ab44, ab45, ab46, ab47

double precision, dimension(:,:), allocatable :: dphidrm, dphidtm, dbdrm, dbdtm, dbetatdrm, dbetatdtm
double precision, dimension(:,:), allocatable :: betahor
double precision, dimension(:,:), allocatable :: capitalomega, uphidn, utup

integer, dimension(:), allocatable :: IPIV44, IPIV45, IPIV46, IPIV47
integer :: jmax, lmax
logical :: istart
double precision :: eps0, prop
double precision :: d1, d2, d3, d4, d5
double precision :: jj2, cc2
character(len=30)     :: outfil, outnum
double precision :: c1, w, omega1, omega2, omega
integer :: j3, ite, jtmax
double precision :: c0, temp1, temp2
double precision :: temp3, temp4
double precision :: tmax, tcut
double precision :: arotlaw
double precision :: nn, c2
double precision :: betamag
integer, dimension(:),   allocatable :: jdiskcut, jdisk
integer :: adjustjdiskcut, jcut, n

double precision, dimension(:), allocatable :: gttf, gtphi, gphiphi, drgttf, ddrgttf, drgtphi, ddrgtphi, &
        drgphiphi, ddrgphiphi, warisco
integer :: iisco
double precision :: risco, riscoc

write(*,*)
write(*,*) '------------------------------------'
write(*,*) '-----------     START    -----------'
write(*,*)

!$omp parallel
write(*,*) 'Running on ', omp_get_num_threads(), ' threads.'
!$omp end parallel

open (unit=165,file="masa.dat",form="formatted",access="append")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   INPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Black hole parameters: $m$ and $a$

kerrm = 1.0d0
kerra = 0.8d0*kerrm

! diskparams initial parameters: These are initial guesses for the values of the angular velocity
! $\Omega$ at the inner ($\Omega_1$) and outer $\Omega_2$ equatorial radii of the disk. The actual
! values are obtained later using the Newton-Raphson scheme.

omega1 = 0.03d0
omega2 = 0.004d0

! Polytropic constant

gam = 4.0d0/3.0d0

! Parameters in the rotation law. This is the black hole spin parameter that appears in the rotation law. 

arotlaw = 0.0d0*kerrm

! Grid resolution:

nr = 800
nt = 200


! Inner and outer radii of the disk
rr1 = 12.0d0
rr2 = 30.0d0

! Density at some particular point of the disk:
!rho0 = 9.0d-5 #Mstar = 0.43
rho0 = 1.2d-4

! Start with the restart files or not
restart = 0

! The number of iteration in the main loop
niter = 10000

! initial cuttoff parameter tcut
pi = acos(-1.0d0)
tcut = 0.25d0*pi

adjustjdiskcut = 2

! magnetic field distribution parameters

nn = 1.0d0
c2 = 0.01d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                           THE END OF INPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Horizon radius

rs = sqrt(kerrm**2 - kerra**2)/2.0d0


! Adding one grid point for each boundary

nrr = nr + 2
ntt = nt + 2
nra = 2
nrb = nr + 1
nta = 2
ntb = nt + 1

jmax = (nrr)*(ntt)
lmax = ntt


! Allocating storage space

allocate(ab44 (7*lmax+1,jmax))
allocate(ab45 (7*lmax+1,jmax))
allocate(ab46 (7*lmax+1,jmax))
allocate(ab47 (7*lmax+1,jmax))

allocate(IPIV44(jmax))
allocate(IPIV45(jmax))
allocate(IPIV46(jmax))
allocate(IPIV47(jmax))

write(*,*) 'size of wsp matrices', 4.0d0*sizeof(ab44)/1.0d9, ' GB'

allocate(axis44      (nrr))
allocate(equator44   (nrr))
allocate(horizon44   (ntt))

allocate(axis45      (nrr))
allocate(equator45   (nrr))
allocate(horizon45   (ntt))

allocate(axis46      (nrr))
allocate(equator46   (nrr))
allocate(horizon46   (ntt))

allocate(axis47      (nrr))
allocate(equator47   (nrr))
allocate(horizon47   (ntt))

allocate(betahor     (ntt,3))

allocate(phi       (nrr,ntt))
allocate(capitalb  (nrr,ntt))
allocate(betat     (nrr,ntt))
allocate(q         (nrr,ntt))
allocate(betak     (nrr,ntt))
allocate(psi       (nrr,ntt))
allocate(sphi      (nrr,ntt))
allocate(scapitalb (nrr,ntt))
allocate(sbetat    (nrr,ntt))
allocate(sq        (nrr,ntt))
allocate(alpha     (nrr,ntt))
allocate(capitalhe (nrr,ntt))
allocate(capitalhf (nrr,ntt))
allocate(b         (nrr,ntt))
allocate(capitala2 (nrr,ntt))


allocate(dphidrm   (nrr,ntt))
allocate(dphidtm   (nrr,ntt))
allocate(dbdrm     (nrr,ntt))
allocate(dbdtm     (nrr,ntt))
allocate(dbetatdrm (nrr,ntt))
allocate(dbetatdtm (nrr,ntt))


allocate(a44       (nrr,ntt))
allocate(b44       (nrr,ntt))
allocate(c44       (nrr,ntt))
allocate(a45       (nrr,ntt))
allocate(b45       (nrr,ntt))
allocate(c45       (nrr,ntt))
allocate(a46       (nrr,ntt))
allocate(b46       (nrr,ntt))
allocate(c46       (nrr,ntt))
allocate(a47       (nrr,ntt))
allocate(b47       (nrr,ntt))
allocate(c47       (nrr,ntt))

allocate(h          (nrr,ntt))
allocate(rho        (nrr,ntt))
allocate(rhocapitalh(nrr,ntt))
allocate(capitalj   (nrr,ntt))
allocate(capitalp   (nrr,ntt))

allocate(storage   (nrr,ntt))
allocate(stor1     (nrr,ntt))
allocate(stor2     (nrr,ntt))
allocate(stor3     (nrr,ntt))
allocate(stor4     (nrr,ntt))
allocate(integrand (nrr,ntt))

allocate(test      (nrr,ntt))
allocate(test2     (nrr,ntt))
allocate(test3     (nrr,ntt))

allocate(r         (nrr))
allocate(storage2  (nrr))
allocate(t         (ntt))
allocate(mu        (ntt))
allocate(storage3  (ntt))

allocate(boundary44 (ntt,1))
allocate(boundary45 (ntt,1))
allocate(boundary46 (ntt,1))
allocate(boundary47 (ntt,1))

allocate(capitalomega (nrr,ntt))
allocate(uphidn       (nrr,ntt))
allocate(utup         (nrr,ntt))

allocate(b2           (nrr,ntt))

allocate(jdisk (nrr))
allocate(jdiskcut (nrr))

allocate(gttf       (nrr))
allocate(gtphi      (nrr))
allocate(gphiphi    (nrr))
allocate(drgttf     (nrr))
allocate(ddrgttf    (nrr))
allocate(drgtphi    (nrr))
allocate(ddrgtphi   (nrr))
allocate(drgphiphi  (nrr))
allocate(ddrgphiphi (nrr))
allocate(warisco    (nrr))


! Constructing the grid

pi = acos(-1.0d0)

rin = rs
dr  = rin/50.0d0
fac = 1.01d0

r(1) = rin

do i = 2, nrr
   r(i) = r(1) + (fac**(i-1) - 1.0d0)*dr/(fac - 1.0d0)
end do

do i = 1, nrr
write(*,*) i, r(i)
end do

! formatted means the file is human readable
open(69, file="all_rvalues.dat", form="formatted")

write(69, *) 'i, r(i), dr(i)'
write(69, *)
write(69, *) 1, rin, 0.0d0
do i = 2, nrr
    write(69,*) i, r(i),  (fac**(i-1) - 1.0d0)*dr/(fac - 1.0d0)
    write(69,*)
end do
close(69)

write(*,*) 'inner boundary:        ', rs
write(*,*) 'outer boundary:        ', r(nrr)
write(*,*) 'outer boundary:        ', r(nrr)/rs, ' rs.'
write(*,*)
write(*,*) '------------------------------------'

dt1 = 0.5d0*pi/real(nt+1,kind(1.0d0))
dmu = 1.0d0/real(nt,kind(1.0d0))

do j = nta, ntb
   mu(j) = 1.0d0 + 0.5d0*dmu - real(j-1,kind(1.0d0))*dmu
end do

! Two types of grid in the angular direction

do j = nta, ntb
   t(j) = acos(mu(j))
!    t(j) = real(j-1,kind(1.0d0))*dt1

  if (t(j) .lt. tcut) then
    jcut = j
  end if

end do

t(nta-1) = 0.0d0
t(ntb+1) = 0.5d0*pi

! Kerr quantities

do i = nra-1, nrb+1
do j = nta-1, ntb+1

! eqs. 40, 39, 38, 36, 37
  rk     = r(i)*(1.0d0 + kerrm/r(i) &
    + (kerrm**2 - kerra**2)/(4.0d0*r(i)**2))
  deltak = rk**2 - 2.0d0*kerrm*rk + kerra**2
  sigmak = rk**2 + kerra**2*cos(t(j))**2

  capitalhe(i,j) = kerrm*kerra*((rk**2 - kerra**2)*sigmak &
    + 2.0d0*rk**2*(rk**2 + kerra**2))
  capitalhe(i,j) = capitalhe(i,j)/sigmak**2

  capitalhf(i,j) = - 2.0d0*kerrm*kerra**3*rk*sqrt(deltak)*cos(t(j))*sin(t(j))**2
  capitalhf(i,j) = capitalhf(i,j)/sigmak**2

! This is the (optional) Kerr metric
! eqs. 32, 33, 35, 34

  alpha(i,j) = sigmak*deltak/((rk**2 + kerra**2)*sigmak &
    + 2.0d0*kerrm*kerra**2*rk*sin(t(j))**2)
  alpha(i,j) = sqrt(alpha(i,j))

  betak(i,j) = - 2.0d0*kerrm*kerra*rk/((rk**2 + kerra**2)*sigmak &
    + 2.0d0*kerrm*kerra**2*rk*sin(t(j))**2)

  q(i,j) = sigmak/sqrt((rk**2 + kerra**2)*sigmak &
    + 2.0d0*kerrm*kerra**2*rk*sin(t(j))**2)
  q(i,j) = log(q(i,j))

  psi(i,j) = rk**2 + kerra**2 + 2.0d0*kerrm*kerra**2*rk*sin(t(j))**2/sigmak
  psi(i,j) = psi(i,j)**0.25d0/sqrt(r(i))

! eq. 42

  phi(i,j) = log(psi(i,j)/(1.0d0 + rs/r(i)))

  capitalb(i,j) = 1.0d0

  b(i,j) = log(capitalb(i,j))

  betat(i,j) = 0.0d0

! eq 48

  integrand(i,j) = 2.0d0*capitalhe(i,j)*capitalb(i,j) &
    *exp(-8.0d0*phi(i,j))*(r(i) - rs)*r(i)**2/(r(i) + rs)**7

  test  (i,j) = betak(i,j)
  test2 (i,j) = phi(i,j)
  test3 (i,j) = q(i,j)

end do
end do

! Start with setting all matter-related terms to zero

rho(:,:)         = 0.0d0
h(:,:)           = 1.0d0
capitalp(:,:)    = 0.0d0
rhocapitalh(:,:) = 0.0d0
capitalj(:,:)    = 0.0d0

! For testing purposes only: save all boundary values from the initial solution

do j = nta-1, ntb+1
  boundary44(j,1) = phi      (nrb+1,j)
  boundary45(j,1) = capitalb (nrb+1,j)
  boundary46(j,1) = betat    (nrb+1,j)
  boundary47(j,1) = q        (nrb+1,j)
end do

do j = nta-1, ntb+1
  horizon44(j) = phi      (nra-1,j)
  horizon45(j) = capitalb (nra-1,j)
  horizon46(j) = betat    (nra-1,j)
  horizon47(j) = q        (nra-1,j)
end do

do i = nra-1, nrb+1
  axis44(i)    = phi      (i,nta-1)
  axis45(i)    = capitalb (i,nta-1)
  axis46(i)    = betat    (i,nta-1)
  axis47(i)    = q        (i,nta-1)

  equator44(i) = phi      (i,ntb+1)
  equator45(i) = capitalb (i,ntb+1)
  equator46(i) = betat    (i,ntb+1)
  equator47(i) = q        (i,ntb+1)
end do


jdiskcut = jcut

! Loading RESTART

if (restart.eq.1) then

open (20,file="RESTART.dat",form="unformatted")

read(20) phi
read(20) sphi
read(20) q
read(20) sq
read(20) capitalb
read(20) scapitalb
read(20) betat
read(20) betak
read(20) h
read(20) rho
read(20) capitalomega
read(20) b
read(20) omega1
read(20) omega2
read(20) kk
read(20) jdisk
read(20) jdiskcut

close(20)

endif

! RESTART loaded

! ----------------------------------------------------------------------------------------
! main iteration loop starts here

do iter = 1, niter

write(*,*)
write(*,*) 'Starting outer loop iteration iter  =', iter
write(*,*)

! compute b, psi, and alpha, alpha = phi\psi (above eq. 18)

  b(:,:) = log(capitalb(:,:))

do i = nra-1, nrb+1
do j = nta-1, ntb+1
  psi(i,j)   = (1.0d0 + rs/r(i))*exp(phi(i,j))
  alpha(i,j) = (1.0d0 - rs/r(i))*exp(-2.0d0*phi(i,j))*capitalb(i,j)/(1.0d0 + rs/r(i))
end do
end do


! These subroutines compute first derivatives of phi, b, betat with respect to
! r and theta. One can implement higher order differentiation in
! subroutines derivr and derivt.

call derivr(phi,   dphidrm,   r)
call derivt(phi,   dphidtm,   t)
call derivr(b,     dbdrm,     r)
call derivt(b,     dbdtm,     t)
call derivr(betat, dbetatdrm, r)
call derivt(betat, dbetatdtm, t)

! Set boundary conditions for each iteration. These values will be used
! by solver subroutines exact. Most of these formulas implement the
! condition corresponding to the vanishing of the gradient normal
! to the boundary.


do i = nra-1, nrb+1

  dt1 = t(ntb-1) - t(ntb-2)
  dt2 = t(ntb) - t(ntb-1)
  dt3 = t(ntb+1) - t(ntb)

!  equator44(i) = + ((dt3*(dt2 + dt3))/(dt1*(dt1 + dt2)*(dt1 + dt2 + dt3)))*phi(i,ntb-2) &
!    - (dt3*(dt1 + dt2 + dt3))/(dt1*dt2*(dt2 + dt3))*phi(i,ntb-1) &
!    + (((dt2 + dt3)*(dt1 + dt2 + dt3))/(dt2*(dt1 + dt2)*dt3))*phi(i,ntb)
!  equator44(i) = equator44(i)/((dt2**2 + 4.0d0*dt2*dt3 + 3.0d0*dt3**2 + dt1*(dt2 + 2.0d0*dt3))/ &
!    (dt3*(dt2 + dt3)*(dt1 + dt2 + dt3)))

! lower order

  equator44(i) = (-dt3**2/(dt2*(dt2 + 2.0d0*dt3)))*phi(i,ntb-1) + ((dt2+dt3)**2/(dt2*(dt2 + 2.0d0*dt3)))*phi(i,ntb)

!  equator45(i) = + ((dt3*(dt2 + dt3))/(dt1*(dt1 + dt2)*(dt1 + dt2 + dt3)))*capitalb(i,ntb-2) &
!    - (dt3*(dt1 + dt2 + dt3))/(dt1*dt2*(dt2 + dt3))*capitalb(i,ntb-1) &
!    + (((dt2 + dt3)*(dt1 + dt2 + dt3))/(dt2*(dt1 + dt2)*dt3))*capitalb(i,ntb)
!  equator45(i) = equator45(i)/((dt2**2 + 4.0d0*dt2*dt3 + 3.0d0*dt3**2 + dt1*(dt2 + 2.0d0*dt3))/ &
!    (dt3*(dt2 + dt3)*(dt1 + dt2 + dt3)))

  equator45(i) = (-dt3**2/(dt2*(dt2 + 2.0d0*dt3)))*capitalb(i,ntb-1) + ((dt2+dt3)**2/(dt2*(dt2 + 2.0d0*dt3)))*capitalb(i,ntb)

!  equator46(i) = + ((dt3*(dt2 + dt3))/(dt1*(dt1 + dt2)*(dt1 + dt2 + dt3)))*betat(i,ntb-2) &
!    - (dt3*(dt1 + dt2 + dt3))/(dt1*dt2*(dt2 + dt3))*betat(i,ntb-1) &
!    + (((dt2 + dt3)*(dt1 + dt2 + dt3))/(dt2*(dt1 + dt2)*dt3))*betat(i,ntb)
!  equator46(i) = equator46(i)/((dt2**2 + 4.0d0*dt2*dt3 + 3.0d0*dt3**2 + dt1*(dt2 + 2.0d0*dt3))/ &
!    (dt3*(dt2 + dt3)*(dt1 + dt2 + dt3)))

  equator46(i) = (-dt3**2/(dt2*(dt2 + 2.0d0*dt3)))*betat(i,ntb-1) + ((dt2+dt3)**2/(dt2*(dt2 + 2.0d0*dt3)))*betat(i,ntb)

!  equator47(i) = + ((dt3*(dt2 + dt3))/(dt1*(dt1 + dt2)*(dt1 + dt2 + dt3)))*q(i,ntb-2) &
!    - (dt3*(dt1 + dt2 + dt3))/(dt1*dt2*(dt2 + dt3))*q(i,ntb-1) &
!    + (((dt2 + dt3)*(dt1 + dt2 + dt3))/(dt2*(dt1 + dt2)*dt3))*q(i,ntb)
!  equator47(i) = equator47(i)/((dt2**2 + 4.0d0*dt2*dt3 + 3.0d0*dt3**2 + dt1*(dt2 + 2.0d0*dt3))/ &
!    (dt3*(dt2 + dt3)*(dt1 + dt2 + dt3)))

! lower order

  equator47(i) = (-dt3**2/(dt2*(dt2 + 2.0d0*dt3)))*q(i,ntb-1) + ((dt2+dt3)**2/(dt2*(dt2 + 2.0d0*dt3)))*q(i,ntb)

  dt1 = t(nta) - t(nta - 1)
  dt2 = t(nta + 1) - t(nta)
  dt3 = t(nta + 2) - t(nta + 1)

!  axis44(i) = -((dt1 + dt2)*(dt1 + dt2 + dt3))/(dt1*dt2*(dt2 + dt3))*phi(i,nta) &
!     +((dt1*(dt1 + dt2 + dt3))/(dt2*(dt1 + dt2)*dt3))*phi(i,nta+1) &
!     -(dt1*(dt1 + dt2))/(dt3*(dt2 + dt3)*(dt1 + dt2 + dt3))*phi(i,nta+2)
!  axis44(i) = axis44(i)/(-(1.0d0/dt1) - 1.0d0/(dt1 + dt2) - 1.0d0/(dt1 + dt2 + dt3))

! lower order

  axis44(i) = ((dt1 + dt2)**2/(dt2*(2*dt1 + dt2)))*phi(i,nta) - (dt1**2/(dt2*(2*dt1 + dt2)))*phi(i,nta+1)

!  axis45(i) = -((dt1 + dt2)*(dt1 + dt2 + dt3))/(dt1*dt2*(dt2 + dt3))*capitalb(i,nta) &
!     +((dt1*(dt1 + dt2 + dt3))/(dt2*(dt1 + dt2)*dt3))*capitalb(i,nta+1) &
!     -(dt1*(dt1 + dt2))/(dt3*(dt2 + dt3)*(dt1 + dt2 + dt3))*capitalb(i,nta+2)
!  axis45(i) = axis45(i)/(-(1.0d0/dt1) - 1.0d0/(dt1 + dt2) - 1.0d0/(dt1 + dt2 + dt3))

  axis45(i) = ((dt1 + dt2)**2/(dt2*(2*dt1 + dt2)))*capitalb(i,nta) - (dt1**2/(dt2*(2*dt1 + dt2)))*capitalb(i,nta+1)

!  axis46(i) = -((dt1 + dt2)*(dt1 + dt2 + dt3))/(dt1*dt2*(dt2 + dt3))*betat(i,nta) &
!     +((dt1*(dt1 + dt2 + dt3))/(dt2*(dt1 + dt2)*dt3))*betat(i,nta+1) &
!     -(dt1*(dt1 + dt2))/(dt3*(dt2 + dt3)*(dt1 + dt2 + dt3))*betat(i,nta+2)
!  axis46(i) = axis46(i)/(-(1.0d0/dt1) - 1.0d0/(dt1 + dt2) - 1.0d0/(dt1 + dt2 + dt3))

  axis46(i) = ((dt1 + dt2)**2/(dt2*(2*dt1 + dt2)))*betat(i,nta) - (dt1**2/(dt2*(2*dt1 + dt2)))*betat(i,nta+1)

! q should vanish at the axis

  axis47(i) = 0.0d0

end do

do j = nta-1, ntb+1
  d1 = r(nra) - r(nra - 1)
  d2 = r(nra + 1) - r(nra)
  d3 = r(nra + 2) - r(nra + 1)

  horizon44(j) = -((d1 + d2)*(d1 + d2 + d3))/(d1*d2*(d2 + d3))*phi(nra,j) &
     +((d1*(d1 + d2 + d3))/(d2*(d1 + d2)*d3))*phi(nra+1,j) &
     -(d1*(d1 + d2))/(d3*(d2 + d3)*(d1 + d2 + d3))*phi(nra+2,j)
  horizon44(j) = horizon44(j)/(-(1.0d0/d1) - 1.0d0/(d1 + d2) - 1.0d0/(d1 + d2 + d3))

!  horizon45(j) = -((d1 + d2)*(d1 + d2 + d3))/(d1*d2*(d2 + d3))*capitalb(nra,j) &
!     +((d1*(d1 + d2 + d3))/(d2*(d1 + d2)*d3))*capitalb(nra+1,j) &
!     -(d1*(d1 + d2))/(d3*(d2 + d3)*(d1 + d2 + d3))*capitalb(nra+2,j)
!  horizon45(j) = horizon45(j)/(-(1.0d0/d1) - 1.0d0/(d1 + d2) - 1.0d0/(d1 + d2 + d3))

  horizon45(j) = ((d1 + d2)**2/(d2*(2*d1 + d2)))*capitalb(nra,j) - (d1**2/(d2*(2*d1 + d2)))*capitalb(nra+1,j)

! betat should vanish at the horizon

  horizon46(j) = 0.0d0

!  horizon47(j) = -((d1 + d2)*(d1 + d2 + d3))/(d1*d2*(d2 + d3))*q(nra,j) &
!     +((d1*(d1 + d2 + d3))/(d2*(d1 + d2)*d3))*q(nra+1,j) &
!     -(d1*(d1 + d2))/(d3*(d2 + d3)*(d1 + d2 + d3))*q(nra+2,j)
!  horizon47(j) = horizon47(j)/(-(1.0d0/d1) - 1.0d0/(d1 + d2) - 1.0d0/(d1 + d2 + d3))

  horizon47(j) = ((d1 + d2)**2/(d2*(2*d1 + d2)))*q(nra,j) - (d1**2/(d2*(2*d1 + d2)))*q(nra+1,j)

end do

! These conditions are used in soubroutine exact to set the values of betat
! for i = nra, nra+1, nra+2. They stem from the condition that first, second and
! third derivatives of betat with respect to r vanish at r = rs.

do j = nta-1, ntb+1

  d1 = r(nra) - r(nra - 1)
  d2 = r(nra + 1) - r(nra)
  d3 = r(nra + 2) - r(nra + 1)
  d4 = r(nra + 3) - r(nra + 2)
  d5 = r(nra + 4) - r(nra + 3)

  betahor(j,1) = (d1**4*(((d2 + d3 + d4 + d5)*betat(nra+3,j))/(d1 + d2 + d3 + d4)**4 - & 
                 ((d2 + d3 + d4)*betat(nra+4,j))/(d1 + d2 + d3 + d4 + d5)**4))/d5 
  betahor(j,2) = ((d1 + d2)**4*(((d3 + d4 + d5)*betat(nra+3,j))/(d1 + d2 + d3 + d4)**4 - &
                 ((d3 + d4)*betat(nra+4,j))/(d1 + d2 + d3 + d4 + d5)**4))/d5
  betahor(j,3) = ((d1 + d2 + d3)**4*(((d4 + d5)*betat(nra+3,j))/(d1 + d2 + d3 + d4)**4 - &
                 (d4*betat(nra+4,j))/(d1 + d2 + d3 + d4 + d5)**4))/d5

end do

! compute the enthalpy within the disk

! set the guess for the inner and outer radii of the disk

r1 = rr1
r2 = rr2

i1 = nra-1
do i = nra-1, nrb+1
  if (r(i) .gt. r1) then
    i1 = i
    exit
  end if
end do

i2 = nra-1
do i = nra-1, nrb+1
  if (r(i) .gt. r2) then
    i2 = i
    exit
  end if
end do

r1    = r(i1)
r2    = r(i2)

alfa1 = alpha(i1,ntb+1)
alfa2 = alpha(i2,ntb+1)
psi1  = psi(i1,ntb+1)
psi2  = psi(i2,ntb+1)
beta1 = betak(i1,ntb+1) + betat(i1,ntb+1)
beta2 = betak(i2,ntb+1) + betat(i2,ntb+1)

call diskparams(c1, w, omega1, omega2, r1, r2, beta1, beta2, alfa1, alfa2, psi1, psi2, arotlaw)

if ((iter .eq. 1) .and. (restart .eq. 0)) then
   capitalomega(:,:) = 0.5d0*(omega1 + omega2)
end if

h(:,:) = 1.0d0
uphidn(:,:) = 0.0d0
utup(:,:) = 1.0d0

do i = nra-1, i1
do j = nta-1, ntb+1
  capitalomega(i,j) = omega1
end do
end do

do i = i1, i2
do j = nta - 1, jdiskcut(i) - 1
   capitalomega(i,j) = omega1
end do
do j = jdiskcut(i), ntb+1

  if (r(i)*sin(t(j)) .ge. r1) then

    call diskomega(capitalomega(i,j), r(i), t(j), alpha(i,j), &
       betat(i,j) + betak(i,j), psi(i,j), w, arotlaw)

    h(i,j) = (c1*sqrt(1.0d0 - arotlaw**2*capitalomega(i,j)**2 - & 
           3.0d0*capitalomega(i,j)**0.6666666666666666d0*(1.0d0 &
      - arotlaw*capitalomega(i,j))**1.3333333333333333d0*W))/ &
       sqrt(alpha(i,j)**2 - (betat(i,j) + betak(i,j))**2*psi(i,j)**4*r(i)**2*sin(t(j))**2 - &
         2.0d0*(betat(i,j)+betak(i,j))*capitalomega(i,j)*psi(i,j)**4*r(i)**2*sin(t(j))**2 - &
      capitalomega(i,j)**2*psi(i,j)**4*r(i)**2*sin(t(j))**2)

    if (h(i,j) .lt. 1.0d0) then
      h(i,j) = 1.0d0
    end if

    uphidn(i,j) = ((betak(i,j) + betat(i,j) + capitalomega(i,j))*psi(i,j)**4*r(i)**2*sin(t(j))**2)/&
      sqrt(alpha(i,j)**2 - (betak(i,j) + betat(i,j) + capitalomega(i,j))**2*psi(i,j)**4*r(i)**2*sin(t(j))**2)

    utup(i,j) = 1.0d0/sqrt(alpha(i,j)**2 - (betak(i,j) + betat(i,j) + &
                capitalomega(i,j))**2*psi(i,j)**4*r(i)**2*sin(t(j))**2)

  end if

end do
end do


tmax = 0.5d0*pi
jtmax = ntb+1
do i = nra, nrb+1
do j = nta, ntb+1
  if (h(i,j) .gt. 1.0d0) then
    if (t(j) .le. tmax) then
      tmax = t(j)
      jtmax = j
    end if
  end if
end do
end do

write(*,*) 'tmax/pi: ', tmax/pi
write(*,*) 'jtmax: ', jtmax, t(jtmax-1)/pi, t(jtmax)/pi, t(jtmax+1)/pi

if ((iter .eq. 1) .and. (restart .eq. 0)) then

! looking for the point of max(rho)

i3 = i1
j3 = nta-1
h0 = h(i1,j3)
do i = i1, i2
do j = nta-1, ntb+1
if (h(i, j) .gt. h0) then
i3 = i
j3 = j
h0 = h(i,j)
end if
end do
end do

! reset r3 to nearest grid point

r3 = r(i3)

write(*,*) 'h0: ', h0, 'r3: ', r3, 'i3: ', i3, 'j3: ', j3

! Compute the polytropic constant (eqs. 4, 5)

kk = (gam - 1.0d0)*(h(i3,j3) - 1.0d0)/(gam*rho0**(gam - 1.0d0))

end if

do i = nra, nrb+1
do j = nta, ntb+1
  if (h(i,j) .gt. 1.0d0) then
    call magbernoulli(h(i,j), alpha(i,j),psi(i,j),r(i),t(j),gam,kk,nn,c2)
  end if
end do
end do

! looking for the point of max(rho)

i3 = i1
j3 = nta-1
h0 = h(i1,j3)
do i = i1, i2
do j = nta-1, ntb+1
if (h(i, j) .gt. h0) then
i3 = i
j3 = j
h0 = h(i,j)
end if
end do
end do

! reset r3 to nearest grid point

r3 = r(i3)

write(*,*) 'h0: ', h0, 'r3: ', r3, 'i3: ', i3, 'j3: ', j3

! compute the range of the disk in j

jdisk(:) = ntb+1
do i = i1, i2
do j = jdiskcut(i), ntb+1
if ((h(i,j) .gt. 1.0d0) .and. (t(j)) .lt. t(jdisk(i))) then
jdisk(i) = j
end if
end do
end do


! Compute the polytropic constant (eqs. 4, 5)

kk = (gam - 1.0d0)*(h(i3,j3) - 1.0d0)/(gam*rho0**(gam - 1.0d0))


! Compute other matter-related quantities

do i = nra, nrb+1
do j = nta, ntb+1

rho(i,j) = ((h(i,j) - 1.0d0)*(gam - 1.0d0)/(kk*gam))**(1.0d0/(gam - 1.0d0))
capitalp(i,j) = kk*rho(i,j)**gam

if (h(i,j) .gt. 1.0d0) then
  b2(i,j) = rho(i,j)*h(i,j)*alpha(i,j)**2*psi(i,j)**4*r(i)**2*sin(t(j))**2
  b2(i,j) = 2.0d0*nn*(b2(i,j) - log(1.0d0 + c2*b2(i,j))/c2)
  b2(i,j) = b2(i,j)/(alpha(i,j)**2*psi(i,j)**4*r(i)**2*sin(t(j))**2)
else
  b2(i,j) = 0.0d0
end if

! compute magnetisation parameter

betamag = capitalp(i3,j3)/(0.5d0*b2(i3,j3))

! eq. 22, 23 but for a different rotation law

!rhocapitalh(i,j) = rho(i,j)*h(i,j) + &
!   &  (rho(i,j)/h(i,j))*jj**2/(psi(i,j)**4*r(i)**2*sin(t(j))**2) &
!   &  - capitalp(i,j)

rhocapitalh(i,j) = rho(i,j)*h(i,j)*alpha(i,j)**2*utup(i,j)**2 - capitalp(i,j) + 0.5d0*b2(i,j)

!capitalj(i,j) = rho(i,j)*jj*sqrt(1.0d0 + jj**2/(h(i,j)**2*psi(i,j)**4*r(i)**2*sin(t(j))**2))

capitalj(i,j) = rho(i,j)*h(i,j)*alpha(i,j)*utup(i,j)*uphidn(i,j)

end do
end do

! compute integrand = r.h.s. of Eq. (48)

do i = nra-1, nrb+1
do j = nta-1, ntb+1

  integrand(i,j) = 2.0d0*capitalhe(i,j)*capitalb(i,j) &
    *exp(-8.0d0*phi(i,j))*(r(i) - rs)*r(i)**2/(r(i) + rs)**7

end do
end do

! compute A^2 = capitala2, eq. (21)

do i = nra-1, nrb+1
do j = nta-1, ntb+1

! coding Eqs. 20,21 and 14,15

   dbetatdr = dbetatdrm(i,j)
   dbetatdt = dbetatdtm(i,j)

   capitala2(i,j) = sin(t(j))**2*(capitalhe(i,j)/r(i)**3 &
    + psi(i,j)**6*r(i)*dbetatdr/(2.0d0*alpha(i,j)))**2

   capitala2(i,j) = capitala2(i,j) + (capitalhf(i,j)/r(i)**3 &
    + psi(i,j)**6*sin(t(j))*dbetatdt/(2.0d0*alpha(i,j)))**2

!   capitala2(i,j) = sin(t(j))**2*(capitalhe(i,j)/r(i)**3)**2
!   capitala2(i,j) = capitala2(i,j) + (capitalhf(i,j)/r(i)**3)**2

end do
end do

! Correcting divergent terms. This is a polynomial extrapolation to r = rs

dr1 = r(nra) - r(nra - 1)
dr2 = r(nra + 1) - r(nra)
dr3 = r(nra + 2) - r(nra + 1)

do j = nta-1, ntb+1

  capitala2(nra-1,j) = ((dr1 + dr2)*dr3*(dr1 + dr2 + dr3)*capitala2(nra,j) - & 
       dr1*(dr2 + dr3)*(dr1 + dr2 + dr3)*capitala2(nra+1,j) + dr1*dr2*(dr1 + dr2)*capitala2(nra+2,j))/ &
       (dr2*dr3*(dr2 + dr3))

end do


! Computing source terms in Eqs. (44), (45), (46) and (47)

! Note that solving Eqs. (44), (45), (46) and (47) actually requires the source terms
! only in the interior of the grid. The values of the source terms at the boundaries
! are only needed by the procedures that integrate the sources in order to obtain
! constants q1, capitalb1, capitalm1, and capitalj1. These values are often multiplied
! by 0, so not always a true value is required.

! computing sphi

do i = nra-1, nrb+1
do j = nta, ntb+1

   dphidr = dphidrm(i,j)
   dphidt = dphidtm(i,j)

   dbdr = dbdrm(i,j)
   dbdt = dbdtm(i,j)

   dbetatdr = dbetatdrm(i,j)
   dbetatdt = dbetatdtm(i,j)

  temp = (r(i) - rs)*dbdr/(r(i)*(r(i) + rs)) &
              + cos(t(j))*dbdt/sin(t(j))/r(i)**2

  sphi(i,j) = - capitala2(i,j)/psi(i,j)**8 - dphidr*dbdr - dphidt*dbdt/r(i)**2
  sphi(i,j) = sphi(i,j) - 0.5d0*temp

!  sphi(i,j) = sphi(i,j) - 2.0d0*pi*exp(2.0d0*q(i,j))*psi(i,j)**4*(rhocapitalh(i,j) &
!     - capitalp(i,j) + rho(i,j)*jj**2/(h(i,j)*psi(i,j)**4*r(i)**2*sin(t(j))**2))

  sphi(i,j) = sphi(i,j) - 2.0d0*pi*exp(2.0d0*q(i,j))*psi(i,j)**4*(rhocapitalh(i,j) &
      - capitalp(i,j) + rho(i,j)*h(i,j)*uphidn(i,j)**2/(psi(i,j)**4*r(i)**2*sin(t(j))**2) &
      - 3.0d0*b2(i,j)/2.0d0)

!   sphi(i,j) = - capitala2(i,j)/psi(i,j)**8

end do
end do

! Correcting divergent terms. Here again a polynomial extrapolation is used.

dt1 = t(nta) - t(nta - 1)
dt2 = t(nta + 1) - t(nta)
dt3 = t(nta + 2) - t(nta + 1)

do i = nra-1, nrb+1

  sphi(i,nta-1) = ((dt1 + dt2)*dt3*(dt1 + dt2 + dt3)*sphi(i,nta) - &
       dt1*(dt2 + dt3)*(dt1 + dt2 + dt3)*sphi(i,nta+1) + dt1*dt2*(dt1 + dt2)*sphi(i,nta+2))/ &
       (dt2*dt3*(dt2 + dt3))

end do

! Computing scapitalb

do i = nra-1, nrb+1
do j = nta-1, ntb+1

  scapitalb(i,j) = 16.0d0*pi*capitalb(i,j)*exp(2.0d0*q(i,j))*psi(i,j)**4*(capitalp(i,j) + 0.5d0*b2(i,j))

end do 
end do

! Computing sbetat

do i = nra-1, nrb+1
do j = nta-1, ntb+1

  dphidr = dphidrm(i,j)
  dphidt = dphidtm(i,j)

  dbdr = dbdrm(i,j)
  dbdt = dbdtm(i,j)

  dbetatdr = dbetatdrm(i,j)
  dbetatdt = dbetatdtm(i,j)

! No rotating matter allowed at the axis

  if (j .eq. nta-1) then
    temp = 0.0d0
  else
    temp = 16.0d0*pi*alpha(i,j)*exp(2.0d0*q(i,j))*capitalj(i,j)/(r(i)**2*sin(t(j))**2)
  end if

  sbetat(i,j)    = temp &
                   - 8.0d0*dphidr*dbetatdr + dbdr*dbetatdr - 8.0d0*dphidt*dbetatdt/r(i)**2 &
                   + dbdt*dbetatdt/r(i)**2

end do
end do

! Computing sq

do i = nra, nrb+1
do j = nta, ntb+1

   dphidr = dphidrm(i,j)
   dphidt = dphidtm(i,j)

   dbdr = dbdrm(i,j)
   dbdt = dbdtm(i,j)

   dbetatdr = dbetatdrm(i,j)
   dbetatdt = dbetatdtm(i,j)

  temp = (r(i) - rs)*dbdr/(r(i)*(r(i) + rs)) &
              + cos(t(j))*dbdt/sin(t(j))/r(i)**2

  sq(i,j) = 3.0d0*capitala2(i,j)/psi(i,j)**8 + 2.0d0*temp &
     + (8.0d0*rs/(r(i)**2 - rs**2) + 4.0d0*(dbdr - dphidr))*dphidr &
     + 4.0d0*dphidt*(dbdt - dphidt)/r(i)**2

!  sq(i,j) = sq(i,j) - 8.0d0*pi*exp(2.0d0*q(i,j))*(psi(i,j)**4*capitalp(i,j) &
!            - rho(i,j)*jj**2/(h(i,j)*r(i)**2*sin(t(j))**2))

   sq(i,j) = sq(i,j) - 8.0d0*pi*exp(2.0d0*q(i,j))*(psi(i,j)**4*capitalp(i,j) &
            - rho(i,j)*h(i,j)*uphidn(i,j)**2/(r(i)**2*sin(t(j))**2) &
            + 3.0d0*psi(i,j)**4*b2(i,j)/2.0d0)

!   dbdr = 0.0d0
!   dbdt = 0.0d0

!   sq(i,j) = 3.0d0*capitala2(i,j)/psi(i,j)**8 + (8.0d0*rs/(r(i)**2 - rs**2) + 4.0d0*(dbdr - dphidr))*dphidr &
!     + 4.0d0*dphidt*(dbdt - dphidt)/r(i)**2

end do
end do

! Correcting divergent terms. Polynomial extrapolation is used.

dr1 = r(nra) - r(nra - 1)
dr2 = r(nra + 1) - r(nra)
dr3 = r(nra + 2) - r(nra + 1)

do j = nta-1, ntb+1

  sq(nra-1,j) = ((dr1 + dr2)*dr3*(dr1 + dr2 + dr3)*sq(nra,j) - & 
       dr1*(dr2 + dr3)*(dr1 + dr2 + dr3)*sq(nra+1,j) + dr1*dr2*(dr1 + dr2)*sq(nra+2,j))/ &
       (dr2*dr3*(dr2 + dr3))

end do

dt1 = t(nta) - t(nta - 1)
dt2 = t(nta + 1) - t(nta)
dt3 = t(nta + 2) - t(nta + 1)

do i = nra-1, nrb+1

  sq(i,nta-1) = ((dt1 + dt2)*dt3*(dt1 + dt2 + dt3)*sq(i,nta) - &
       dt1*(dt2 + dt3)*(dt1 + dt2 + dt3)*sq(i,nta+1) + dt1*dt2*(dt1 + dt2)*sq(i,nta+2))/ &
       (dt2*dt3*(dt2 + dt3))

end do


! Computing coefficients in boundary conditions

! Compute capitalm1

do i = nra-1, nrb+1
do j = nta-1, ntb+1
  storage(i,j) = (r(i)**2 - rs**2)*sin(t(j))*sphi(i,j)
end do
end do

call integration2d (storage, capitalm1, r, t, storage2)

capitalm1 = -2.0d0*capitalm1

! Compute capitalj1

do i = nra-1, nrb+1
do j = nta, ntb+1
!  storage(i,j) = r(i)**2*sin(t(j))*rho(i,j)*psi(i,j)**6*exp(2.0d0*q(i,j))*jj*sqrt(1.0d0 &
!    + jj**2/(h(i,j)**2*psi(i,j)**4*r(i)**2*sin(t(j))**2))

storage(i,j) = r(i)**2*sin(t(j))*rho(i,j)*h(i,j)*alpha(i,j)* &
               utup(i,j)*psi(i,j)**6*exp(2.0d0*q(i,j))*uphidn(i,j)  

end do
end do

do i = nra-1, nrb+1
  storage(i,nta-1) = 0.0d0
end do

call integration2d (storage, capitalj1, r, t, storage2)

capitalj1 = 4.0d0*pi*capitalj1

! Compute capitalb1

do i = nra-1, nrb+1
do j = nta-1, ntb+1
  storage(i,j) = (r(i)**2 - rs**2)**2*sin(t(j))**2*scapitalb(i,j)/r(i)
end do
end do

call integration2d (storage, capitalb1, r, t, storage2)
 
capitalb1 = 2.0d0*capitalb1/pi

! Compute q1

do j = nta-1, ntb+1
  storage3(j) = cos(2.0d0*t(j))*q(nra-1,j)
end do

call integration1d (storage3, temp, t)

do i = nra-1, nrb+1
do j = nta-1, ntb+1
  storage(i,j) = r(i)**3*cos(2.0d0*t(j))*sq(i,j)
end do
end do

call integration2d (storage, q1, r, t, storage2)

q1 = (2.0d0/pi)*q1 - (4.0d0*rs**2/pi)*temp

write(*,*) 'r1:       ', r1, 'i1:   ', i1
write(*,*) 'r2:       ', r2, 'i2:   ', i2
write(*,*) 'r1C:      ', r1*psi(i1,ntb+1)**2
write(*,*) 'r2C:      ', r2*psi(i2,ntb+1)**2
!write(*,*) 'j:        ', jj, 'jj2:  ', jj2
!write(*,*) 'C1:       ', cc, 'C2:   ', cc2
write(*,*) 'k:        ', kk
write(*,*) 'epsilon_0 ', (h(i3,ntb+1) - 1.0d0)/gam
write(*,*) 'M1:       ', capitalm1
write(*,*) 'J1:       ', capitalj1
write(*,*) 'B1:       ', capitalb1
write(*,*) 'q1        ', q1
write(*,*) 'stala 1d  ', -(4.0d0*rs**2/pi)*temp

! This is not necessary for boundary conditions, but I would like to have M* computed here

do i = nra-1, nrb+1
do j = nta, ntb+1
!  storage(i,j) = r(i)**2*sin(t(j))*rho(i,j)*psi(i,j)**6*exp(2.0d0*q(i,j))*sqrt(1.0d0 &
!    + jj**2/(h(i,j)**2*psi(i,j)**4*r(i)**2*sin(t(j))**2))

   storage(i,j) = r(i)**2*sin(t(j))*rho(i,j)*psi(i,j)**6*exp(2.0d0*q(i,j))*alpha(i,j)*utup(i,j)
end do
end do

do i = nra-1, nrb+1
  storage(i,nta-1) = 0.0d0
end do

call integration2d (storage, capitalmast, r, t, storage2)
capitalmast = 4.0d0*pi*capitalmast

write(165,*) capitalmast, kk

! ...and the other expression for M*, suggested by prof. Malec (without \alpha and u^t in \rho* (eq. 90) )
! \alpha*u^t could be obtained from eq. 6 & eq. 9

do i = nra-1, nrb+1
do j = nta, ntb+1
storage(i,j) = r(i)**2*sin(t(j)) * rho(i,j)*psi(i,j)**6*exp(2.0d0*q(i,j))
end do
end do

do i = nra-1, nrb+1
storage(i,nta-1) = 0.0d0
end do

call integration2d (storage, capitalmast2, r, t, storage2)
capitalmast2 = 4.0d0*pi*capitalmast2

! J* (eq. 89)

do i = nra-1, nrb+1
do j = nta, ntb+1
storage(i,j) = r(i)**2*sin(t(j)) * alpha(i,j)* rho(i,j)*psi(i,j)**6*exp(2.0d0*q(i,j)) * h(i,j)*uphidn(i,j)*utup(i,j)
end do
end do

do i = nra-1, nrb+1
storage(i,nta-1) = 0.0d0
end do

call integration2d (storage, capitaljast, r, t, storage2)
capitaljast = 4.0d0*pi*capitaljast

! also not necessary: compute other values characterizing the solution

! computing kappa from Eq. (76)
! this should be a constant

do j = nta-1, ntb+1
  storage3(j) = capitalB(nra-1,j)*exp(-4.0d0*phi(nra-1,j) - q(nra-1,j))/(8.0d0*rs)
end do

! take a mean value

call integration1d(storage3, kappa, t)
kappa = kappa/(pi/2.0d0)

! computing A_H (eq.72)

do j = nta-1, ntb + 1
  storage3(j) = exp(4.0d0*phi(nra-1,j) + q(nra-1,j))*sin(t(j))
end do

call integration1d (storage3, ah, t)
ah = ah*rs**2*16.0d0*4.0d0*pi

! computing c_p (eq. 73)

do j = nta-1, ntb + 1
storage3(j) = exp(2.0d0*phi(nra-1,j) + q(nra-1,j))
end do

call integration1d (storage3, cp, t)
cp = cp * 4.0d0 * 4.0d0 * rs

! computing c_e (eq. 74)

ce = 2.0d0 * pi * rs * psi(nra-1,ntb)**2

! computing M_C (eq. 84)

mc = ce / (4.0d0 * pi)


! compute M_irr and M_BH (eq. 87)

mirr = sqrt(ah/(16.0d0*pi))

mbh = mirr*sqrt(1.0d0 + kerrm**2*kerra**2/(4.0d0*mirr**4))

! compute Omega_H, J_H (eq. 75 & eq. 71)

do j = nta-1, ntb+1
  storage3(j) = -betak(nra-1,j)
end do

! take a mean value

call integration1d(storage3, omegah, t)
omegah = omegah/(pi/2.0d0)

jh = kerrm*kerra

! compute M_H

mh = kappa*ah/(4.0d0*pi) + 2.0d0*omegah*jh

! compute M_H once again

do j = nta-1, ntb+1
  storage3(j) = psi(nra-1,j)**4*capitalb(nra-1,j)*exp(-4.0d0*phi(nra-1,j))*sin(t(j))
end do

call integration1d(storage3, mh2, t)

mh2 = mh2*rs/8.0d0 + 2.0d0*omegah*jh

! compute M_ADM

madm = sqrt(kerrm**2 - kerra**2) + capitalm1

! compute M_T

storage(:,:) = 0.0d0

do i = nra, nrb + 1
  do j = nta, ntb + 1
    storage(i,j) = - 0.5d0*rho(i,j)*h(i,j) + 2.0d0*capitalp(i,j) + rhocapitalh(i,j) &
                   - (betat(i,j) + betak(i,j))*capitalj(i,j)/alpha(i,j)
    storage(i,j) = storage(i,j)*alpha(i,j)*psi(i,j)**6*exp(2.0d0*q(i,j))*sin(t(j))*r(i)**2
    storage(i,j) = 8.0d0*pi*storage(i,j)
  end do
end do

call integration2d (storage, mt, r, t, storage2)


write(*,*) 'Mstar:           ', capitalmast
write(*,*) 'Mstar2:          ', capitalmast2
write(*,*) 'M_ADM:           ', madm
write(*,*) 'M_H:             ', mh
write(*,*) 'M_H2:            ', mh2
write(*,*) 'M_H - MH2:       ', mh - mh2
write(*,*) 'M_BH:            ', mbh
write(*,*) 'M_irr:           ', mirr
write(*,*) 'M_T:             ', mt
write(*,*) 'Jstar:           ', capitaljast, capitaljast/sqrt(r2)
write(*,*) 'J_H:             ', jh
write(*,*) 'Omega_H:         ', omegah
write(*,*) 'kappa:           ', kappa
write(*,*) 'ah:              ', ah
write(*,*) 'Accuracy test:   ', abs((madm - mh  - mt)/madm)
write(*,*) 'Accuracy test 2: ', abs((madm - mh2 - mt)/madm)
write(*,*) 'mbh, w^2, madm:  ', mbh, w**(3.0d0/2.0d0), madm   
write(*,*) 'betamag:         ', betamag    

! computing ISCO

do i = nra-1, nrb+1
  gttf(i) = exp(4*phi(i,1 + ntb))*(betak(i,1 + ntb) + betat(i,1 + ntb))**2*(1 + rs/r(i))**4*r(i)**2 &
     -(capitalb(i,1 + ntb)**2*(-rs + r(i))**2)/(exp(4.0d0*phi(i,1 + ntb))*(rs + r(i))**2)
  gtphi(i) = exp(4*phi(i,1 + ntb))*(betak(i,1 + ntb) + betat(i,1 + ntb))*(1 + rs/r(i))**4*r(i)**2
  gphiphi(i) = exp(4*phi(i,1 + ntb))*(1 + rs/r(i))**4*r(i)**2
end do

call derivr1(gttf, drgttf, r)
call derivr1(drgttf, ddrgttf, r)
call derivr1(gtphi, drgtphi, r)
call derivr1(drgtphi, ddrgtphi, r)
call derivr1(gphiphi, drgphiphi, r)
call derivr1(drgphiphi, ddrgphiphi, r)

do i = nra-1, nrb+1
warisco(i)=(2.0d0*(drgtphi(i)**2 - drgphiphi(i)*drgttf(i))* &
       (-2.0d0*drgtphi(i)**2*gphiphi(i) + drgphiphi(i)*drgttf(i)*gphiphi(i) + &
       2.0d0*drgphiphi(i)*drgtphi(i)*gtphi(i) + &
       2.0d0*Sqrt(drgtphi(i)**2 - drgphiphi(i)*drgttf(i))*(drgtphi(i)*gphiphi(i) - drgphiphi(i)*gtphi(i)) - &
       drgphiphi(i)**2*gttf(i)) - (ddrgttf(i)*drgphiphi(i)**2 + &
       2.0d0*ddrgtphi(i)*drgphiphi(i)*(-drgtphi(i) + Sqrt(drgtphi(i)**2 - drgphiphi(i)*drgttf(i))) + &
       ddrgphiphi(i)*(2.0d0*drgtphi(i)**2 - drgphiphi(i)*drgttf(i) - &
       2.0d0*drgtphi(i)*Sqrt(drgtphi(i)**2 - drgphiphi(i)*drgttf(i))))*(gtphi(i)**2 - &
       gphiphi(i)*gttf(i)))
end do

!varisco = 0

do i = nra-1, nrb+1
  if (warisco(nrb+1-i) .ge. 0) then
    iisco = nrb+1-i
  else
 !   varisco = 1
    exit
  end if
end do

risco = r(iisco)
riscoc = risco*psi(iisco,ntb+1)**2

write(*,*) 'risco:          ', risco, iisco
write(*,*) 'riscoc:         ', riscoc





! Setting outer boundary values

do j = nta-1, ntb+1
   boundary44(j,1) = capitalm1/(2.0d0*r(nrb+1))
   boundary45(j,1) = 1.0d0 - capitalb1/r(nrb+1)**2
   boundary46(j,1) = -2.0d0*capitalj1/r(nrb+1)**3
   boundary47(j,1) = q1*sin(t(j))**2/r(nrb+1)**2
end do

! Computing coefficients for Eqs. (44), (45), (46) and (47)

do i = nra, nrb
do j = nta, ntb

  a44(i,j) = 2.0d0*r(i)/(r(i)**2 - rs**2)
  b44(i,j) = 1.0d0/r(i)**2
  c44(i,j) = cos(t(j))/(sin(t(j))*r(i)**2)

  a45(i,j) = (3.0d0*r(i)**2 + rs**2)/(r(i)*(r(i)**2 - rs**2))
  b45(i,j) = 1.0d0/r(i)**2
  c45(i,j) = 2.0d0*cos(t(j))/(sin(t(j))*r(i)**2)

  a46(i,j) = (4.0d0*r(i)**2 - 8.0d0*rs*r(i) &
             + 2.0d0*rs**2)/(r(i)*(r(i)**2 - rs**2))
  b46(i,j) = 1.0d0/r(i)**2
  c46(i,j) = (3.0d0*cos(t(j))/sin(t(j)))/r(i)**2

  a47(i,j) = 1.0d0/r(i)
  b47(i,j) = 1.0d0/r(i)**2
  c47(i,j) = 0.0d0

end do
end do

! This output serves for testing purposes. It is a possibility to store actual
! values of different quantities at each iteration.

do i = nra-1, nrb+1
do j = nta-1, ntb+1
  temp1 = -alpha(i,j)**2 + r(i)**2*(betat(i,j)+betak(i,j))**2*psi(i,j)**4*sin(t(j))**2
  temp2 = r(i)**2*(betat(i,j)+betak(i,j))*psi(i,j)**4*sin(t(j))**2
  temp3 = r(i)**2*psi(i,j)**4*sin(t(j))**2
  temp4 = (arotlaw**2*capitalomega(i,j)**1.3333333333333333d0 + &
     &    (1.0d0 - arotlaw*capitalomega(i,j))**0.3333333333333333d0*W - &
     &    3.0d0*arotlaw*capitalomega(i,j)*(1.0d0 - arotlaw*capitalomega(i,j))**0.3333333333333333d0*W)/ &
     &  (capitalomega(i,j)**0.3333333333333333d0* &
     &    (1.0d0 - arotlaw**2*capitalomega(i,j)**2 - &
     &      3.0d0*capitalomega(i,j)**0.6666666666666666d0* &
     &       (1.0d0 - arotlaw*capitalomega(i,j))**1.3333333333333333d0*W))
  test(i,j) = temp4*(-temp1 - 2.0d0*temp2*capitalomega(i,j) - temp3*capitalomega(i,j)**2) &
    - temp3*capitalomega(i,j) - temp2
  if (h(i,j) .le. 1.0d0) then
    test(i,j) = 0.0d0
  end if
end do
end do

if (mod(iter,100) .eq. 0) then

write(outnum,10) iter
outnum = adjustl(outnum)
outfil = 'wynik'//trim(outnum)//'.dat'

open(107,file=outfil,form="formatted")

!do i = nra-1, nrb+1
!do j = nta-1, ntb+1
!  write(107,*) r(i)*sin(t(j)), r(i)*cos(t(j)), capitalomega(i,j), h(i,j), r(i)*sin(t(j))*psi(i,j)**2
!end do
!  write(107,*)
!end do

!do i = nra-1, nrb+1
!  write(107,*) r(i), r(i)*psi(i,ntb+1)**2, warisco(i)
!end do

do i = nra-1, nrb+1
  write(107,*) r(i), r(i)*psi(i,ntb+1)**2, capitalp(i,ntb+1), 0.5d0*b2(i,ntb+1)
end do

close(107)

end if

! Solving Eqs. (44-47)

if (iter .eq. 1) then
 istart = .true.
else
 istart = .false.
endif

!$omp parallel
!$omp sections
!$omp section
call exact(phi, sphi, &
             a44, b44, &
             c44, r, t, &
             boundary44(:,1), nr+2, nt+2, ab44, IPIV44, istart,'p', &
             axis44, equator44, horizon44, betahor)
!$omp section
call exact(capitalb, scapitalb, &
             a45, b45, &
             c45, r, t, &
             boundary45(:,1), nr+2, nt+2, ab45, IPIV45, istart,'b', &
             axis45, equator45, horizon45, betahor)
!$omp section
call exact(betat, sbetat, &
             a46, b46, &
             c46, r, t, &
             boundary46(:,1), nr+2, nt+2, ab46, IPIV46, istart,'s', &
             axis46, equator46, horizon46, betahor)
!$omp section
call exact(q, sq, &
             a47, b47, &
             c47, r, t, &
             boundary47(:,1), nr+2, nt+2, ab47, IPIV47, istart,'q', &
             axis47, equator47, horizon47, betahor)
!$omp end sections
!$omp end parallel

! Solve Eq. (48)

call solvebetak(betak, integrand, r)

! adjusting the area in which capitalomega is computed

select case (adjustjdiskcut)

case(1:)
do i = i1, i2
if (jdiskcut(i) .gt. jdisk(i)-adjustjdiskcut) then
do n = 1, jdiskcut(i)+adjustjdiskcut-jdisk(i)
if(r(i)*sin(t(jdiskcut(i)-n)) .ge. r1) then
capitalomega(i,jdiskcut(i)-n) = capitalomega(i,jdiskcut(i))
end if
end do
end if
end do

jdiskcut(:) = jdisk(:) - adjustjdiskcut

case(:-1)

do i = i1, i2
if (jdiskcut(i) .gt. maxval(jdisk)+adjustjdiskcut) then
do n = 1, jdiskcut(i)-adjustjdiskcut-maxval(jdisk)
if(r(i)*sin(t(jdiskcut(i)-n)) .ge. r1) then
capitalomega(i,jdiskcut(i)-n) = capitalomega(i,jdiskcut(i))
end if
end do
end if
end do

jdiskcut(:) = maxval(jdisk) + adjustjdiskcut

end select



end do

! Main iteration loop ends here
! ----------------------------------------------------------------------------------------

open(22, file = 'x.dat', form="formatted", position='append')

write(22,*) 'input:           '
write(22,*) 'M_C:             ', kerrm
write(22,*) 'a:               ', kerra
write(22,*) 'omega1:          ', omega1
write(22,*) 'omega2:          ', omega2
write(22,*) 'gamma:           ', gam
write(22,*) 'grid (nr x nt):  ', nr, nt
write(22,*) 'r1:              ', rr1
write(22,*) 'r2:              ', rr2
write(22,*) 'rho0:            ', rho0
write(22,*) 'restart?:        ', restart
write(22,*) 'niter:           ', niter
write(22,*)
write(22,*) 'output:          '
write(22,*) 'r1:              ', r1, 'i1:   ', i1
write(22,*) 'r2:              ', r2, 'i2:   ', i2
write(22,*) 'k:               ', kk
write(22,*) 'epsilon_0        ', (h(i3,ntb+1) - 1.0d0)/gam
write(22,*) 'M1:              ', capitalm1
write(22,*) 'J1:              ', capitalj1
write(22,*) 'B1:              ', capitalb1
write(22,*) 'q1               ', q1
write(22,*) 'stala 1d         ', -(4.0d0*rs**2/pi)*temp
write(22,*) 'Mstar:           ', capitalmast
write(22,*) 'Mstar2:          ', capitalmast2
write(22,*) 'M_ADM:           ', madm
write(22,*) 'M_H:             ', mh
write(22,*) 'M_H2:            ', mh2
write(22,*) 'M_H - MH2:       ', mh - mh2
write(22,*) 'M_BH:            ', mbh
write(22,*) 'M_irr:           ', mirr
write(22,*) 'M_T:             ', mt
write(22,*) 'Jstar:           ', capitaljast
write(22,*) 'J_H:             ', jh
write(22,*) 'Omega_H:         ', omegah
write(22,*) 'kappa:           ', kappa
write(22,*) 'ah:              ', ah
write(22,*) 'Accuracy test:   ', abs((madm - mh  - mt)/madm)
write(22,*) 'Accuracy test 2: ', abs((madm - mh2 - mt)/madm)
write(22,*) '----------------------------------------------------------------'
close(22)

open (20,file="RESTART.dat",form="unformatted")

write(20) phi
write(20) sphi
write(20) q
write(20) sq
write(20) capitalb
write(20) scapitalb
write(20) betat
write(20) betak
write(20) h
write(20) rho
write(20) capitalomega
write(20) b
write(20) omega1
write(20) omega2
write(20) kk
write(20) jdisk
write(20) jdiskcut

close(20)

! Here we write the results to files

write(*,*)
write(*,*) 'End of iterations'
write(*,*)
write(*,*) '------------------------------------'
write(*,*)

! values for calculating ALL necessary variables
open(117, file="all_vars.dat", form="formatted")

! Loop over i and j
do i = nra-1, nrb+1
    do j = nta-1, ntb+1
        ! Write values
        write(117,*) r(i)*sin(t(j)), r(i)*cos(t(j)), r(i), t(j), alpha(i,j), psi(i,j), q(i,j), &
                    betak(i,j), betat(i,j), dbetatdrm(i,j), dbetatdtm(i,j), &
                    capitalhe(i,j), capitalhf(i,j), &
                    capitalomega(i,j), b2(i,j), rho(i,j) 
                    
    end do
    ! Write a blank line after each row
    write(117,*)
end do

! Close file
close(117)



! Deallocating storage space

!deallocate(wsp44 )
!deallocate(wsp45 )
!deallocate(wsp46 )
!deallocate(wsp47 )
!deallocate(IPIV44)

deallocate(ab44)
deallocate(ab45)
deallocate(ab46)
deallocate(ab47)

deallocate(IPIV45)
deallocate(IPIV46)
deallocate(IPIV47)

deallocate(axis44 )
deallocate(equator44)
deallocate(horizon44)

deallocate(axis45 )
deallocate(equator45)
deallocate(horizon45)

deallocate(axis46 )
deallocate(equator46)
deallocate(horizon46)

deallocate(axis47 )
deallocate(equator47)
deallocate(horizon47)

deallocate(betahor  )

deallocate(phi      )
deallocate(capitalb )
deallocate(betat    )
deallocate(q        )
deallocate(betak    )
deallocate(psi      )
deallocate(sphi     )
deallocate(scapitalb)
deallocate(sbetat   )
deallocate(sq       )
deallocate(alpha    )
deallocate(capitalhe)
deallocate(capitalhf)
deallocate(b        )
deallocate(capitala2)

deallocate(dphidrm  )
deallocate(dphidtm  )
deallocate(dbdrm    )
deallocate(dbdtm    )
deallocate(dbetatdrm)
deallocate(dbetatdtm)

deallocate(a44      )
deallocate(b44      )
deallocate(c44      )
deallocate(a45      )
deallocate(b45      )
deallocate(c45      )
deallocate(a46      )
deallocate(b46      )
deallocate(c46      )
deallocate(a47      )
deallocate(b47      )
deallocate(c47      )

deallocate(h        )
deallocate(rho        )
deallocate(rhocapitalh)
deallocate(capitalj   )
deallocate(capitalp   )


deallocate(storage  )
deallocate(integrand)
deallocate(test     )
deallocate(test2    )
deallocate(test3    )

deallocate(stor1    )
deallocate(stor2    )
deallocate(stor3    )
deallocate(stor4    )

deallocate(r       )
deallocate(storage2)
deallocate(t       )
deallocate(mu      )
deallocate(storage3)

deallocate(boundary44 )
deallocate(boundary45 )
deallocate(boundary46 )
deallocate(boundary47 )

deallocate(capitalomega)
deallocate(uphidn      )
deallocate(utup        )

deallocate(b2          )

deallocate(jdisk)
deallocate(jdiskcut)

deallocate(gttf       )
deallocate(gtphi      )
deallocate(gphiphi    )
deallocate(drgttf     )
deallocate(ddrgttf    )
deallocate(drgtphi    )
deallocate(ddrgtphi   )
deallocate(drgphiphi  )
deallocate(ddrgphiphi )
deallocate(warisco    )

close(165)

write(*,*) 'FINISHED!'

10  format (i5.5)

end program torusy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This subroutine solves Eq. (48)

subroutine solvebetak(betak, integrand, r)
implicit none
double precision, dimension(:,:) :: betak, integrand
double precision, dimension(:)   :: r

! Internal variables

double precision :: dr1, limit
integer :: i, j, lbr, ubr, lbt, ubt

lbr = lbound (betak,1)
ubr = ubound (betak,1)
lbt = lbound (betak,2)
ubt = ubound (betak,2)

do j = lbt, ubt
  betak(lbr,j) = 0.0d0
  do i = lbr+1, ubr
    dr1 = r(i) - r(i-1)
    betak(i,j) = betak(i-1,j) + dr1*0.5d0*(integrand(i-1,j) + integrand(i,j))
  end do
  limit = betak(ubr,j)
  do i = lbr, ubr
    betak(i,j) = betak(i,j) - limit
  end do
end do

end subroutine solvebetak

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine integration2d (integrand, outcome, r, t, storage)
implicit none
double precision, dimension(:,:) :: integrand
double precision :: outcome
double precision, dimension(:)   :: r, t, storage

! Internal variables

double precision :: dr1, dt1
integer :: i, j, lbr, ubr, lbt, ubt
double precision :: d0, d1, d2, f1, f2, f3

! auxiliary variables in the Kahan summation algorithm
double precision :: c, y, tt

lbr = lbound (integrand,1)
ubr = ubound (integrand,1)
lbt = lbound (integrand,2)
ubt = ubound (integrand,2)

do i = lbr, ubr
storage(i) = 0.0d0
c = 0.0d0
do j = lbt+1, ubt
    dt1 = t(j) - t(j-1)
    y = 0.5d0*dt1*(integrand(i,j) + integrand(i,j-1)) - c
    tt = storage(i) + y
    c = (tt - storage(i)) - y
    storage(i) = tt
end do
end do

outcome = 0.0d0
c = 0.0d0
do i = lbr+1, ubr
  dr1 = r(i) - r(i-1)
  y = 0.5d0*dr1*(storage(i) + storage(i-1)) - c
  tt = outcome + y
  c = (tt - outcome) - y
  outcome = tt
end do

end subroutine integration2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine integration1d (integrand, outcome, t)
implicit none
double precision, dimension(:) :: integrand
double precision :: outcome
double precision, dimension(:) :: t

! Internal variables

double precision :: d1, f1, f2, f3
integer :: j, lbt, ubt

! auxiliary variables in the Kahan summation algorithm
double precision :: c, y, tt

lbt = lbound (integrand,1)
ubt = ubound (integrand,1)

outcome = 0.0d0
c = 0.0d0

do j = lbt+1, ubt
    d1 = t(j) - t(j-1)
    y = 0.5d0*d1*(integrand(j) + integrand(j-1)) - c
    tt = outcome + y
    c = (tt - outcome) - y
    outcome = tt
end do

end subroutine integration1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine exact(q, sq, a, b, c, r, t, boundary, kmax, lmax, ab, IPIV, istart, &
  equ, axis, equator, horizon, betahor)
implicit none

double precision, dimension(:,:) :: q, sq, ab
double precision, dimension(:,:) :: a, b, c
double precision, dimension(:)   :: r, t, boundary
integer :: kmax, lmax
integer, dimension (:) :: IPIV
logical :: istart
character :: equ
double precision, dimension(:)   :: axis, equator, horizon
double precision, dimension(:,:) :: betahor

! Internal variables

double precision, dimension (:), allocatable :: rhs
integer :: INFO1, INFO2
integer :: jmax
integer :: k, l, jj
double precision :: dphi, dr, f
double precision :: d1, d2, d3, d4
integer :: kin
double precision :: dr1, dr2, dr3, dt1, dt2, dt3
integer :: version

version = 1

jmax = kmax*lmax

allocate (rhs (jmax))

if (istart) then

ab(:,:) = 0.0d0
 
rhs (:)   = 0.0d0

! the boundary condition at r = rin

do l = 1, lmax

if (equ .eq. 's') then

   k = 1
   jj = (k - 1)*lmax + l 

!  numery k, l

   ab(5*lmax + 1, jj) = 1.0d0

   rhs(jj) = 0.0d0

   k = 2
   jj = (k - 1)*lmax + l

   d1 = r(2) - r(1)
   d2 = r(3) - r(2)
   d3 = r(4) - r(3)
   d4 = r(5) - r(4)

! this condition can be also implemented differently

   select case (version)

     case (1)

! numery (jj,jj)
       ab(5*lmax + 1, jj) = ((d1 + d2)*(d1 + d2 + d3)*(d1 + d2 + d3 + d4))/ &
        (d1*d2*(d2 + d3)*(d2 + d3 + d4))
! numery (jj,jj+lmax)
       ab(4*lmax + 1, jj + lmax) = -((d1*(d1 + d2 + d3)*(d1 + d2 + d3 + d4))/ &
        (d2*(d1 + d2)*d3*(d3 + d4)))
! numery (jj,jj+2*lmax)
       ab(3*lmax + 1, jj + 2*lmax) = (d1*(d1 + d2)*(d1 + d2 + d3 + d4))/ &
        (d3*(d2 + d3)*(d1 + d2 + d3)*d4)
! numery (jj,jj+3*lmax)
       ab(2*lmax + 1, jj + 3*lmax) = -((d1*(d1 + d2)*(d1 + d2 + d3))/ &
        (d4*(d3 + d4)*(d2 + d3 + d4)*(d1 + d2 + d3 + d4)))

       rhs(jj) = 0.0d0

     case (2)

       ab(5*lmax + 1, jj) = 1.0d0  
       rhs(jj) = betahor(l,1)
 
   end select

   k = 3
   jj = (k - 1)*lmax + l

   d1 = r(2) - r(1)
   d2 = r(3) - r(2)
   d3 = r(4) - r(3)
   d4 = r(5) - r(4)

  select case (version)

    case (1)

! numery (jj,jj - lmax)
      ab(6*lmax + 1, jj - lmax) = (-2.0d0*(3.0d0*d1**2 + 3.0d0*d2**2 + d3*(d3 + d4) + 2.0d0*d2*(2.0d0*d3 + d4) + &
        &  2.0d0*d1*(3.0d0*d2 + 2.0d0*d3 + d4)))/(d1*d2*(d2 + d3)*(d2 + d3 + d4))
! numery (jj,jj)
      ab(5*lmax + 1, jj) = (2.0d0*(3.0d0*d1**2 + (d2 + d3)*(d2 + d3 + d4) + 2.0d0*d1*(2.0d0*d2 + 2.0d0*d3 + d4)))/ &
        &  (d2*(d1 + d2)*d3*(d3 + d4))
! numery (jj,jj + lmax)
      ab(4*lmax + 1, jj + lmax) = (-2.0d0*(3.0d0*d1**2 + d2*(d2 + d3 + d4) + 2.0d0*d1*(2.0d0*d2 + d3 + d4)))/ &
        &  (d3*(d2 + d3)*(d1 + d2 + d3)*d4)
! numery (jj,jj + 2*lmax)
      ab(3*lmax + 1, jj + 2*lmax) = (2.0d0*(3.0d0*d1**2 + d2*(d2 + d3) + 2.0d0*d1*(2*d2 + d3)))/ &
        &  (d4*(d3 + d4)*(d2 + d3 + d4)*(d1 + d2 + d3 + d4))

      rhs(jj) = 0.0d0

    case (2)

      ab(5*lmax + 1, jj) = 1.0d0
      rhs(jj) = betahor(l,2)

  end select

   k = 4
   jj = (k - 1)*lmax + l

   d1 = r(2) - r(1)
   d2 = r(3) - r(2)
   d3 = r(4) - r(3)
   d4 = r(5) - r(4)

   select case (version)

     case (1)

! numery (jj,jj - 2*lmax)
       ab(7*lmax + 1, jj - 2*lmax) = (6.0d0*(3.0d0*d1 + 3.0d0*d2 + 2.0d0*d3 + d4))/ &
         (d1*d2*(d2 + d3)*(d2 + d3 + d4))
! numery (jj,jj - lmax)
       ab(6*lmax + 1, jj - lmax) = (-6.0d0*(3.0d0*d1 + 2.0d0*d2 + 2*d3 + d4))/ &
         (d2*(d1 + d2)*d3*(d3 + d4))
! numery (jj,jj)
       ab(5*lmax + 1, jj) = (6.0d0*(3.0d0*d1 + 2.0d0*d2 + d3 + d4))/ &
         (d3*(d2 + d3)*(d1 + d2 + d3)*d4)
! numery (jj,jj + lmax)
       ab(4*lmax + 1, jj + lmax) = (-6.0d0*(3.0d0*d1 + 2.0d0*d2 + d3))/ &
         (d4*(d3 + d4)*(d2 + d3 + d4)*(d1 + d2 + d3 + d4))

       rhs(jj) = 0.0d0

     case (2)
   
       ab(5*lmax + 1, jj) = 1.0d0
       rhs(jj) = betahor(l,3)

   end select

else

   k = 1

   d1 = r(2) - r(1)
   d2 = r(3) - r(2)

   jj = (k - 1)*lmax + l 

   ab(5*lmax + 1, jj) = 1.0d0

   rhs(jj) = horizon(l)

 end if

 k = kmax
 jj = (k - 1)*lmax + l

! numbers k, l

 ab(5*lmax + 1, jj) = 1.0d0

 rhs(jj) = boundary(l)
end do


if (equ .eq. 's') then
  kin = 5
else
  kin = 2
end if

do k = kin, kmax -1

! the boundary condition ar t = 0

 l = 1
 jj = (k - 1)*lmax + l

! pochodna znika na osi

 ab(5*lmax + 1, jj) = 1.0d0

 rhs(jj) = axis(k)

! warunek na rowniku

 l = lmax
 jj = (k - 1)*lmax + l

 ab(5*lmax + 1, jj) = 1.0d0
 rhs(jj) = equator(k)

end do

! uzupelniamy wspolczynniki odpowiadajace glownym rownaniom

do l = 2, lmax - 1
  do k = kin, kmax - 1

  dr1 = r(k) - r(k-1)
  dr2 = r(k+1) - r(k)
  dr3 = dr1 + dr2

  dt1 = t(l) - t(l-1)
  dt2 = t(l+1) - t(l)
  dt3 = dt1 + dt2

  jj = (k - 1)*lmax + l

!  numery k-1,l

  ab(6*lmax + 1, jj - lmax) = (2.0d0 - dr2*a(k,l))/(dr1*dr3)

!  numery k,l-1

  ab(5*lmax + 2, jj - 1) =  (2.0d0*b(k,l) - dt2*c(k,l))/(dt1*dt3)

!  numery k,l

   ab(5*lmax + 1, jj) = - 2.0d0/(dr1*dr2) + a(k,l)/dr1 - a(k,l)/dr2 - (2.0d0*b(k,l))/(dt1*dt2) + c(k,l)/dt1 - c(k,l)/dt2

!  numery k, l+1

  ab(5*lmax,jj + 1) = (2.0d0*b(k,l) + dt1*c(k,l))/(dt2*dt3)

!  numery k+1,l

  ab(4*lmax + 1, jj + lmax) = (2.0d0 + dr1*a(k,l))/(dr2*dr3)

  rhs(jj) = sq(k,l)

  end do
end do

  call DGBTRF(jmax,jmax,2*lmax,3*lmax,ab,7*lmax+1,IPIV,INFO1)
  call DGBTRS('N',jmax,2*lmax,3*lmax,1,ab,7*lmax+1,IPIV,rhs,jmax,INFO2)

  if (INFO1 .ne. 0 .or. INFO2 .ne. 0) then
    write(*,*) 'Error in lapack DGETRF: ', INFO1, &
      ' DGETRS: ', INFO2
    stop
  end if

do k = 1, kmax
  do l = 1, lmax
    jj = (k - 1)*lmax + l
    q(k,l) = rhs(jj)
  end do
end do

else

rhs (:)   = 0.0d0

do l = 1, lmax
 k = 1
 jj = (k - 1)*lmax + l 

   rhs(jj) = horizon(l)

 if (equ .eq. 's') then
   k = 2
   jj = (k - 1)*lmax + l

   select case (version)

     case (1)

       rhs(jj) = 0.0d0

     case (2)

       rhs(jj) = betahor(l,1)

   end select

   k = 3
   jj = (k - 1)*lmax + l

   select case (version)

   case (1)

     rhs(jj) = 0.0d0

   case (2)

     rhs(jj) = betahor(l,2)

   end select

   k = 4
   jj = (k - 1)*lmax + l

   select case (version)

   case (1)

     rhs(jj) = 0.0d0

   case (2)

     rhs(jj) = betahor(l,3)

   end select

 end if

 k = kmax
 jj = (k - 1)*lmax + l

 rhs(jj) = boundary(l)
end do

if (equ .eq. 's') then
  kin = 5
else
  kin = 2
end if

do k = kin, kmax -1
 l = 1
 jj = (k - 1)*lmax + l
 rhs(jj) = axis(k)
 l = lmax
 jj = (k - 1)*lmax + l
 rhs(jj) = equator(k)
end do

do l = 2, lmax - 1
  do k = kin, kmax - 1

  jj = (k - 1)*lmax + l
  rhs(jj) = sq(k,l)

  end do
end do

  call DGBTRS('N',jmax,2*lmax,3*lmax,1,ab,7*lmax+1,IPIV,rhs,jmax,INFO2)

  if (INFO2 .ne. 0) then
    write(*,*) 'Error in lapack DGETRS: ', INFO2
    stop
  end if

do k = 1, kmax
  do l = 1, lmax
    jj = (k - 1)*lmax + l
    q(k,l) = rhs(jj)
  end do
end do

end if


deallocate (rhs)

end subroutine exact


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine derivr(f, dfdr, r)
implicit none
double precision, dimension(:,:) :: f, dfdr
double precision, dimension(:)   :: r

! Internal variables

double precision :: dr1, dr2, dr3, dr4
integer :: i, j, lbr, ubr, lbt, ubt

lbr = lbound (f,1)
ubr = ubound (f,1)
lbt = lbound (f,2)
ubt = ubound (f,2)

do j = lbt, ubt

  dfdr(lbr,j) = 0.0d0

  do i = lbr + 1, ubr - 1 

    dr1 = r(i) - r(i-1)
    dr2 = r(i+1) - r(i)
    dr3 = dr1 + dr2

    dfdr(i,j) = - f(i-1,j)*dr2/(dr1*dr3) &
       + f(i,j)*(dr2 - dr1)/(dr1*dr2) + f(i+1,j)*dr1/(dr2*dr3)

  end do

  dr1 = r(ubr - 2) - r(ubr - 3)
  dr2 = r(ubr - 1) - r(ubr - 2)
  dr3 = r(ubr) - r(ubr - 1)

  dfdr(ubr,j) = -((dr3*(dr2 + dr3))/(dr1*(dr1 + dr2)*(dr1 + dr2 + dr3)))*f(ubr-3,j) &
    + ((dr3*(dr1 + dr2 + dr3))/(dr1*dr2*(dr2 + dr3)))*f(ubr-2,j) &
    - (((dr2 + dr3)*(dr1 + dr2 + dr3))/(dr2*(dr1 + dr2)*dr3))*f(ubr-1,j) &
    + ((dr2**2 + 4.0d0*dr2*dr3 + 3.0d0*dr3**2 + dr1*(dr2 + 2.0d0*dr3))/ &
    (dr3*(dr2 + dr3)*(dr1 + dr2 + dr3)))*f(ubr,j)

end do

end subroutine derivr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine derivt(f, dfdt, t)
implicit none
double precision, dimension(:,:) :: f, dfdt
double precision, dimension(:)   :: t

! Internal variables

double precision :: dt1, dt2, dt3
integer :: i, j, lbr, ubr, lbt, ubt

lbr = lbound (f,1)
ubr = ubound (f,1)
lbt = lbound (f,2)
ubt = ubound (f,2)

do i = lbr, ubr

  dfdt(i,lbt) = 0.0d0

  do j = lbt + 1, ubt - 1

    dt1 = t(j) - t(j-1)
    dt2 = t(j+1) - t(j)
    dt3 = dt1 + dt2

    dfdt(i,j) = - f(i,j-1)*dt2/(dt1*dt3) &
       + f(i,j)*(dt2 - dt1)/(dt1*dt2) + f(i,j+1)*dt1/(dt2*dt3)

  end do

  dfdt(i,ubt) = 0.0d0

end do

end subroutine derivt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diskparams(c1, w, omega1, omega2, r1, r2, beta1, beta2, alpha1, alpha2, psi1, psi2, a)
implicit none
double precision :: c1, w, omega1, omega2, r1, r2, beta1, beta2, alpha1, alpha2, psi1, psi2, a

! Internal variables

double precision :: param
double precision, dimension (2,2) :: jak
double precision, dimension (2)   :: eqs, sol
integer, dimension (2) :: IPIV
integer :: INFO, it, nit

! Uwaga. Tutaj W = w^(4/3)

nit = 800

param = 1.0d-15

do it = 1, nit

jak(1,1) =         -(psi1**4*r1**2) + (2.0d0*(beta1 + omega1)*psi1**4*r1**2* &
     &     (a**2*alpha2**2*((1.0d0 - a*omega1)**0.3333333333333333d0*(-1.0d0 + 3.0d0*a*omega1)* &
     &           omega2**1.3333333333333333d0 + &
     &          omega1**1.3333333333333333d0*(1.0d0 - a*omega2)**0.3333333333333333d0 - &
     &          3.0d0*a*omega1**1.3333333333333333d0*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0 &
     &          ) + (beta2 + omega2)* &
     &        ((1.0d0 - 3.0d0*a*omega1)*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**0.3333333333333333d0 + &
     &          a**2*beta2*(1.0d0 - 3.0d0*a*omega1)*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**1.3333333333333333d0 - &
     &          a**2*beta2*omega1**1.3333333333333333d0* &
     &           (1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a**2*(2.0d0 + 3.0d0*a*beta2)*omega1**1.3333333333333333d0*omega2* &
     &           (1.0d0 - a*omega2)**0.3333333333333333d0)*psi2**4*r2**2))/ &
     &   (omega1**0.3333333333333333d0*(-1.0d0 + a*omega1)* &
     &     (alpha2**2*((1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a*(omega1 - 3.0d0*omega2)*(1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a**2*(3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &              omega2**1.3333333333333333d0 - &
     &             3.0d0*omega1*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0)) - &
     &       (beta2 + omega2)*((1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &          a*omega1*(1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &          3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**0.3333333333333333d0*(1.0d0 + a**2*beta2*omega2))*psi2**4*r2**2)) + &
     &  (a*(alpha1**2 - (beta1 + omega1)**2*psi1**4*r1**2)* &
     &     (a**2*alpha2**2*((1.0d0 - a*omega1)**0.3333333333333333d0*(-1.0d0 + 3.0d0*a*omega1)* &
     &           omega2**1.3333333333333333d0 + &
     &          omega1**1.3333333333333333d0*(1.0d0 - a*omega2)**0.3333333333333333d0 - &
     &          3.0d0*a*omega1**1.3333333333333333d0*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0 &
     &          ) + (beta2 + omega2)* &
     &        ((1.0d0 - 3.0d0*a*omega1)*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**0.3333333333333333d0 + &
     &          a**2*beta2*(1.0d0 - 3.0d0*a*omega1)*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**1.3333333333333333d0 - &
     &          a**2*beta2*omega1**1.3333333333333333d0* &
     &           (1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a**2*(2.0d0 + 3.0d0*a*beta2)*omega1**1.3333333333333333d0*omega2* &
     &           (1.0d0 - a*omega2)**0.3333333333333333d0)*psi2**4*r2**2))/ &
     &   (omega1**0.3333333333333333d0*(-1.0d0 + a*omega1)**2* &
     &     (alpha2**2*((1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a*(omega1 - 3.0d0*omega2)*(1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a**2*(3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &              omega2**1.3333333333333333d0 - &
     &             3.0d0*omega1*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0)) - &
     &       (beta2 + omega2)*((1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &          a*omega1*(1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &          3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**0.3333333333333333d0*(1.0d0 + a**2*beta2*omega2))*psi2**4*r2**2)) + &
     &  ((alpha1**2 - (beta1 + omega1)**2*psi1**4*r1**2)* &
     &     (a**2*alpha2**2*((1.0d0 - a*omega1)**0.3333333333333333d0*(-1.0d0 + 3.0d0*a*omega1)* &
     &           omega2**1.3333333333333333d0 + &
     &          omega1**1.3333333333333333d0*(1.0d0 - a*omega2)**0.3333333333333333d0 - &
     &          3.0d0*a*omega1**1.3333333333333333d0*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0 &
     &          ) + (beta2 + omega2)* &
     &        ((1.0d0 - 3.0d0*a*omega1)*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**0.3333333333333333d0 + &
     &          a**2*beta2*(1.0d0 - 3.0d0*a*omega1)*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**1.3333333333333333d0 - &
     &          a**2*beta2*omega1**1.3333333333333333d0* &
     &           (1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a**2*(2.0d0 + 3.0d0*a*beta2)*omega1**1.3333333333333333d0*omega2* &
     &           (1.0d0 - a*omega2)**0.3333333333333333d0)*psi2**4*r2**2))/ &
     &   (3.0d0*omega1**1.3333333333333333d0*(-1.0d0 + a*omega1)* &
     &     (alpha2**2*((1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a*(omega1 - 3.0d0*omega2)*(1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a**2*(3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &              omega2**1.3333333333333333d0 - &
     &             3.0d0*omega1*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0)) - &
     &       (beta2 + omega2)*((1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &          a*omega1*(1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &          3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**0.3333333333333333d0*(1.0d0 + a**2*beta2*omega2))*psi2**4*r2**2)) + &
     &  (2.0d0*a*(alpha1**2 - (beta1 + omega1)**2*psi1**4*r1**2)* &
     &     (a*alpha2**2*(5.0d0*a*omega2**1.3333333333333333d0 - &
     &          6.0d0*a**2*omega1*omega2**1.3333333333333333d0 + &
     &          2.0d0*omega1**0.3333333333333333d0*(1.0d0 - a*omega1)**0.6666666666666666d0* &
     &           (1.0d0 - a*omega2)**0.3333333333333333d0 - &
     &          6.0d0*a*omega1**0.3333333333333333d0*(1.0d0 - a*omega1)**0.6666666666666666d0* &
     &           omega2*(1.0d0 - a*omega2)**0.3333333333333333d0) + &
     &       (beta2 + omega2)*((-5.0d0 + 6.0d0*a*omega1)*omega2**0.3333333333333333d0 + &
     &          a**2*beta2*(-5.0d0 + 6.0d0*a*omega1)*omega2**1.3333333333333333d0 - &
     &          2.0d0*a*beta2*omega1**0.3333333333333333d0*(1.0d0 - a*omega1)**0.6666666666666666d0* &
     &           (1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          2.0d0*a*(2.0d0 + 3.0d0*a*beta2)*omega1**0.3333333333333333d0* &
     &           (1.0d0 - a*omega1)**0.6666666666666666d0*omega2* &
     &           (1.0d0 - a*omega2)**0.3333333333333333d0)*psi2**4*r2**2))/ &
     &   (3.0d0*omega1**0.3333333333333333d0*(1.0d0 - a*omega1)**1.6666666666666667d0* &
     &     (alpha2**2*((1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a*(omega1 - 3.0d0*omega2)*(1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a**2*(3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &              omega2**1.3333333333333333d0 - &
     &             3.0d0*omega1*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0)) - &
     &       (beta2 + omega2)*((1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &          a*omega1*(1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &          3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**0.3333333333333333d0*(1.0d0 + a**2*beta2*omega2))*psi2**4*r2**2)) + &
     &  ((alpha1**2 - (beta1 + omega1)**2*psi1**4*r1**2)* &
     &     (a**2*alpha2**2*((1.0d0 - a*omega1)**0.3333333333333333d0*(-1.0d0 + 3.0d0*a*omega1)* &
     &           omega2**1.3333333333333333d0 + &
     &          omega1**1.3333333333333333d0*(1.0d0 - a*omega2)**0.3333333333333333d0 - &
     &          3.0d0*a*omega1**1.3333333333333333d0*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0 &
     &          ) + (beta2 + omega2)* &
     &        ((1.0d0 - 3.0d0*a*omega1)*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**0.3333333333333333d0 + &
     &          a**2*beta2*(1.0d0 - 3.0d0*a*omega1)*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**1.3333333333333333d0 - &
     &          a**2*beta2*omega1**1.3333333333333333d0* &
     &           (1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a**2*(2.0d0 + 3.0d0*a*beta2)*omega1**1.3333333333333333d0*omega2* &
     &           (1.0d0 - a*omega2)**0.3333333333333333d0)*psi2**4*r2**2)* &
     &     (a*alpha2**2*(-((a*(-2.0d0 + 3.0d0*a*omega1)*omega2**1.3333333333333333d0)/ &
     &             (omega1**0.3333333333333333d0*(1.0d0 - a*omega1)**0.6666666666666666d0)) + &
     &          (1.0d0 - a*omega2)**0.3333333333333333d0 - &
     &          3.0d0*a*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0) - &
     &       (beta2 + omega2)*(a*(1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) - &
     &          (a*omega1**0.6666666666666666d0*omega2**0.3333333333333333d0* &
     &             (1.0d0 + a**2*beta2*omega2))/(1.0d0 - a*omega1)**0.6666666666666666d0 + &
     &          (2.0d0*(1.0d0 - a*omega1)**0.3333333333333333d0*omega2**0.3333333333333333d0* &
     &             (1.0d0 + a**2*beta2*omega2))/omega1**0.3333333333333333d0)*psi2**4*r2**2))/ &
     &   (omega1**0.3333333333333333d0*(-1.0d0 + a*omega1)* &
     &     (-(alpha2**2*((1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &             a*(omega1 - 3.0d0*omega2)*(1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &             a**2*(3.0d0*omega1**0.6666666666666666d0* &
     &                 (1.0d0 - a*omega1)**0.3333333333333333d0*omega2**1.3333333333333333d0 - &
     &                3.0d0*omega1*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0))) + &
     &        (beta2 + omega2)*((1.0d0 - a*omega2)**0.3333333333333333d0* &
     &            (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &           a*omega1*(1.0d0 - a*omega2)**0.3333333333333333d0* &
     &            (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &           3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &            omega2**0.3333333333333333d0*(1.0d0 + a**2*beta2*omega2))*psi2**4*r2**2)**2)


jak(1,2) =  -((-1.0d0 + 2.0d0*a*omega1)*(alpha1**2 - (beta1 + omega1)**2*psi1**4*r1**2)* &
     &     (alpha2**2 - (beta2 + omega2)**2*psi2**4*r2**2)* &
     &     (2.0d0*a**2*alpha2**2*omega2*(-2.0d0 + 3.0d0*a*omega2) + &
     &       (2.0d0*omega2*(2.0d0 - 3.0d0*a*omega2) + 2.0d0*a**2*beta2**2*omega2*(2.0d0 - 3.0d0*a*omega2) + &
     &          beta2*(1.0d0 + 6.0d0*a*omega2 - 11.0d0*a**2*omega2**2))*psi2**4*r2**2))/ &
     &  (3.0d0*omega1**0.3333333333333333d0*(1.0d0 - a*omega1)**0.6666666666666666d0* &
     &    omega2**0.6666666666666666d0*(1.0d0 - a*omega2)**0.6666666666666666d0* &
     &    (-(alpha2**2*((1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &            a*(omega1 - 3.0d0*omega2)*(1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &            a**2*(3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &                omega2**1.3333333333333333d0 - &
     &               3.0d0*omega1*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0))) + &
     &       (beta2 + omega2)*((1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &          a*omega1*(1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &          3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**0.3333333333333333d0*(1.0d0 + a**2*beta2*omega2))*psi2**4*r2**2)**2) 

jak(2,1) = -(((beta1 + omega1)*psi1**4*r1**2)/ &
     &     sqrt(alpha1**2 - (beta1 + omega1)**2*psi1**4*r1**2)) - &
     &  (sqrt(((3.0d0 - 9.0d0*a*omega2 + 6.0d0*a**2*omega2**2)* &
     &         (-alpha2**2 + (beta2 + omega2)**2*psi2**4*r2**2))/ &
     &       (alpha2**2*(-1.0d0 + 3.0d0*a*omega2) - &
     &         (beta2 + omega2)*(2.0d0*omega2 + beta2*(-1.0d0 + 3.0d0*a*omega2))*psi2**4*r2**2))* &
     &     (-2.0d0*a*omega1*(1.0d0 - a*omega1)**0.3333333333333333d0*omega2**0.3333333333333333d0* &
     &        (3.0d0*a**2*alpha2**2*omega2 - &
     &          3.0d0*(beta2 + omega2)*(1.0d0 + a**2*beta2*omega2)*psi2**4*r2**2) + &
     &       (1.0d0 - a*omega1)**1.3333333333333333d0*omega2**0.3333333333333333d0* &
     &        (3.0d0*a**2*alpha2**2*omega2 - &
     &          3.0d0*(beta2 + omega2)*(1.0d0 + a**2*beta2*omega2)*psi2**4*r2**2) + &
     &       3.0d0*a**2*omega1**1.3333333333333333d0*(1.0d0 - a*omega2)**0.3333333333333333d0* &
     &        (alpha2**2*(-1.0d0 + 3.0d0*a*omega2) - &
     &          (beta2 + omega2)*(2.0d0*omega2 + beta2*(-1.0d0 + 3.0d0*a*omega2))*psi2**4*r2**2)))/ &
     &   (omega1**0.3333333333333333d0*(1.0d0 - a*omega2)**0.3333333333333333d0* &
     &     (3.0d0 - 9.0d0*a*omega2 + 6.0d0*a**2*omega2**2)* &
     &     sqrt(alpha2**2 - (beta2 + omega2)**2*psi2**4*r2**2)* &
     &     sqrt(3.0d0 - 3.0d0*a**2*omega1**2 - &
     &       (3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**1.3333333333333333d0* &
     &          omega2**0.3333333333333333d0* &
     &          (3.0d0*a**2*alpha2**2*omega2 - &
     &            3.0d0*(beta2 + omega2)*(1.0d0 + a**2*beta2*omega2)*psi2**4*r2**2))/ &
     &        ((1.0d0 - a*omega2)**0.3333333333333333d0* &
     &          (alpha2**2*(-1.0d0 + 3.0d0*a*omega2) - &
     &            (beta2 + omega2)*(2.0d0*omega2 + beta2*(-1.0d0 + 3.0d0*a*omega2))*psi2**4*r2**2))))


jak(2,2) = ((-1.0d0 + a*omega1)*(omega2**0.6666666666666666d0*(1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &      a*omega1*omega2**0.6666666666666666d0*(1.0d0 - a*omega2)**0.3333333333333333d0 - &
     &      omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0*(1.0d0 + a*omega2))* &
     &    sqrt(-(((1.0d0 - 3.0d0*a*omega2 + 2.0d0*a**2*omega2**2)* &
     &          (-alpha2**2 + (beta2 + omega2)**2*psi2**4*r2**2))/ &
     &        (alpha2**2*(1.0d0 - 3.0d0*a*omega2) + &
     &          (beta2 + omega2)*(2.0d0*omega2 + beta2*(-1.0d0 + 3.0d0*a*omega2))*psi2**4*r2**2)))* &
     &    (2.0d0*a**2*alpha2**2*omega2*(-2.0d0 + 3.0d0*a*omega2) + &
     &      (2.0d0*omega2*(2.0d0 - 3.0d0*a*omega2) + 2.0d0*a**2*beta2**2*omega2*(2.0d0 - 3.0d0*a*omega2) + &
     &         beta2*(1.0d0 + 6.0d0*a*omega2 - 11.0d0*a**2*omega2**2))*psi2**4*r2**2))/ &
     &  (2.0d0*omega2**0.6666666666666666d0*(1.0d0 - 2.0d0*a*omega2)**2* &
     &    (1.0d0 - a*omega2)**2.3333333333333335d0* &
     &    sqrt(alpha2**2 - (beta2 + omega2)**2*psi2**4*r2**2)* &
     &    sqrt(((-1.0d0 + a*omega1)*(alpha2**2* &
     &           (1.0d0 + a*(omega1 - 4.0d0*omega2) + 3.0d0*a**3*omega1*omega2**2 + &
     &             a**2*omega2*(-4.0d0*omega1 + 3.0d0*omega2 + &
     &                3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &                 omega2**0.3333333333333333d0*(1.0d0 - a*omega2)**0.6666666666666666d0)) - &
     &          (beta2 + omega2)*(-2.0d0*(1.0d0 + a*omega1)*omega2 + &
     &             2.0d0*a*(1.0d0 + a*omega1)*omega2**2 + &
     &             3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &              omega2**0.3333333333333333d0*(1.0d0 - a*omega2)**0.6666666666666666d0 + &
     &             beta2*(1.0d0 + a*(omega1 - 4.0d0*omega2) + 3.0d0*a**3*omega1*omega2**2 + &
     &                a**2*omega2*(-4.0d0*omega1 + 3.0d0*omega2 + &
     &                   3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &                    omega2**0.3333333333333333d0*(1.0d0 - a*omega2)**0.6666666666666666d0)))* &
     &           psi2**4*r2**2))/ &
     &      ((-1.0d0 + a*omega2)*(alpha2**2*(1.0d0 - 3.0d0*a*omega2) + &
     &          (beta2 + omega2)*(2.0d0*omega2 + beta2*(-1.0d0 + 3.0d0*a*omega2))*psi2**4*r2**2))))


eqs(1) =  -(beta1*psi1**4*r1**2) - omega1*psi1**4*r1**2 - &
     &  ((alpha1**2 - (beta1 + omega1)**2*psi1**4*r1**2)* &
     &     (a**2*alpha2**2*((1.0d0 - a*omega1)**0.3333333333333333d0*(-1.0d0 + 3.0d0*a*omega1)* &
     &           omega2**1.3333333333333333d0 + &
     &          omega1**1.3333333333333333d0*(1.0d0 - a*omega2)**0.3333333333333333d0 - &
     &          3.0d0*a*omega1**1.3333333333333333d0*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0) &
     &        + (beta2 + omega2)*((1.0d0 - 3.0d0*a*omega1)*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**0.3333333333333333d0 + &
     &          a**2*beta2*(1.0d0 - 3.0d0*a*omega1)*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**1.3333333333333333d0 - &
     &          a**2*beta2*omega1**1.3333333333333333d0*(1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a**2*(2.0d0 + 3.0d0*a*beta2)*omega1**1.3333333333333333d0*omega2* &
     &           (1.0d0 - a*omega2)**0.3333333333333333d0)*psi2**4*r2**2))/ &
     &   (omega1**0.3333333333333333d0*(-1.0d0 + a*omega1)* &
     &     (alpha2**2*((1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a*(omega1 - 3.0d0*omega2)*(1.0d0 - a*omega2)**0.3333333333333333d0 + &
     &          a**2*(3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &              omega2**1.3333333333333333d0 - &
     &             3.0d0*omega1*omega2*(1.0d0 - a*omega2)**0.3333333333333333d0)) - &
     &       (beta2 + omega2)*((1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &          a*omega1*(1.0d0 - a*omega2)**0.3333333333333333d0* &
     &           (beta2 - 2.0d0*omega2 - 3.0d0*a*beta2*omega2) + &
     &          3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**0.3333333333333333d0* &
     &           omega2**0.3333333333333333d0*(1.0d0 + a**2*beta2*omega2))*psi2**4*r2**2))


eqs(2) = sqrt(alpha1**2 - (beta1 + omega1)**2*psi1**4*r1**2) - &
     &  (sqrt(alpha2**2 - (beta2 + omega2)**2*psi2**4*r2**2)* &
     &     sqrt(3.0d0 - 3.0d0*a**2*omega1**2 - &
     &       (3.0d0*omega1**0.6666666666666666d0*(1.0d0 - a*omega1)**1.3333333333333333d0* &
     &          omega2**0.3333333333333333d0* &
     &          (3.0d0*a**2*alpha2**2*omega2 - &
     &            3.0d0*(beta2 + omega2)*(1.0d0 + a**2*beta2*omega2)*psi2**4*r2**2))/ &
     &        ((1.0d0 - a*omega2)**0.3333333333333333d0* &
     &          (alpha2**2*(-1.0d0 + 3.0d0*a*omega2) - &
     &            (beta2 + omega2)*(2.0d0*omega2 + beta2*(-1.0d0 + 3.0d0*a*omega2))*psi2**4*r2**2))))/ &
     &   sqrt(((3.0d0 - 9.0d0*a*omega2 + 6.0d0*a**2*omega2**2)* &
     &       (-alpha2**2 + (beta2 + omega2)**2*psi2**4*r2**2))/ &
     &     (alpha2**2*(-1.0d0 + 3.0d0*a*omega2) - &
     &       (beta2 + omega2)*(2*omega2 + beta2*(-1.0d0 + 3.0d0*a*omega2))*psi2**4*r2**2))


sol(:) = - eqs(:)
call DGESV(2, 1, jak, 2, IPIV, sol, 2, INFO)

omega1 = omega1 + sol(1)
omega2 = omega2 + sol(2)

write(*,*) omega1, omega2

if ((abs(sol(1)) .lt. param) .and. (abs(sol(2)) .lt. param)) then
     exit
end if

end do

if (it .ge. nit) then
     write(*,*) 'diskparams: BRAK ZBIEZNOSCI'
     stop
end if


w = (omega2**0.3333333333333333d0*((beta2 + omega2)*psi2**4*r2**2 + &
     &      a**2*omega2*(-alpha2**2 + beta2*(beta2 + omega2)*psi2**4*r2**2)))/ &
     &  ((1.0d0 - a*omega2)**0.3333333333333333d0* &
     &    (alpha2**2*(1.0d0 - 3.0d0*a*omega2) + &
     &      (beta2 + omega2)*(2.0d0*omega2 + beta2*(-1.0d0 + 3.0d0*a*omega2))*psi2**4*r2**2))

c1 = sqrt(alpha2**2 - (beta2 + omega2)**2*psi2**4*r2**2)/ &
     &  sqrt(-((-1.0d0 + a*omega2)*(1.0d0 + a*omega2 - &
     &        3.0d0*omega2**0.6666666666666666d0*(1.0d0 - a*omega2)**0.3333333333333333d0*W)))

write(*,*) 'w:', w, 'C1:', c1 

end subroutine diskparams

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diskomega(omega, r, t, alpha, beta, psi, w, a)
implicit none
double precision :: omega, r, t, alpha, beta, psi, w, a

! Internal variables

integer :: it, nit
double precision :: popr, param

nit = 800

param = 1.0d-15

do it = 1, nit

popr = (-(alpha**2*(a**2*omega**1.3333333333333333d0 + &
     &         (1.0d0 - a*omega)**0.3333333333333333d0*W - &
     &         3.0d0*a*omega*(1.0d0 - a*omega)**0.3333333333333333d0*W)) + &
     &    (beta + omega)*psi**4*r**2* &
     &     (omega**0.3333333333333333d0 + a**2*beta*omega**1.3333333333333333d0 + &
     &       beta*(1.0d0 - a*omega)**0.3333333333333333d0*W - &
     &       (2.0d0 + 3.0d0*a*beta)*omega*(1.0d0 - a*omega)**0.3333333333333333d0*W)*sin(t)**2)/ &
     &  (omega**0.3333333333333333d0*(-1.0d0 + a*omega)* &
     &    (1.0d0 + a*omega - 3.0d0*omega**0.6666666666666666d0*(1.0d0 - a*omega)**0.3333333333333333d0* &
     &       W)*(-(psi**4*r**2*sin(t)**2) + &
     &      (2.0d0*(beta + omega)*psi**4*r**2* &
     &         (a**2*omega**1.3333333333333333d0 + (1.0d0 - a*omega)**0.3333333333333333d0*W - &
     &           3.0d0*a*omega*(1.0d0 - a*omega)**0.3333333333333333d0*W)*sin(t)**2)/ &
     &       (omega**0.3333333333333333d0* &
     &         (-1.0d0 + a**2*omega**2 + &
     &           3.0d0*omega**0.6666666666666666d0*(1.0d0 - a*omega)**1.3333333333333333d0*W)) - &
     &      ((a**2*omega**1.3333333333333333d0 + (1.0d0 - a*omega)**0.3333333333333333d0*W - &
     &           3.0d0*a*omega*(1.0d0 - a*omega)**0.3333333333333333d0*W)* &
     &         (-alpha**2 + (beta + omega)**2*psi**4*r**2*sin(t)**2))/ &
     &       (3.0d0*omega**1.3333333333333333d0*(-1.0d0 + a*omega)* &
     &         (1.0d0 + a*omega - 3.0d0*omega**0.6666666666666666d0* &
     &            (1.0d0 - a*omega)**0.3333333333333333d0*W)) - &
     &      (2.0d0*(a**2*omega**1.3333333333333333d0 + (1.0d0 - a*omega)**0.3333333333333333d0*W - &
     &            3.0d0*a*omega*(1.0d0 - a*omega)**0.3333333333333333d0*W)**2* &
     &         (-alpha**2 + (beta + omega)**2*psi**4*r**2*sin(t)**2))/ &
     &       (omega**0.6666666666666666d0*(-1.0d0 + a*omega)**2* &
     &         (1.0d0 + a*omega - 3.0d0*omega**0.6666666666666666d0* &
     &             (1.0d0 - a*omega)**0.3333333333333333d0*W)**2) + &
     &      (2.0d0*a*(2.0d0*a*omega**0.3333333333333333d0*(1.0d0 - a*omega)**0.6666666666666666d0 + &
     &           (-5.0d0 + 6.0d0*a*omega)*W)* &
     &         (-alpha**2 + (beta + omega)**2*psi**4*r**2*sin(t)**2))/ &
     &       (3.0d0*omega**0.3333333333333333d0*(1.0d0 - a*omega)**1.6666666666666667d0* &
     &         (-1.0d0 - a*omega + 3.0d0*omega**0.6666666666666666d0* &
     &            (1.0d0 - a*omega)**0.3333333333333333d0*W))))

omega = omega - popr

if (abs(popr) .lt. param) then
     exit
end if

end do

if (it .ge. nit) then
     write(*,*) 'diskomega: BRAK ZBIEZNOSCI'
     stop
end if

end subroutine diskomega

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine magbernoulli(h, alpha, psi, r, t, gam, kk, nn, c2)
implicit none
double precision :: h, alpha, psi, r, t, gam, kk, nn, c2

! Internal variables

integer :: it, nit
double precision :: popr, param, f

nit = 800
param = 1.0d-15

f = h

do it = 1, nit

popr =   (gam*(((-1.0d0 + gam)*(-1.0d0 + h))/(gam*kk))**(1.0d0 + 1.0d0/(1.0d0 - gam))*kk* &
     &    (1.0d0 + alpha**2*c2*h*(((-1.0d0 + gam)*(-1.0d0 + h))/(gam*kk))**(1.0d0/(-1.0d0 + gam))* &
     &        psi**4*r**2*Sin(t)**2)**(1.0d0 - nn)* &
     &    (-f + h*(1.0d0 + alpha**2*c2*h* &
     &           (((-1.0d0 + gam)*(-1.0d0 + h))/(gam*kk))**(1.0d0/(-1.0d0 + gam))*psi**4*r**2* &
     &           Sin(t)**2)**nn))/ &
     &  (gam*(((-1.0d0 + gam)*(-1.0d0 + h))/(gam*kk))**(1.0d0 + 1.0d0/(1.0d0 - gam))*kk + &
     &    alpha**2*c2*h*(-((-1.0d0 + gam)*(1.0d0 + nn)) + h*(-1.0d0 + gam + gam*nn))*psi**4* &
     &     r**2*Sin(t)**2)

h = h - popr

if (abs(popr) .lt. param) then
     exit
end if

end do

if (it .ge. nit) then
     write(*,*) 'magbernoulli: BRAK ZBIEZNOSCI'
     stop
end if

end subroutine magbernoulli

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine derivr1(f, dfdr, r)
implicit none
double precision, dimension(:)   :: f, dfdr
double precision, dimension(:)   :: r

! Internal variables

double precision :: dr1, dr2, dr3, dr4
integer :: i, lbr, ubr

lbr = lbound (f,1)
ubr = ubound (f,1)

  dfdr(lbr) = 0.0d0

  do i = lbr + 1, ubr - 1

    dr1 = r(i) - r(i-1)
    dr2 = r(i+1) - r(i)
    dr3 = dr1 + dr2

    dfdr(i) = - f(i-1)*dr2/(dr1*dr3) &
       + f(i)*(dr2 - dr1)/(dr1*dr2) + f(i+1)*dr1/(dr2*dr3)

  end do

  dr1 = r(ubr - 2) - r(ubr - 3)
  dr2 = r(ubr - 1) - r(ubr - 2)
  dr3 = r(ubr) - r(ubr - 1)

  dfdr(ubr) = -((dr3*(dr2 + dr3))/(dr1*(dr1 + dr2)*(dr1 + dr2 + dr3)))*f(ubr-3) &
    + ((dr3*(dr1 + dr2 + dr3))/(dr1*dr2*(dr2 + dr3)))*f(ubr-2) &
    - (((dr2 + dr3)*(dr1 + dr2 + dr3))/(dr2*(dr1 + dr2)*dr3))*f(ubr-1) &
    + ((dr2**2 + 4.0d0*dr2*dr3 + 3.0d0*dr3**2 + dr1*(dr2 + 2.0d0*dr3))/ &
    (dr3*(dr2 + dr3)*(dr1 + dr2 + dr3)))*f(ubr)


end subroutine derivr1
