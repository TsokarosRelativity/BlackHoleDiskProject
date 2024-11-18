module graysmod
    use iso_fortran_env, only: real64, int64
  contains
    subroutine cartesian_to_spherical(cart_coords, sph_coords, n_points)
        implicit none
        ! Input/output variables
        integer, intent(in) :: n_points
        real(real64), intent(in) :: cart_coords(n_points, 3)  ! Changed from real(int64)
        real(real64), intent(out) :: sph_coords(n_points, 3)  ! Changed from real(int64)
        
        ! Local variables
        integer :: i
        real(real64) :: x, y, z, r  ! Changed from real(int64)
        real(real64), parameter :: eps = 1.0d-14  ! Changed from real(int64)
        
        do i = 1, n_points
            x = cart_coords(i, 1)
            y = cart_coords(i, 2)
            z = cart_coords(i, 3)
            
            ! Calculate radius (r)
            r = sqrt(x**2 + y**2 + z**2)
            sph_coords(i, 1) = r
            
            ! Calculate theta (polar angle)
            if (r /= 0.0d0) then
                ! Clip z/r to [-1,1] range to avoid numerical issues
                sph_coords(i, 2) = acos(max(min(z/r, 1.0d0), -1.0d0))
            else
                sph_coords(i, 2) = 0.0d0
            end if
            
            ! Calculate phi (azimuthal angle)
            sph_coords(i, 3) = atan2(y, x)
            
            ! Zero out very small values (less than eps)
            if (abs(sph_coords(i, 1)) < eps) sph_coords(i, 1) = 0.0d0
            if (abs(sph_coords(i, 2)) < eps) sph_coords(i, 2) = 0.0d0
            if (abs(sph_coords(i, 3)) < eps) sph_coords(i, 3) = 0.0d0
        end do
    end subroutine cartesian_to_spherical
end module graysmod

program main
  use graysmod
  implicit none
  real(real64) :: cart_coords(16, 3), spherical_coords(16, 3)
  integer :: i, j, k, ios, n
  
  ! Initialize cart_coords to zero to avoid undefined values
  cart_coords = 0.0d0
  
  print *, "i:"
  do i = 1, 16, 4
    print *, i
  end do
  print *, "j:"
  do j = 4, 7
    print *, j
  end do 
  print *, "k:"
  do k = 5, 8
    print *, k
  end do
  
  print *, "cartesian grids:"
  n = 1
  ! Add bounds checking to prevent segmentation fault
  do i = 1, 16, 5
    do j = 4, 7
      do k = 5, 8
        if (n <= 16) then  ! Add bounds check
            cart_coords(n, 1) = real(i, real64)
            cart_coords(n, 2) = real(j, real64)
            cart_coords(n, 3) = real(k, real64)
            n = n + 1
        end if
      end do
    end do
  end do
  
  call cartesian_to_spherical(cart_coords, spherical_coords, 16)
  
  ! Write coordinates with specific format for better readability
  open(unit=10, file='spherical_coords.txt', status='unknown', action='write', iostat=ios)
  if (ios /= 0) then
    print *, "Error opening file"
    stop
  end if
  do i = 1, 16
    write(10, '(3ES15.6)') spherical_coords(i, 1), spherical_coords(i, 2), spherical_coords(i, 3)
  end do
  close(10)
  
  open(unit=22, file='cartesian_coords.txt', status='unknown', action='write', iostat=ios)
  if (ios /= 0) then
    print *, "Error opening file"
    stop
  end if
  do i = 1, 16
    write(22, '(3ES15.6)') cart_coords(i, 1), cart_coords(i, 2), cart_coords(i, 3)
  end do
  close(22)
end program main
