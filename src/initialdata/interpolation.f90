module interpolation_module
  use iso_fortran_env, only: real64, int64
  use m_KdTree
  implicit none

  ! Global variables
  real(real64), allocatable :: rad_arr(:), theta_arr(:), phi_arr(:)
  real(real64), allocatable :: data_3d(:,:)
  type(KdTree) :: tree_grid
  type(KdTree) :: tree_rad
  type(KdTree) :: tree_theta
  type(KdTree) :: tree_phi
  type :: gridparams 
     real(real64) :: xmin, ymin, zmin, dz, dy, dx
     integer(int64) :: nz, ny, nx
     character(len=40) :: outputfile
     real(real64), allocatable :: gridcoords(:,:), sphericalgrid(:,:)
  end type gridparams 

contains

  
  subroutine load_data()
    character(len=100) :: filename
    integer :: io_status, ind, rows
    ! Implement data loading here
    print *, "Loading in data..."
    filename = "df3d.txt"
    open(unit=8, file=filename, status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
      print *, "Error opening file: ", filename
      exit
    endif 
    
    rows = 64600901 !matrix dimensions -> (64600901, 27) 
    allocate(data_3d(rows, 27), stat=io_status)
    if (io_status /= 0) then 
      print *, "Error allocating data"
      exit
    end if
    
    do ind = 1, row
      read(8, '(27ES16.6)', stat=io_status) data_3d(ind, :)
      if (io_status /= 0) then
        print *, "Error reading data on line " , ind
        exit
      end if
    end do

    print *, "Data loaded successfully!"
  end subroutine load_data

  subroutine create_kdtree(kdtg, kdtr, kdtt, kdtp)
    type(KdTree), intent(out) :: kdtg
    type(KdTree), intent(out) :: kdtr
    type(KdTree), intent(out) :: kdtt
    type(KdTree), intent(out) :: kdtp
    print *, "Creating KD-tree..."
    kdtg = KdTree(rad_arr, theta_arr, phi_arr)
    kdtr = KdTree(rad_arr)
    kdtt = KdTree(theta_arr)
    kdtp = KdTree(phi_arr)
    print *, "KD-tree created successfully!"
  end subroutine create_kdtree_grid

  function get_nearest_point_index(tree, r, theta, phi) result(ind)
    real(real64), intent(in) :: r, theta, phi
    type(KdTree), intent(in) :: tree
    integer :: ind
    ind = search%nearest(tree, rad_arr, theta_arr, phi_arr, r, theta, phi) 
  end function get_nearest_point_index

  function interpolationpoint(point) result(vals)
    real(real64), intent(in) :: point(3)
    real(real64) :: vals(24)
    real(real64) :: rp, thetap, phip
    integer :: rid, thetad, phid
    logical :: rmatch, thetamatch, phimatch
    
    rp = point(1)
    thetap = point(2)
    phip = point(3)

    ! TODO: Implement searchsorted for rid, thetad, phid
    ! TODO: Implement isclose for rmatch, thetamatch, phimatch

    ! Implement the different cases for interpolation
    if (rmatch .and. thetamatch .and. phimatch) then
      ! Case 1: Exact match
      ! TODO: Implement
    else if (rmatch .and. thetamatch .and. .not. phimatch) then
      ! Case 2: Interpolate along phi
      ! TODO: Implement
    else if (rmatch .and. .not. thetamatch .and. phimatch) then
      ! Case 3: Interpolate along theta
      ! TODO: Implement
    else if (.not. rmatch .and. thetamatch .and. phimatch) then
      ! Case 4: Interpolate along r
      ! TODO: Implement
    else if (rmatch .and. .not. thetamatch .and. .not. phimatch) then
      ! Case 5: Interpolate along theta and phi
      ! TODO: Implement
    else if (.not. rmatch .and. thetamatch .and. .not. phimatch) then
      ! Case 6: Interpolate along r and phi
      ! TODO: Implement
    else if (.not. rmatch .and. .not. thetamatch .and. phimatch) then
      ! Case 7: Interpolate along r and theta
      ! TODO: Implement
    else
      ! Case 8: Interpolate along all dimensions
      ! TODO: Implement
    end if
  end function interpolationpoint

  subroutine process_line(line, grids)
    integer(int64) :: i, j, k, l, io_status, n_parts, ntot, ind
    character(len=250), intent(in) :: line
    character(len=250) :: filename, errmsg, tmp, prefix
    character(len=250), dimension(11) :: parts
    real(real64) :: xmax, ymax, zmax
    
    type(gridparams), dimension(1280), intent(out) :: grids
    
    filename = "grids_bh_disk_patrik"
    prefix = "CTS_bin-proc"
    
    open(unit=9, file=filename, status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
        print *, "Error opening file: ", filename
        stop
    end if
    
    do i = 1, 1280 
        read(9, '(A)', iostat=io_status) line
        if (io_status /= 0) then
            print *, "error reading data on line ", i
            exit
        end if 
        
        ! Split the line into parts
        tmp = ""
        n_parts = 1
        parts = ""
        do j = 1, len(trim(line))
            if (line(j:j) == ' ') then
                if (len(trim(tmp)) > 0) then
                    parts(n_parts) = tmp
                    n_parts = n_parts + 1
                    tmp = ""
                end if
            else
                tmp = trim(tmp) // line(j:j)
            end if
        end do
        parts(n_parts) = trim(tmp)
        
        ! Read the values with error checking
        read(parts(2), *, iostat=io_status) grids(i)%xmin
        if (io_status /= 0) then
            print *, "Error reading xmin on line ", i
            exit
        end if
        read(parts(3), *, iostat=io_status) grids(i)%ymin

        if (io_status /= 0) then
            print *, "Error reading ymin on line", i
            exit
        end if
        read(parts(4), *, iostat=io_status) grids(i)%zmin
        
        if (io_status /= 0) then
            print *, "Error reading zmin on line", i
            exit
        end if
        read(parts(5), *, iostat=io_status) grids(i)%dx
        
        if (io_status /= 0) then
            print *, "Error reading dx on line", i
            exit
        end if
        read(parts(6), *, iostat=io_status) grids(i)%dy
        
        if (io_status /= 0) then
            print *, "Error reading dy on line", i
            exit
        end if
        read(parts(7), *, iostat=io_status) grids(i)%dz
        
        if (io_status /= 0) then
            print *, "Error reading dz on line", i
            exit
        end if
        read(parts(8), *, iostat=io_status) grids(i)%nx
        
        if (io_status /= 0) then
            print *, "Error reading nx on line", i
            exit
        end if
        read(parts(9), *, iostat=io_status) grids(i)%ny
        
        if (io_status /= 0) then
            print *, "Error reading ny on line", i
            exit
        end if
        read(parts(10), *, iostat=io_status) grids(i)%nz
        
        if (io_status /= 0) then
            print *, "Error reading nz on line", i
            exit
        end if
        grids(i)%outputfile = prefix // trim(parts(11)) // ".d"
        
        ! Calculate total size needed
        ntot = grids(i)%nx * grids(i)%ny * grids(i)%nz
        
        ! Allocate with size checking
        if (ntot <= 0) then
            print *, "Invalid grid dimensions:", grids(i)%nx, grids(i)%ny, grids(i)%nz
            exit
        end if
        
        allocate(grids(i)%gridcoords(ind, 3), stat=io_status)
        allocate(grids(i)%gridcoords(ind, 3), stat=io_status)
        if (io_status /= 0) then
            print *, "Error allocating grid"
            exit
        end if
        
        print *, "Grid allocated successfully with size:", ntot 
        
        ! Fill the grid using integer indices instead of real-valued loops
        ind = 1
        do k = 0, grids(i)%nz - 1
            do j = 0, grids(i)%ny - 1
                do l = 0, grids(i)%nx - 1
                    if (ind > size(grids(i)%gridcoords, 1)) then
                        print *, "Index out of bounds:", ind
                        exit
                    end if
                    
                    grids(i)%gridcoords(ind, 1) = grids(i)%xmin + l * grids(i)%dx
                    grids(i)%gridcoords(ind, 2) = grids(i)%ymin + j * grids(i)%dy
                    grids(i)%gridcoords(ind, 3) = grids(i)%zmin + k * grids(i)%dz
                    ind = ind + 1
                end do
            end do
        end do
        call cartesian_to_spherical(grids(i)%gridcoords, grids(i)%sphericalgrid, ntot) 
    end do
    
    
    close(9)

    ! -> checked with processline.f90 
  end subroutine process_line 

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

  subroutine routine(line)
    character(len=*), intent(in) :: line
    real(real64) :: gridparams(10)
    character(len=100) :: output_file
    real(real64), allocatable :: cartgrid(:,:), sphercalgrid(:,:), interpolated_data(:,:)
    integer :: nx, ny, nz, ntot, i
    
    call process_line(line, gridparams, output_file)
    call cartesiantospherical(cartgrid, sphercalgrid)
    
    allocate(interpolated_data(size(cartgrid, 1), 24))
    
    !$OMP PARALLEL DO
    do i = 1, size(sphercalgrid, 1)
      interpolated_data(i,:) = interpolationpoint(sphercalgrid(i,:))
    end do
    !$OMP END PARALLEL DO
    
    ntot = nx * ny * nz
    
    ! Write output to binary file
    open(unit=10, file=trim(output_file), form='unformatted', access='stream')
    write(10) real(nx, real64), real(ny, real64), real(nz, real64), real(ntot, real64)
    write(10) cartgrid
    write(10) interpolated_data
    close(10)
  end subroutine routine

end module interpolation_module

program main
  use interpolation_module
  implicit none
  
  character(len=200) :: line
  integer :: io_status
  
  call load_data()
  call create_kdtree()
  
  open(unit=11, file='/data/sjammi6/thesisproject/grids_bh_disk_patrik', status='old', action='read')
  
  do
    read(11, '(A)', iostat=io_status) line
    if (io_status /= 0) exit
    call routine(line)
  end do
  
  close(11)
  
  print *, "FINISHED INTERPOLATION :)"
end program main
