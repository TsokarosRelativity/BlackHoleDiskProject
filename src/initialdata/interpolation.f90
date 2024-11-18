module interpolation_module
  use iso_fortran_env, only: real64, int64
  implicit none

  ! Global variables
  real(real64), allocatable :: rad_arr(:), theta_arr(:), phi_arr(:)
  real(real64), allocatable :: data_3d(:,:)
  type :: gridparams 
     real(real64) :: xmin, ymin, zmin, dz, dy, dx
     integer(int64) :: nz, ny, nx
     character(len=40) :: outputfile
     real(real64), allocatable :: gridcoords(:,:), sphericalgrid(:,:)
  end type gridparams 
  type(gridparams), dimension(1280) :: grids

contains

! tested and works  
  subroutine load_data()
    character(len=100) :: filename
    integer :: io_status, ind, rows
    ! Implement data loading here
    print *, "Loading in 3D data..."
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

    print *, "3D Data loaded successfully!"
    print *, "------------------------------"
    print *, "Loading in Rad Array..."
    filename = "radarr.txt"
    open(unit=7, file=filename, status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
      print *, "Error opening file: ", filename
      exit
    endif 
    
    rows = 1604 !matrix dimensions -> (1604, ) 
    allocate(rad_arr(rows), stat=io_status)
    if (io_status /= 0) then 
      print *, "Error allocating data"
      exit
    end if
    
    do ind = 1, row
      read(7, '(ES16.6)', stat=io_status) rad_arr(ind)
      if (io_status /= 0) then
        print *, "Error reading data on line " , ind
        exit
      end if
    end do

    print *, "Rad Array Data loaded successfully!"

    print *, "------------------------------"
    print *, "Loading in Theta Array..."
    filename = "thetaarr.txt"
    open(unit=10, file=filename, status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
      print *, "Error opening file: ", filename
      exit
    endif 
    
    rows = 403 !matrix dimensions -> (403, ) 
    allocate(theta_arr(rows), stat=io_status)
    if (io_status /= 0) then 
      print *, "Error allocating data"
      exit
    end if
    
    do ind = 1, row
      read(10, '(ES16.6)', stat=io_status) theta_arr(ind)
      if (io_status /= 0) then
        print *, "Error reading data on line " , ind
        exit
      end if
    end do

    print *, "Theta Array Data loaded successfully!"

    print *, "------------------------------"
    print *, "Loading in Phi Array..."
    filename = "phiarr.txt"
    open(unit=11, file=filename, status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
      print *, "Error opening file: ", filename
      exit
    endif 
    
    rows = 100 !matrix dimensions -> (100, ) 
    allocate(phi_arr(rows), stat=io_status)
    if (io_status /= 0) then 
      print *, "Error allocating data"
      exit
    end if
    
    do ind = 1, row
      read(11, '(ES16.6)', stat=io_status) phi_arr(ind)
      if (io_status /= 0) then
        print *, "Error reading data on line " , ind
        exit
      end if
    end do

    print *, "Phi Array Data loaded successfully!"

  end subroutine load_data

  subroutine distanceEuclidian(p1, p2, d)
    real(real64), intent(in) :: p1, p2
    real(real64), intent(out) :: d
    d = sqrt((p1 - p2)**2)
  end subroutine distanceEuclidian

  subroutine distanceEuclidian2d(p1, p2, d)
    real(real64), intent(in) :: p1, p2
    real(real64), intent(out) :: d
    d = sqrt((p1 - p2)**2)
  end subroutine distanceEuclidian2d 

  subroutine distanceEuclidian3d(p1, p2, d)
    real(real64), intent(in) :: p1(3), p2(3)
    real(real64), intent(out) :: d

    d = sqrt((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2) 

  end subroutine distanceEuclidian3d 

  subroutine nearestneighbors(pdata, p, n)
    real(real64), intent(in) :: pdata(:), p
    integer, intent(out) :: n
    integer :: i
    real(real64) :: d, dmin
    
    dmin = huge(dmin)
    do i = 1, size(pdata, 1)
      call distanceEuclidian(pdata(i), p, d)
      if (d < dmin) then
        dmin = d
        n = i
      end if
    end do

  end subroutine nearestneighbors

  subroutine nearestneighbors2d(pdata, p, n)
    real(real64), intent(in) :: pdata(:,:), p(2)
    integer, intent(out) :: n
    integer :: i
    real(real64) :: d, dmin
    
    dmin = huge(dmin)
    do i = 1, size(pdata, 1)
      call distanceEuclidian2d(pdata(i, :), p, d)
      if (d < dmin) then
        dmin = d
        n = i
      end if
    end do

  end subroutine nearestneighbors2d

  subroutine nearestneighbors3d(pdata, p, n)
    real(real64), intent(in) :: pdata(:,:), p(3)
    integer, intent(out) :: n
    integer :: i
    real(real64) :: d, dmin
    
    dmin = huge(dmin)
    do i = 1, size(pdata, 1)
      call distanceEuclidian3d(pdata(i, :), p, d)
      if (d < dmin) then
        dmin = d
        n = i
      end if
    end do

  end subroutine nearestneighbors3d

  subroutine process_line()
    integer(int64) :: i, j, k, l, io_status, n_parts, ntot, ind
    character(len=250)  :: line
    character(len=250) :: filename, errmsg, tmp, prefix
    character(len=250), dimension(11) :: parts
    real(real64) :: xmax, ymax, zmax
    
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
        read(parts(2), '(ES16.6)', iostat=io_status) grids(i)%xmin
        if (io_status /= 0) then
            print *, "Error reading xmin on line ", i
            exit
        end if
        read(parts(3), '(ES16.6)', iostat=io_status) grids(i)%ymin

        if (io_status /= 0) then
            print *, "Error reading ymin on line", i
            exit
        end if
        read(parts(4), '(ES16.6)', iostat=io_status) grids(i)%zmin
        
        if (io_status /= 0) then
            print *, "Error reading zmin on line", i
            exit
        end if
        read(parts(5), '(ES16.6)', iostat=io_status) grids(i)%dx
        
        if (io_status /= 0) then
            print *, "Error reading dx on line", i
            exit
        end if
        read(parts(6), '(ES16.6)', iostat=io_status) grids(i)%dy
        
        if (io_status /= 0) then
            print *, "Error reading dy on line", i
            exit
        end if
        read(parts(7), '(ES16.6)', iostat=io_status) grids(i)%dz
        
        if (io_status /= 0) then
            print *, "Error reading dz on line", i
            exit
        end if
        read(parts(8), '(I4)', iostat=io_status) grids(i)%nx
        
        if (io_status /= 0) then
            print *, "Error reading nx on line", i
            exit
        end if
        read(parts(9), '(I4)', iostat=io_status) grids(i)%ny
        
        if (io_status /= 0) then
            print *, "Error reading ny on line", i
            exit
        end if
        read(parts(10), '(I4)', iostat=io_status) grids(i)%nz
        
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
        
        allocate(grids(i)%gridcoords(ntot, 27), stat=io_status)
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
                        print *, "Index out of bounds:", ind, "which is bigger than ntot: ", ntot
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

  end subroutine process_line 

  subroutine cartesian_to_spherical(cart_coords, sph_coords, n_points)
      implicit none
      ! Input/output variables
      integer, intent(in) :: n_points
      real(real64), intent(in) :: cart_coords(n_points, 27)  ! Changed from real(int64)
      real(real64), intent(out) :: sph_coords(n_points, 27)  ! Changed from real(int64)
      
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

  
  subroutine interpolatepoint(point, interpolateddata)
    real(real64), dimension(3), intent(in) :: point
    real(real64), dimension(27), intent(out) :: interpolateddata
    integer :: rid, tid, pid, ind, i, j, k, a, c
    integer, dimension(5) :: ps, ts, rs
    real(real64) :: rp, thetap, phip, d1(5, 27), d2(27, 25), d3(27, 125), tmppoint(3), tmp, w1(25), w2(125) 
    real(real64), parameter :: r0 = 0.0124_real64
    logical :: rmatch, thetamatch, phimatch    
    
    rp = point(1)
    thetap = point(2)
    phip = point(3)

    call nearestneighbors(rad_arr, rp, rid)
    call nearestneighbors(theta_arr, thetap, tid)
    call nearestneighbors(phi_arr, phip, pid) 
    if (rp == rad_arr(rid)) then
      rmatch = .true.
    else
      rmatch = .false.
      if (rp > rad_arr(rid)) then
        rid = rid + 1
      end if
    end if

    if (thetap == theta_arr(tid)) then
      thetamatch = .true.
    else 
      thetamatch = .false.
      if (thetap > theta_arr(tid)) then 
        tid = tid + 1
        if (tid > size(theta_arr)) then
          tid = 1
        end if
      end if
    end if 

    if (phip == phi_arr(pid)) then
      phimatch = .true.
    else 
      phimatch = .false.
      if (phip > phi_arr(pid)) then 
        pid = pid + 1
        if (pid > size(phi_arr)) then
          pid = 1
        end if
      end if
    end if


    ! case 1: exact match
    if (rmatch .and. thetamatch .and. phimatch) then
      ind = 0
      call nearestneighbors3d(data_3d, point, ind)
      interpolateddata = data_3d(ind, :)
      return
    end if

    interpolateddata(1) = rp
    interpolateddata(2) = thetap
    interpolateddata(3) = phip
    ! case 2: phi not match
    if (rmatch .and. thetamatch .and. .not. phimatch) then
      if (pid < 3) then 
        ps = [3, size(phi_arr), 1, 2, size(phi_arr) - 1]
      else if (pid > size(phi_arr) - 2) then
        ps = [1, size(phi_arr) - 3, size(phi_arr) - 2, size(phi_arr) - 1, size(phi_arr)]
      else 
        ps = [pid - 2, pid - 1, pid, pid + 1, pid + 2]
      end if
      
      do i = 1, 5
        tmppoint(1) = rp
        tmppoint(2) = thetap
        tmppoint(3) = phi_arr(ps(i))
        call nearestneighbors3d(data_3d, tmppoint, ind)
        d1(i, :) = data_3d(ind, :)
      end do
      
      do i = 4, 27
        call lagcheby1_interp_1d(5, phi_arr(ps), d1(:, i), 1, phip, interpolateddata(i))
      end do
      return
    end if

    ! case 3: theta not match
    if (rmatch .and. .not. thetamatch .and. phimatch) then
      if (tid < 3) then 
        ts = [1, 2, 3, 4, 5]
      else if (tid > size(theta_arr) - 3) then
        ts = [size(theta_arr) - 4, size(theta_arr) - 3, size(theta_arr) - 2, size(theta_arr) - 1, size(theta_arr)]
      else 
        ts = [tid - 2, tid - 1, tid, tid + 1, tid + 2]
      end if
      ! interpolate
      do i = 1, 5
        tmppoint(1) = rp
        tmppoint(2) = theta_arr(ts(i))
        tmppoint(3) = phip
        call nearestneighbors3d(data_3d, tmppoint, ind)
        d1(i, :) = data_3d(ind, :)
      end do

      do i = 4, 27
        call lagcheby1_interp_1d(5, theta_arr(ts), d1(:, i), 1, thetap, interpolateddata(i))
      end do
      return
    end if


    ! case 4: r not match
    if (.not. rmatch .and. thetamatch .and. phimatch) then
      if (rid < 3) then 
        rs = [1, 2, 3, 4, 5]
      else if (rid > size(rad_arr) - 3) then
        rs = [size(rad_arr) - 4, size(rad_arr) - 3, size(rad_arr) - 2, size(rad_arr) - 1, size(rad_arr)]
      else 
        rs = [rid - 2, rid - 1, rid, rid + 1, rid + 2]
      end if
      
      do i = 1, 5
        tmppoint(1) = rad_arr(rs(i))
        tmppoint(2) = thetap
        tmppoint(3) = phip

        call nearestneighbors3d(data_3d, tmppoint, ind)
        d1(i, :) = data_3d(ind, :)
      end do
      do i = 4, 27
        call lagcheby1_interp_1d(5, rad_arr(rs), d1(:, i), 1, rp, interpolateddata(i))
      end do
      return
      ! interpolate
    end if

    ! case 5: theta, phi not match
    if (rmatch .and. .not. thetamatch .and. .not. phimatch) then
      if (tid < 3) then 
        ts = [1, 2, 3, 4, 5]
      else if (tid > size(theta_arr) - 3) then
        ts = [size(theta_arr) - 4, size(theta_arr) - 3, size(theta_arr) - 2, size(theta_arr) - 1, size(theta_arr)]
      else 
        ts = [tid - 2, tid - 1, tid, tid + 1, tid + 2]
      end if
      if (pid < 3) then 
        ps = [3, size(phi_arr), 1, 2, size(phi_arr) - 1]
      else if (pid > size(phi_arr) - 2) then
        ps = [1, size(phi_arr) - 3, size(phi_arr) - 2, size(phi_arr) - 1, size(phi_arr)]
      else 
        ps = [pid - 2, pid - 1, pid, pid + 1, pid + 2]
      end if

      c = 1
      do i = 1, 5
        do j = 1, 5
          tmppoint(1) = rp
          tmppoint(2) = theta_arr(ts(i))
          tmppoint(3) = phi_arr(ps(j))
          call nearestneighbors3d(data_3d, tmppoint, ind)
          d2 = TRANSPOSE(data_3d(ind, :))
          c = c + 1
        end do
      end do

      do k = 4, 27
        call rbf_weight(3, 25, d2(:3,:), r0, phi1, d2(k,:), w1)
        call rbf_interp_nd(3, 25, d2(:3, :), r0, phi1, w1, 1, TRANSPOSE(interpolateddata)(:3, :), interpolateddata(i))
      end do
      return
    end if

    ! case 6: r, phi not match
    if (.not. rmatch .and. thetamatch .and. .not. phimatch) then
      if (rid < 3) then 
        rs = [1, 2, 3, 4, 5]
      else if (rid > size(rad_arr) - 3) then
        rs = [size(rad_arr) - 4, size(rad_arr) - 3, size(rad_arr) - 2, size(rad_arr) - 1, size(rad_arr)]
      else 
        rs = [rid - 2, rid - 1, rid, rid + 1, rid + 2]
      end if
      if (pid < 3) then 
        ps = [3, size(phi_arr), 1, 2, size(phi_arr) - 1]
      else if (pid > size(phi_arr) - 2) then
        ps = [1, size(phi_arr) - 3, size(phi_arr) - 2, size(phi_arr) - 1, size(phi_arr)]
      else 
        ps = [pid - 2, pid - 1, pid, pid + 1, pid + 2]
      end if
      ! interpolate
      c = 1
      do i = 1, 5
        do j = 1, 5
          tmppoint(1) = rad_arr(rs(i))
          tmppoint(2) = thetap
          tmppoint(3) = phi_arr(ps(j))
          !tmppoint = [rad_arr(rs(i)), thetap, phi_arr(ps(j))]
          call nearestneighbors3d(data_3d, tmppoint, ind)
          d2 = TRANSPOSE(data_3d(ind, :))
          c = c + 1
        end do
      end do
      do k = 4, 27
        call rbf_weight(3, 25, d2(:3,:), r0, phi1, d2(k,:), w1)
        call rbf_interp_nd(3, 25, d2(:3, :), r0, phi1, w1, 1, TRANSPOSE(interpolateddata)(:3, :), interpolateddata(i))
      end do
      return
    end if


    ! case 7: r, theta not match
    if (.not. rmatch .and. .not. thetamatch .and. phimatch) then
      if (rid < 3) then 
        rs = [1, 2, 3, 4, 5]
      else if (rid > size(rad_arr) - 3) then
        rs = [size(rad_arr) - 4, size(rad_arr) - 3, size(rad_arr) - 2, size(rad_arr) - 1, size(rad_arr)]
      else 
        rs = [rid - 2, rid - 1, rid, rid + 1, rid + 2]
      end if
      if (tid < 3) then 
        ts = [1, 2, 3, 4, 5]
      else if (tid > size(theta_arr) - 3) then
        ts = [size(theta_arr) - 4, size(theta_arr) - 3, size(theta_arr) - 2, size(theta_arr) - 1, size(theta_arr)]
      else 
        ts = [tid - 2, tid - 1, tid, tid + 1, tid + 2]
      end if
      c = 1
      do i = 1, 5
        do j = 1, 5
          tmppoint(1) = rad_arr(rs(i))
          tmppoint(2) = theta_arr(ts(j))
          tmppoint(3) = phip
          !tmppoint = [rad_arr(rs(i)), theta_arr(ts(j)), phip]
          call nearestneighbors3d(data_3d, tmppoint, ind)
          d2 = TRANSPOSE(data_3d(ind, :))
          c = c + 1
        end do
      end do
!
      do k = 4, 27
        call rbf_weight(3, 25, d2(:3,:), r0, phi1, d2(k,:), w1)
        call rbf_interp_nd(3, 25, d2(:3, :), r0, phi1, w1, 1, TRANSPOSE(interpolateddata)(:3, :), interpolateddata(i))
      end do
      return
    end if

    ! case 8: r, theta, phi not match
    if (.not. rmatch .and. .not. thetamatch .and. .not. phimatch) then
      if (pid < 3) then 
        ps = [3, size(phi_arr), 1, 2, size(phi_arr) - 1]
      else if (pid > size(phi_arr) - 2) then
        ps = [1, size(phi_arr) - 3, size(phi_arr) - 2, size(phi_arr) - 1, size(phi_arr)]
      else 
        ps = [pid - 2, pid - 1, pid, pid + 1, pid + 2]
      end if
      if (rid < 3) then 
        rs = [1, 2, 3, 4, 5]
      else if (rid > size(rad_arr) - 3) then
        rs = [size(rad_arr) - 4, size(rad_arr) - 3, size(rad_arr) - 2, size(rad_arr) - 1, size(rad_arr)]
      else 
        rs = [rid - 2, rid - 1, rid, rid + 1, rid + 2]
      end if
      if (tid < 3) then 
        ts = [1, 2, 3, 4, 5]
      else if (tid > size(theta_arr) - 3) then
        ts = [size(theta_arr) - 4, size(theta_arr) - 3, size(theta_arr) - 2, size(theta_arr) - 1, size(theta_arr)]
      else 
        ts = [tid - 2, tid - 1, tid, tid + 1, tid + 2]
      end if

      c = 1
      do i = 1, 5
        do j = 1, 5
          do k = 1, 5
            tmppoint(1) = rad_arr(rs(i))
            tmppoint(2) = theta_arr(ts(j))
            tmppoint(3) = phi_arr(ps(k))
!            tmppoint = [rad_arr(rs(i)), theta_arr(ts(j)), phi_arr(ps(k))]
            call nearestneighbors3d(data_3d, tmppoint, ind)
            d3 = TRANSPOSE(data_3d(ind, :))
            c = c + 1
          end do
        end do
      end do
      do k = 4, 27
        call rbf_weight(3, 125, d2(:3,:), r0, phi3, d2(k,:), w2)
        call rbf_interp_nd(3, 125, d2(:3, :), r0, phi3, w2, 1, TRANSPOSE(interpolateddata)(:3, :), interpolateddata(i))
      end do
      return
    end if
  end subroutine interpolatepoint

  subroutine interpolategrid(grid)
    type(gridparams), intent(inout) :: grid
    integer :: i, npoints    
    real(real64) :: interpolateddata(27)

    npoints = grid%nx * grid%ny * grid%nz
    do i = 1, npoints
      call interpolatepoint(grid%sphericalgrid(i, :), interpolateddata)
      grid%sphericalgrid(i, :) = interpolateddata 
      grid%gridcoords(i, 4:) = interpolateddata(4:)
    end do
  end subroutine interpolategrid

  subroutine gridwriter(grid, iou)
    type(gridparams), intent(in) :: grid
    integer :: ntot
    integer, intent(in) :: iou
    
    ntot = grid%nx * grid%ny * grid%nz
    open(newunit=iou, file=trum(grid%outputfile), status='replace', form='unformatted', access='stream', action='write')
    
    write(iou) real(grid%nx, kind=real64), real(grid%ny, kind=real64), real(grid%nz, kind=real64), real(ntot, kind=real64) 

    write(iou) grid%gridcoords

    close(iou)
  end subroutine
  subroutine routine()
    integer :: i
    call read_data()
    call process_line()
    do i = 1, 1280 
      call interpolategrid(grids(i))
      call gridwriter(grids(i), i + 1000)
    end do
  end subroutine routine
end module interpolation_module

program main 
  use interpolation_module
  implicit none
  call routine()
end program main
! read data -> populated 3d data, rad arr, phi arr, theta arr
! process line -> populate grids
! interpolate point -> takes in spherical point and interpolates data
! interpolate grid -> interpolates points on grid
! write grid -> done
