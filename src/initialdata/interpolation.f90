module interpolation_module
  use iso_fortran_env, only: real64, int32
  implicit none

  ! Global variables
  real(real64), allocatable :: rad_arr(:), theta_arr(:), phi_arr(:)
  real(real64), allocatable :: data_3d(:,:)

  type :: gridparams 
     real(real64) :: xmin, ymin, zmin, dz, dy, dx
     integer(int64) :: nz, ny, nx, MPI_ID
     character(len=40) :: outputfile
     real(real64), allocatable :: grid(:,:)
  end type gridparams 

contains

  subroutine load_data()
    ! Implement data loading here
    print *, "Loading in data..."
    ! TODO: Implement HDF5 reading for df_3d
    ! TODO: Implement sorting and unique operations on df_3d
    ! TODO: Initialize rad_arr, theta_arr, phi_arr
    print *, "Data loaded successfully!"
  end subroutine load_data

  subroutine create_kdtree()
    print *, "Creating KD-tree..."
    ! TODO: Implement KD-tree creation
    print *, "KD-tree created successfully!"
  end subroutine create_kdtree

  function get_nearest_point_index(r, theta, phi) result(index)
    real(real64), intent(in) :: r, theta, phi
    integer :: index
    ! TODO: Implement KD-tree query
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
    character(len=250), intent(in) :: line
    type(gridparams), dimension(1280), intent(out) :: grids
    character(len=100) :: tmp, prefix, filename
    character(len=50), dimension(11) :: parts
    integer ::, i, j, n_parts, ios

    prefix = "CTs_bin-proc"
    filename = "grids_bh_disk_patrick"

    open(unit=9, file=filename, status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
      print *, "Error opening file: ", filename
      stop
    end if
    
    do i = 1, 1280 
      read(9, '(A)', iostat=io_status) line
      if (io_status /= 0) then
        print *, "error reading data"
      end if 
      print *, "line ", i, "data: "
      print *, line
      print *, "line length:"
      print *, len(trim(line))
      print *, line(1:len(trim(line)))
      tmp = ""
      n_parts = 1
      parts = ""
      do j = 1, len(trim(line))
        !print *, line(j:j)
        if (line(j:j) == ' ') then
          if (len(trim(tmp)) > 0) then
            !print *, "tmp: ", trim(tmp)
            parts(n_parts) = tmp
            n_parts = n_parts + 1
            tmp = ""
          end if
        else
          tmp = trim(tmp) // line(j:j)
          !print *, trim(tmp)
        end if
      end do
      parts(n_parts) = trim(tmp)
      grids(i)%xmin = real(parts(2), real64)
      grids(i)%ymin = real(parts(3), real64)
      grids(i)%zmin = real(parts(4), real64)
      grids(i)%dx = real(parts(5), real64)
      grids(i)%dy = real(parts(6), real64)
      grids(i)%dz = real(parts(7), real64)
      grids(i)%nx = int(parts(8), int64)
      grids(i)%ny = int(parts(9), int64)
      grids(i)%nz = int(parts(10), int64)
      grids(i)%outputfile = prefix // trim(parts(11)) // ".d"
      
    end do
    close(9)

    ! -> checked with processline.f90 
  end subroutine process_line

  subroutine gridmaker(gridparams, gridcoords, nx, ny, nz)
    real(real64), intent(in) :: gridparams(10)
    real(real64), allocatable, intent(out) :: gridcoords(:,:)
    integer, intent(out) :: nx, ny, nz
    ! TODO: Implement grid creation
  end subroutine gridmaker

  subroutine cartesiantospherical(cart, spher)
    real(real64), intent(in) :: cart(:,:)
    real(real64), allocatable, intent(out) :: spher(:,:)
    ! TODO: Implement coordinate transformation
  end subroutine cartesiantospherical

  subroutine routine(line)
    character(len=*), intent(in) :: line
    real(real64) :: gridparams(10)
    character(len=100) :: output_file
    real(real64), allocatable :: cartgrid(:,:), sphercalgrid(:,:), interpolated_data(:,:)
    integer :: nx, ny, nz, ntot, i
    
    call process_line(line, gridparams, output_file)
    call gridmaker(gridparams, cartgrid, nx, ny, nz)
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
