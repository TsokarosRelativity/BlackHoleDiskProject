program processline 
  use iso_fortran_env, only: real64, int64
  implicit none
  integer(int64) :: i, j, io_status, n_parts
  real(real64), allocatable :: new_grid(:)
  character(len=250) :: line, filename, errmsg
  character(len=40) :: tmp, prefix
  !character(len=40), dimension(11) :: parts(:)
  character(len=250), dimension(11) :: parts
  character(len=40), allocatable :: output_file
 
  type :: gridparams 
     real(real64) :: xmin, ymin, zmin, dz, dy, dx
     integer(int64) :: nz, ny, nx
     character(len=40) :: outputfile
     real(real64), allocatable :: grid(:,:)
  end type gridparams 

  type(gridparams), dimension(5) :: grids
  filename = "grids_bh_disk_patrik"
  prefix = "CTS_bin-proc"

  open(unit=9, file=filename, status='old', action='read', iostat=io_status)
  if (io_status /= 0) then
    print *, "Error opening file: ", filename
    stop
  end if
  
  do i = 1, 5 
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

end program processline
