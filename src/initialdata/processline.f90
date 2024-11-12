program processline 
    use graysmod
    implicit none
    integer(int64) :: i, j, k, l, io_status, n_parts, ind
    character(len=250) :: line, filename, errmsg
    character(len=40) :: tmp, prefix
    character(len=250), dimension(11) :: parts
    character(len=40), allocatable :: output_file
    real(real64) :: xmax, ymax, zmax
    
    type :: gridparams 
        real(real64) :: xmin, ymin, zmin, dz, dy, dx
        integer(int64) :: nz, ny, nx
        character(len=40) :: outputfile
        real(real64), allocatable :: gridcoords(:,:)
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
            print *, "Error reading xmin"
            exit
        end if
        read(parts(3), *, iostat=io_status) grids(i)%ymin
        read(parts(4), *, iostat=io_status) grids(i)%zmin
        read(parts(5), *, iostat=io_status) grids(i)%dx
        read(parts(6), *, iostat=io_status) grids(i)%dy
        read(parts(7), *, iostat=io_status) grids(i)%dz
        read(parts(8), *, iostat=io_status) grids(i)%nx
        read(parts(9), *, iostat=io_status) grids(i)%ny
        read(parts(10), *, iostat=io_status) grids(i)%nz
        
        grids(i)%outputfile = prefix // trim(parts(11)) // ".d"
        
        ! Calculate total size needed
        ind = grids(i)%nx * grids(i)%ny * grids(i)%nz
        
        ! Allocate with size checking
        if (ind <= 0) then
            print *, "Invalid grid dimensions:", grids(i)%nx, grids(i)%ny, grids(i)%nz
            exit
        end if
        
        allocate(grids(i)%gridcoords(ind, 3), stat=io_status)
        if (io_status /= 0) then
            print *, "Error allocating grid"
            exit
        end if
        
        print *, "Grid allocated successfully with size:", ind
        
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
        
        ! Print first few points for verification
        print *, "Printing first few grid points for verification:"
        do ind = 1, min(10, size(grids(i)%gridcoords, 1))
            print *, "grid point index ", ind
            print *, grids(i)%gridcoords(ind, :)
        end do
        
    end do
    
    ! Clean up
    do i = 1, 5
        if (allocated(grids(i)%gridcoords)) then
            deallocate(grids(i)%gridcoords)
        end if
    end do
    
    close(9)
end program processline

