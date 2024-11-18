module comp
  use iso_fortran_env, only : real64, int64
  implicit none

contains
  function dist(p1, p2) result(d)
    real(real64), intent(in) :: p1(3), p2(3)
    real(real64) :: d

    d = sqrt((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2)

  end function dist

  function nn(pdata, q) result(i)
    real(real64), intent(in) :: pdata(:,:), q(3)
    integer(int64) :: i, n 
    real(real64) :: d, mindist

    mindist = huge(mindist)
    do n = 1, size(pdata, 1)
      d = dist(pdata(n, :), q)
      if (d < mindist) then
        mindist = d
        i = n
      end if
    end do

  end function nn
end module comp

program nearestneighbors
  use iso_fortran_env, only : real64, int64
  use comp
  implicit none
  real(real64) :: pointdata(300, 3)
  integer(int64) :: i, j, k, ind, n
  real(real64) :: query(3)

  n = 1
  print *, "populating point data"
  do i = 1, 10
    do j = 1, 100, 10
      do k = 30, 90, 30
        pointdata(n, 1) = real(i, real64)
        pointdata(n, 2) = real(j, real64)
        pointdata(n, 3) = real(k, real64)
        print *, n, pointdata(n, 1), pointdata(n, 2), pointdata(n, 3)
        n = n + 1
      end do
    end do 
  end do
  
  query(1) = real(5, real64)
  query(2) = real(50, real64)
  query(3) = real(60, real64)
  ind = nn(pointdata, query)
  print *, "nearest neighbor to query point"
  print *, ind, pointdata(ind, :)
  print *, "query point"
  print *, query


!  do i = 1, 10
!    print *, i
!  end do 
!
!  print *, "------"
!  do j = 1, 100, 10
!    print *, j
!  end do
!
!  print *, "------"
!  do k = 30, 90, 30
!    print *, k
!  end do 
  

end program nearestneighbors
