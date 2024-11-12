program datareader
  use iso_fortran_env, only: real64, int64
  implicit none

  real(real64) :: data(16, 3)
  integer(int64) :: i, j
  character(len=200) :: filename, line, outputfile, ios

  open(unit=10, file='data.txt, status='old, action='read')




end program datareader
