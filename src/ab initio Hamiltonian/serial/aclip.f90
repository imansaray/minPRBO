subroutine zfrxy( z, x, y, allow, n )
  implicit none
  !
  integer :: n
  real( kind = kind( 1.0d0 ) ) :: x( n ), y( n ), allow( n )
  complex( kind = kind( 1.0d0 ) ) :: z( n )
  !
  complex( kind = kind( 1.0d0 ) ) :: rm1
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  z( : ) = allow( : ) * ( x( : ) + rm1 * y( : ) )
  !
  return
end subroutine zfrxy
!
subroutine ztoxy( z, x, y, allow, n )
  implicit none
  !
  integer :: n
  real( kind = kind( 1.0d0 ) ) :: x( n ), y( n ), allow( n )
  complex( kind = kind( 1.0d0 ) ) :: z( n )
  !
  complex( kind = kind( 1.0d0 ) ) :: rm1
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  x( : ) = allow( : ) * z( : )
  y( : ) = -rm1 * allow( : ) * z( : )
  !
  return
end subroutine ztoxy
