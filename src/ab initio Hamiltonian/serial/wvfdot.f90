function wvfdot( w, z, allow, n )
  implicit none
  !
  integer :: n
  double precision :: allow( n )
  double complex :: wvfdot, w( n ), z( n )
  !
  wvfdot = sum( conjg( w( : ) ) * z( : ) * allow( : ) )
!  wvfdot = sum( conjg( w( : ) ) * z( : ) )
  return
end function wvfdot
