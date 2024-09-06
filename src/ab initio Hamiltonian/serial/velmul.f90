subroutine velmul( vec, mul, n )
  implicit none
  !
  integer :: n
  real( kind = kind( 1.0d0 ) ) :: vec( n ), mul( n )
  !
  vec( : ) = vec( : ) * mul( : )
  !
  return
end subroutine velmul
