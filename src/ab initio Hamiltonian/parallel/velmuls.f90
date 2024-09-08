subroutine velmuls( vec, v2, mul, n, nn, ii )
  implicit none
  !
  integer :: n, nn
  integer :: ii( nn )
  real( kind = kind( 1.0d0 ) ) :: vec( n ), mul( nn ), v2( nn )
  !
  integer :: i
  !
  do i = 1, nn
     v2( i ) = vec( ii( i ) ) * mul( i )
  end do
  vec( : ) = 0
  do i = 1, nn
     vec( ii( i ) ) = v2( i )
  end do
  !
  return
end subroutine velmuls
