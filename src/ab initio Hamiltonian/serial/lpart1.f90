subroutine lpart1( nx, nk, nv, nc, u, c, a )
  implicit none
  !
  integer :: nx, nk, nv, nc
  complex( kind = kind( 1.0d0 ) ) :: u( nv + nc, nk, nx ), c( nv, nc, nk )
  complex( kind = kind( 1.0d0 ) ) :: a( nv, nk, nx )
  !
  integer :: ix, ik, j
  !
  a=0
  do ix = 1, nx
     do ik = 1, nk
        do j = 1, nc
           a( :, ik, ix ) = a( :, ik, ix ) + u( nv + j, ik, ix ) * c( :, j, ik )
        end do
     end do
  end do
! write ( 6, * ) 'here we are in part 1'
  !
  return
end subroutine lpart1
