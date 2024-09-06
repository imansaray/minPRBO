subroutine lpart3( nk, ixl, ixh, nx, nv, nc, u, c, e, bt )
  implicit none
  !
  integer :: nk, ixl, ixh, nx, nv, nc
  complex( kind = kind( 1.0d0 ) ) :: u( nv + nc, nk, nx ), c( nv, nc, nk )
  complex( kind = kind( 1.0d0 ) ) :: e( 1 : nx, ixl : ixh, nk ), bt( nv, nx )
  !
  integer :: ik, ix, j
  real( kind = kind( 1.0d0 ) ) :: ur, ui
  complex( kind = kind( 1.0d0 ) ) :: rm1
  !
! write ( 6, * ) ' in lpart3', ixl, ixh
  rm1 = -1
  rm1 = sqrt( rm1 )
  do ik = 1, nk
     bt = 0
     do ix = ixl, ixh
        do j = 1, nx
           bt( 1 : nv, ix ) = bt( 1 : nv, ix ) + e( j, ix, ik ) * u( 1 : nv, ik, j )
        end do
     end do
     do j = 1, nc
        do ix = ixl, ixh
           ur = u( j + nv, ik, ix )
           ui = -rm1 * u( j + nv, ik, ix )
           c( :, j, ik ) = c( :, j, ik ) + ( ur - rm1 * ui ) * bt( :, ix )
        end do
     end do
  end do
  !
  return
end subroutine lpart3
