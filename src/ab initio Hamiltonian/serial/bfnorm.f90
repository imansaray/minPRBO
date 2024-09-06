subroutine bfnorm( psi, eband, nk, nb, nbc, nbv )
  implicit none
  !
  integer :: nk, nb, nbc, nbv
  !
  complex( kind = kind( 1.0d0 ) ) :: psi( nbv, nbc, nk )
  real( kind = kind( 1.0d0 ) ) :: eband( nb, nk )
  !
  integer :: i, j, ik
  real( kind = kind( 1.0d0 ) ) :: ec, ev
  ! 
  do i = 1, nbv
     do j = 1, nbc
        do ik = 1, nk
            ec = eband( j + nbv, ik )
            ev = eband( i, ik )
            psi( i, j, ik ) = psi( i, j, ik ) * sqrt( 2.0d0 * abs( ec - ev ) )
        end do
     end do
  end do
  !
  return
end subroutine bfnorm
