subroutine uphase( u, nb, nk, ng, ngx, ngy, ngz, xk )
  implicit none
  !
  integer :: nb, nk, ng, ngx, ngy, ngz
  real( kind = kind( 1.0d0 ) ) :: xk( nk, 3 )
  complex( kind = kind( 1.0d0 ) ) :: u( nb, nk, ng )
  !
  integer :: ix, ik
  real( kind = kind( 1.0d0 ) ) :: cc, ss, arg
  complex( kind = kind( 1.0d0 ) ) :: rm1
  real( kind = kind( 1.0d0 ) ) :: x( 3, ng )
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  call formx( ngx, ngy, ngz, x )
  do ix = 1, ng
     do ik = 1, nk
        arg = sum( x( :, ix ) * xk( ik, : ) )
        cc = dcos( arg )
        ss = dsin( arg )
        u( :, ik, ix ) = u( :, ik, ix ) * ( cc + rm1 * ss )
     end do
  end do
  !
  return
end subroutine uphase
