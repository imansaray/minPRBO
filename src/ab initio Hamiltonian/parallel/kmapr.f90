subroutine kmapr( nk, nkx, nky, nkz, xk, kk )
  implicit none
  !
  integer nk, nkx, nky, nkz
  integer :: kk( nk, 3 )
  double precision :: xk( nk, 3 )
  !
  integer :: j, ik, ikx, iky, ikz, ick, ikv( 3 ), nkv( 3 ), kl( 3 ), kh( 3 )
  double precision :: c( 3 ), xck( 3 ), pi
  !
  pi = 4.0d0 * atan( 1.0d0 )
  open( unit=99, file='ladcap', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) kl( : ), kh( : )
  close ( unit=99 )
  open( unit=99, file='kpts.dat', form='formatted', status='unknown' )
  rewind 99
  nkv( 1 ) = nkx; nkv( 2 ) = nky; nkv( 3 ) = nkz
  c( : ) = 2.0d0 * pi / dble( nkv( : ) )
  ik = 0
  do ikx = 1, nkx
     do iky = 1, nky
        do ikz = 1, nkz
           ikv( 1 ) = ikx - 1; ikv( 2 ) = iky - 1; ikv( 3 ) = ikz - 1
           ik = ik + 1
           xk( ik, : ) = c( : ) * dble( ikv( : ) )
           kk( ik, : ) = ikv( : )
           read ( 99, * ) ick, xck( : )
           if ( ick .ne. ik ) stop 'bad ick'
           ! mapping k points onto (approximately) zero-centered grid
           do j = 1, 3
              if ( abs( xck( j ) - xk( ik, j ) ) .gt. 0.0001d0 ) stop 'bad kpt'
              if ( ikv( j ) .gt. kh( j ) ) kk( ik, j ) = ikv( j ) - nkv( j )
           end do
        end do
     end do
  end do
  close( unit=99 )
  !
  return
end subroutine kmapr
