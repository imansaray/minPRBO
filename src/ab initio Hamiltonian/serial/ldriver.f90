subroutine smlsauber( c, u, kk, nk, nkx, nky, nkz, ng, nb, nc, nv, nkret, kret, ladprec, mr, mmr )
  implicit none
  !
  integer :: nk, nkx, nky, nkz, ng, nb, nc, nv, nkret
  integer :: kk( nk, 3 ), kret( nkret )
  real( kind = kind( 1.0 ) ) :: mr( nkret, ng, ng )
  real( kind = kind( 1.0d0 ) ) :: mmr( nkret, ng, ng )
  complex( kind = kind( 1.0d0 ) ) :: u( nb, nk, ng ), c( nv, nc, nk )
  character * 6 :: ladprec
  !
  integer :: lx, hx, ly, hy, lz, hz, idwrk, xstep, ixl, ixh, nn
  real( kind = kind( 1.0d0 ) ), allocatable :: fr( :, :, : ), fi( :, :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: wrk( : ), vv( : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: e( :, :, : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: aa( :, :, : ), bt( :, : )
  !
  open( unit=99, file='ladcap', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) lx, ly, lz, hx, hy, hz
  close( unit=99 )
  nn = ( 1 + hx - lx ) * ( 1 + hy - ly ) * ( 1 + hz - lz ) ! alloc keyword
  call sizereport( 8 * nn, 'fr........' ); allocate( fr( lx : hx, ly : hy, lz : hz ) )
  call sizereport( 8 * nn, 'fi........' ); allocate( fi( lx : hx, ly : hy, lz : hz ) )
  idwrk = 2 * max( nkx * ( nkx + 1 ), nky * ( nky + 1 ), nkz * ( nkz + 1 ) )
  call sizereport( 8 * idwrk, 'wrk.......' ); allocate( wrk( idwrk ) )
  call sizereport( 8 * nkret, 'vv........' ); allocate( vv( nkret ) )
  call sizereport( 16 * nv * nk * ng, 'aa........' ); allocate( aa( nv, nk, ng ) )
  call sizereport( 16 * nv * ng, 'bt........' ); allocate( bt( nv, ng ) )
  open( unit=99, file='stepper.h', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) xstep
  close( unit=99 )
  call lpart1( ng, nk, nv, nc, u, c, aa )
  c = 0
  do ixl = 1, ng, xstep
     ixh = ixl + xstep - 1
     if ( ixh .gt. ng ) ixh = ng
     open( unit=99, file='ladstep', form='formatted', status='unknown' )
     rewind 99
     write ( 99, '(4i8)' ) ixl, ixh, ng, xstep
     close( unit=99 )
     call sizereport( 16 * ng * xstep * nk, 'e.........' ); allocate( e( 1 : ng, ixl : ixh, nk ) )
     call lpart2( ng, ixl, ixh, nk, nv, nc, nkret, lx, hx, ly, hy, lz, hz, idwrk, kk, kret, &
          ladprec, mr, mmr, vv, fr, fi, wrk, u, aa, e )
     call lpart3( nk, ixl, ixh, ng, nv, nc, u, c, e, bt )
     deallocate( e ); call sizereport( 0, 'e.........' )
  end do
  deallocate( vv ); call sizereport( 0, 'vv........' )
  deallocate( fr ); call sizereport( 0, 'fr........' )
  deallocate( fi ); call sizereport( 0, 'fi........' )
  deallocate( wrk ); call sizereport( 0, 'wrk.......' )
  deallocate( aa ); call sizereport( 0, 'aa........' )
  deallocate( bt ); call sizereport( 0, 'bt........' )
  !
  return
end subroutine smlsauber
