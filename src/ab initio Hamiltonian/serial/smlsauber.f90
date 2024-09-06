subroutine smlsauber( c, u, kk, nk, nkx, nky, nkz, ng, nb, nc, nv, nkret, kret, mr )
  implicit none
  !
  integer :: nk, nkx, nky, nkz, ng, nb, nc, nv, nkret
  integer :: kk( nk, 3 ), kret( nkret )
  real( kind = kind( 1.0d0 ) ) :: mr( nkret, ng, ng )
  complex( kind = kind( 1.0d0 ) ) :: u( nb, nk, ng ), c( nv, nc, nk )
  !
  integer :: i, j, ix, iy, ik, lx, hx, ly, hy, lz, hz, idwrk, xstep, ixl, ixh
  real( kind = kind( 1.0d0 ) ) :: tr, ti, ur, ui
  complex( kind = kind( 1.0d0 ) ) :: r, rm1
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: fr( :, :, : ), fi( :, :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: wrk( : ), vv( : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: e( :, :, : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: a( :, :, : ), bt( :, : )
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  open( unit=99, file='ladcap', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) lx, ly, lz, hx, hy, hz
  close( unit=99 )
  allocate( fr( lx : hx, ly : hy, lz : hz ), fi( lx : hx, ly : hy, lz : hz ) )
  idwrk = 2 * max( nkx * ( nkx + 1 ), nky * ( nky + 1 ), nkz * ( nkz + 1 ) )
  allocate( wrk( idwrk ), vv( nkret ) )
  allocate( a( nv, nk, ng ), bt( nv, ng ) )
  !
  open( unit=99, file='stepper.h', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) xstep
  close( unit=99 )
  !
  a=0
  do ix = 1, ng
     do ik = 1, nk
        do j = 1, nc
           a( :, ik, ix ) = a( :, ik, ix ) + u( nv + j, ik, ix ) * c( :, j, ik )
        end do
     end do
  end do
  c = 0
  do ixl = 1, ng, xstep
     ixh = ixl + xstep - 1
     if ( ixh .gt. ng ) ixh = ng
     allocate( e( 1 : ng, ixl : ixh, nk ) )
     do ix = ixl, ixh
!IM        open( unit=99, file='lrprog', form='formatted', status='unknown' )
!IM        rewind 99
!IM        write ( 99, '(2i6)' ) ix, ng
!IM        close( unit=99 )
        do iy = 1, ng
           do ik = 1, nk
              r = 0
              do i = 1, nv
                 ur = u( i, ik, iy )
                 ui = -rm1 * u( i, ik, iy )
                 r = r + ( ur - rm1 * ui ) * a( i, ik, ix )
              end do
              fr( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) = r
              fi( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) = -rm1 * r
           end do
           call ctrdfft( fr, fi, lx, hx, ly, hy, lz, hz, -1, wrk, idwrk )
           call velmuls( fr, vv, mr( 1, iy, ix ), nk, nkret, kret )
           call velmuls( fi, vv, mr( 1, iy, ix ), nk, nkret, kret )
           call ctrdfft( fr, fi, lx, hx, ly, hy, lz, hz, +1, wrk, idwrk )
           do ik = 1, nk
              tr = fr( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) )
              ti = fi( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) )
              e( iy, ix, ik ) = tr + rm1 * ti
           end do
        end do
     end do
     do ik = 1, nk
        bt = 0
        do ix = ixl, ixh
           do j = 1, ng
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
     deallocate( e )
  end do
  !
  deallocate( fr, fi, wrk, a, bt )
  !
  return
end subroutine smlsauber
