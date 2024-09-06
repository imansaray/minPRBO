subroutine lpart2( nx, ixl, ixh, nk, nv, nc, nkret, lx, hx, ly, hy, lz, hz, idwrk, kk, kret, &
     ladprec, mr, mmr, vv, fr, fi, wrk, u, a, e )
  implicit none
  !
  integer :: nx, ixl, ixh, nk, nv, nc, nkret, lx, hx, ly, hy, lz, hz, idwrk, ii
  integer :: kk( nk, 3 ), kret( nkret )
  character( len=6 ) :: ladprec
  real( kind = kind( 1.0 ) ) :: mr( nkret, nx, nx )
  real( kind = kind( 1.0d0 ) ) :: mmr( nkret, nx, nx )
  real( kind = kind( 1.0d0 ) ) :: vv( nkret )
  real( kind = kind( 1.0d0 ) ) :: fr( lx : hx, ly : hy, lz : hz ), fi( lx : hx, ly : hy, lz : hz ), wrk( idwrk )
  complex( kind = kind( 1.0d0 ) ) :: u( nv + nc, nk, nx ), a( nv, nk, nx ), e( 1 : nx, ixl : ixh, nk )
  !
  integer :: ix, iy, ik
  real( kind = kind( 1.0d0 ) ) :: re, im
  complex( kind = kind( 1.0d0 ) ) :: rm1, z
  real( kind = kind( 1.0d0 ) ), allocatable :: muse( : )
  !
  integer * 8 :: plan1, plan2
  integer :: n1, n2, n3
  complex( kind = kind( 1.0d0 ) ), allocatable ::  var( :, :, : )
  !
  call sizereport( 8 * nkret, 'muse......' ); allocate( muse( nkret ) )
  rm1 = -1
  rm1 = sqrt( rm1 )
  n1 = hx - lx + 1
  n2 = hy - ly + 1
  n3 = hz - lz + 1
  allocate( var( n1, n2, n3 ) )
  call dfftw_plan_dft_3d( plan1, n1, n2, n3, var, var, -1, 32 )
  call dfftw_plan_dft_3d( plan2, n1, n2, n3, var, var, +1, 32 )
  ii = 0
  do ix = ixl, ixh
     do iy = 1, nx
        do ik = 1, nk
           z = dot_product( u( 1 : nv, ik, iy ), a( :, ik, ix ) )
           fr( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) = z
           fi( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) = -rm1 * z
        end do
        select case( ladprec )
        case( 'single' )
           do ik = 1, nkret
              muse( ik ) = mr( ik, iy, ix )
           end do
        case( 'double' )
           do ik = 1, nkret
              muse( ik ) = mmr( ik, iy, ix )
           end do
        end select 
!       call ctrdfft( fr, fi, lx, hx, ly, hy, lz, hz, -1, wrk, idwrk )
        call ctrdfftw( fr, fi, lx, hx, ly, hy, lz, hz, -1, var, plan1, rm1 )
!       call velmuls( fr, vv, mr( 1, iy, ix ), nk, nkret, kret )
!       call velmuls( fi, vv, mr( 1, iy, ix ), nk, nkret, kret )
        call velmuls( fr, vv, muse, nk, nkret, kret )
        call velmuls( fi, vv, muse, nk, nkret, kret )
!       call ctrdfft( fr, fi, lx, hx, ly, hy, lz, hz, +1, wrk, idwrk )
        call ctrdfftw( fr, fi, lx, hx, ly, hy, lz, hz, +1, var, plan2, rm1 )
        do ik = 1, nk
           re = fr( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) )
           im = fi( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) )
           e( iy, ix, ik ) = re + rm1 * im
        end do
     end do
     ii = ii + 1
     if ( ii .eq. 20 ) then
        open( unit=99, file='lprog', form='formatted', status='unknown' )
        rewind 99
        write ( 99, '(3i8)' ) ix, ixl, ixh
        close( unit=99 )
        ii = 0
     end if
  end do
  deallocate( muse ); call sizereport( 0, 'muse......' )
  !
  return
end subroutine lpart2
