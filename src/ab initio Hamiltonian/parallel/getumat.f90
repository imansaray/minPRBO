subroutine getumat_orig( nb, nk, ng, u )
  implicit none
  !
  integer, parameter :: u2dat = 35
  !
  integer nb, nk, ng
  double complex u( nb, nk, ng )
  double complex, allocatable :: utmp(:)
  !
  integer ik, ib, ig, idm1, idm2, idm3
  double precision unorm
  !
  allocate( utmp( ng ) )
  open( unit=u2dat, file='u2.dat', form='unformatted', status='unknown' )
  rewind u2dat
  do ik = 1, nk
     do ib = 1, nb
        do ig = 1, ng
           read ( u2dat ) idm1, idm2, idm3, utmp( ig )
        end do
        unorm = sqrt( dble( dot_product( utmp, utmp ) ) )
        unorm = 1.0d0 / unorm
        utmp( : ) = utmp( : ) * unorm
        do ig = 1, ng
           u( ib, ik, ig ) = utmp( ig )
        end do
     end do
  end do
  close( unit=u2dat )
  !
  deallocate( utmp )
  !
  return
end subroutine getumat_orig

subroutine getumat( nb,ilftl,ilfth,irgtl,irgth,nk,ng,u )
  ! For reading the u2.dat from OCEAN using the wantLegacy=.True.      
  implicit none
  !
  integer, parameter :: u2dat = 35, conu2dat = 36
  !
  integer nb, nk, ng,ilftl,ilfth,irgtl,irgth
  double complex u( nb, nk, ng )
  double complex, allocatable :: utmp(:),valutmp(:),conutmp(:)
  !
  integer ik, ib, ig, idm1, idm2, idm3
  double precision unorm
  !
  allocate( utmp( ng ) )
  ! old way with all the u fcns (unphased) in u2.dat
  open( unit=u2dat, file='u2.dat', form='unformatted', status='unknown' )
  rewind u2dat
  do ik = 1, nk

  ! doing the valence bands first (ilftl:ilfth)
     do ib = ilftl, ilfth
        utmp = 0.0
        do ig = 1, ng
           read ( u2dat ) idm1, idm2, idm3, utmp( ig )
        end do
        unorm = sqrt( dble( dot_product( utmp, utmp ) ) )
        unorm = (1.0d0 / unorm) !/(1/10029.4937)
        utmp( : ) = utmp( : ) * unorm
        do ig = 1, ng
           u( ib, ik, ig ) = utmp( ig )
        end do
     end do


  ! doing the conduction bands first (irgtl:irgth)
     do ib = irgtl, irgth
        utmp = 0.0
        do ig = 1, ng
           read ( u2dat ) idm1, idm2, idm3, utmp( ig )
        end do
        unorm = sqrt( dble( dot_product( utmp, utmp ) ) )
        unorm = (1.0d0 / unorm)!/(1/10029.4937)
        utmp( : ) = utmp( : ) * unorm
        do ig = 1, ng
           u( ib, ik, ig ) = utmp( ig )
        end do
     end do


  end do
  close( unit=u2dat )
  !
  deallocate( utmp )
  !
  return
end subroutine getumat
