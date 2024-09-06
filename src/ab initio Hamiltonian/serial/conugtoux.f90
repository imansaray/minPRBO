 program conugtoux
   implicit none
   !
   integer, parameter :: u1dat = 22, u1xdat = 23
   integer :: nk, nb, nx( 3 ), ng, i, nwfile, iwfile,nspinor,iband
   integer, allocatable :: g( :, : ), ikl( :, : )
   real( kind = kind( 1.0d0 ) ), allocatable :: zr( :, : ), zi( :, : ),zlin(:)
   character * 20, allocatable :: wfnam( : )
   !
   open( unit=99, file='kandb.h', form='formatted', status='unknown' )
   call igetval( nk, 'nk   ' )
   call igetval( nb, 'nb   '  )
   call igetval( nx( 1 ), 'ngx  ' )
   call igetval( nx( 2 ), 'ngy  ' )
   call igetval( nx( 3 ), 'ngz  ' )
   close( unit=99 )

!  Not needed since ABINIT outputs single wavefxn file now, v9.8.3   
!   open( unit=99, file='masterwfile', form='formatted', status='unknown' )
!   rewind 99
!   read ( 99, * ) nwfile
!   close( unit=99 )
!   allocate( wfnam( nwfile ), ikl( 2, nwfile ) )
!   open( unit=99, file='listwfile', form='formatted', status='unknown' )
!   rewind 99
!   do i = 1, nwfile
!      read ( 99, * ) ikl( 1, i ), wfnam( i )
!      if ( i .gt. 1 ) ikl( 2, i - 1 ) = ikl( 1, i ) - 1
!   end do
!   close( unit=99 )
!   ikl( 2, nwfile ) = nk
   !
   open( unit=u1dat, file='u1.dat', form='unformatted', status='unknown' )
   open( unit=u1xdat, file='u1x.dat', form='unformatted', status='unknown' )
   rewind u1dat
   rewind u1xdat

   do i = 1, nk
       read(u1dat) ng,nspinor,nb
       write(6,*) " ik = ",i
       write(6,*) " "
       write(6,*) "npw =",ng,"nspinor=",nspinor,"nband=",nb
       allocate (g(ng,3))
       allocate (zr(ng,nb),zi(ng,nb))
       allocate(zlin(2*ng*nspinor))
       read(u1dat) g    
       write(6,*) "size gvec array =",size(g),"shape=",shape(g) 
       write(6,*) "shape zr, zi,zlin",shape(zr),shape(zi),shape(zlin)
       write(6,*) " "
       write(6,*) " "
       do iband = 1, nb
          read(u1dat) zlin
          zr(:,iband) = zlin(1:ng)
          zi(:,iband) = zlin(ng+1:2*ng*nspinor)
       enddo
       call gentoreal( nx, nb, zr, zi, ng, g, u1xdat, i .eq. 1 )
       deallocate( g, zr, zi,zlin )
   end do
   close( unit=u1dat )
   !
!   deallocate( wfnam, ikl )
   !
 end program conugtoux
