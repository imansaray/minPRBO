program setup1
  implicit none
  ! changed qtransf to qinb to match the OCEAN formulas
  integer, parameter :: kptsdat = 29, cpbd = 48, pdadat = 72
  integer, parameter :: enku= 73, melu = 74
  real( kind = kind( 1.0d0 ) ), parameter :: eryd = 13.6057d0, a0 = 0.529177d0
  character * 7, parameter :: u7 = 'unknown'
  character * 9, parameter :: f9 = 'formatted'
  !
  integer :: i, j, il, ir, ik, ix, iy, iz, band
  integer :: ilftl, ilfth, irgtl, irgth, nk, nkx, nky, nkz
  integer :: nb, nbl, nbh, nband, nbv, nbc, ngx, ngy, ngz, nx
  real( kind = kind( 1.0d0 ) ) :: minccap, cx, cy, cz, efermi, pi
  real( kind = kind( 1.0d0 ) ) :: bmet( 3, 3 ), bvec( 3, 3 ), avec( 3, 3 )
  real( kind = kind( 1.0d0 ) ) :: celvol
  real( kind = kind( 1.0d0 ) ) :: qinb( 3 ), tmp( 3 ), kxv( 3 ),inv_qlength
  logical :: new
  character * 80 :: fstr
  !
  integer, allocatable :: occ( : , : )
  real( kind = kind( 1.0d0 ) ), allocatable :: wv( : ), wc( : ), wtmp(:)
!  complex( kind = kind( 1.0d0 ) ), allocatable :: pdota( :, :, : ) !, using the approx polarization operator
  complex( kind = kind( 1.0d0 ) ), allocatable :: pdota( :, : )
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  open( unit=99, file='efermiinrydberg.ipt', form=f9, status=u7 )
  rewind 99
  read ( 99, * ) efermi
  close( unit=99 )
  !
  open( unit=99, file='avecsinbohr.ipt', form=f9, status=u7 )
  rewind 99
  read ( 99, * ) avec( :, : )
  close( unit=99 )
  call getrecvec( avec, bvec )
  do i = 1, 3
     do j = 1, 3
        bmet( i, j ) = dot_product( bvec( :, i ), bvec( :, j ) )
     end do
  end do
  !
  open( unit=99, file='brange.ipt', form=f9, status=u7 )
  rewind 99
  read ( 99, * ) ilftl, ilfth, irgtl, irgth
  close( unit=99 )
  open( unit=99, file='kmesh.ipt', form=f9, status=u7 )
  rewind 99
  read ( 99, * ) nkx, nky, nkz
  close( unit=99 )
  open( unit=99, file='xmesh.ipt', form=f9, status=u7 )
  rewind 99
  read ( 99, * ) ngx, ngy, ngz
  close( unit=99 )
  open( unit=99, file='qinunitsofbvectors.ipt', form=f9, status=u7 )
  rewind 99
  read ( 99, * ) qinb( : )
  close( unit=99 )
  open( unit=99, file='bvecs', form=f9, status=u7 )
  rewind 99
  read ( 99, * ) bvec( :, : )
  close( unit=99 )

  nk = nkx * nky * nkz
  nbl = ilftl
  nbh = irgth
  nbv = 1 + ilfth - ilftl
  nbc = 1 + irgth - irgtl
!  nb = nbv + nbc ! only if there's no overlap/coupling. Need to include the bwflg here
  !nband = nb
  nb = irgth
  nband = irgth ! more general, or get this from the DFT results.
  
  allocate( occ( nk, nband ) )
  occ( :, : ) = 0
!  allocate( pdota( 0 : 3, ilftl : ilfth, irgtl : irgth ) )
  allocate( pdota( ilftl : ilfth, irgtl : irgth ) )
  open( unit=99, file='kandb.h', form=f9, status=u7 )
  rewind 99
  call iputval( ilftl, 'ilftl' )
  call iputval( ilfth, 'ilfth' )
  call iputval( irgtl, 'irgtl' )
  call iputval( irgth, 'irgth' )
  call iputval( nkx, 'nkx  ' )
  call iputval( nky, 'nky  ' )
  call iputval( nkz, 'nkz  ' )
  call iputval( nk, 'nk   ' )
  call iputval( nb, 'nb   ' )
  call iputval( nbl, 'nbl  ' )
  call iputval( nbh, 'nbh  ' )
  call iputval( nband, 'nband' )
  call iputval( nbv, 'nbv  ' )
  call iputval( nbc, 'nbc  ' )
  nx = ngx * ngy * ngz
  call iputval( ngx, 'ngx  ' )
  call iputval( ngy, 'ngy  ' )
  call iputval( ngz, 'ngz  ' )
  call iputval( nx, 'ng   ' )
  call iputval( nx, 'nx   ' )
  close( unit=99 )
  !
  cx = 2.0d0 * pi / dble( nkx )
  cy = 2.0d0 * pi / dble( nky )
  cz = 2.0d0 * pi / dble( nkz )
  open( unit=kptsdat, file='kpts.dat', form=f9, status=u7 )
  rewind kptsdat
  ik = 0
  do ix = 1, nkx
     do iy = 1, nky
        do iz = 1, nkz
           ik = ik + 1
           kxv( 1 ) = cx * dble( ix - 1 )
           kxv( 2 ) = cy * dble( iy - 1 )
           kxv( 3 ) = cz * dble( iz - 1 )
           write( kptsdat, '(1i4,3f10.6)' ) ik, kxv( : )
        end do
     end do
  end do
  close( unit=kptsdat )
  !
  open( unit=99, file='vecy', form=f9, status='unknown')
  rewind 99
  fstr = '(3f15.10,2i5)'
  do i = 1, 3
     do j = 1, 3
        write ( 99, fstr ) bmet( i, j ), bvec( i, j ), avec( i, j ), i, j
     end do
  end do
  write( 99, * ) 'each lines reads:'
  write( 99, * ) 'b_i dot b_j, i component of b_j, i component of a_j, i, j'
  write( 99, * ) 'b-vectors are in inverse bohrs; a-vectors are in bohrs'
  close( unit=99 )
  !
  call cross( avec( :, 1 ), avec( :, 2 ), tmp )
  celvol = abs( dot_product( tmp, avec( :, 3 ) ) )
  open( unit=99, file='omega.h', form=f9, status=u7 )
  rewind 99
  call rputval( celvol, 'volbo' )
  celvol = celvol * a0 ** 3
  call rputval( celvol, 'volan' )
  call rputval( qinb( 1 ), 'qtra1' )
  call rputval( qinb( 2 ), 'qtra2' )
  call rputval( qinb( 3 ), 'qtra3' )
  close( unit=99 )
  !
   
  ! Calculating the inverse q length to normalize ptmels.dat 
  inv_qlength = (qinb(1)*bvec(1,1) +qinb(2)*bvec(1,2)+qinb(3)*bvec(1,3))**2 &
               + (qinb(1)*bvec(2,1) +qinb(2)*bvec(2,2)+qinb(3)*bvec(2,3))**2 &
               + (qinb(1)*bvec(3,1) +qinb(2)*bvec(3,2)+qinb(3)*bvec(3,3))**2 

  inv_qlength = 1.0/dsqrt(inv_qlength)

   write(6,*) " "
   write(6,*) " inverse q length = ",inv_qlength
   write(6,*) " "

  allocate( wtmp(ilftl:irgth),wv( ilftl : ilfth ), wc( irgtl : irgth ) )
  !
  open( unit=cpbd, file='cpbd', status=u7 ) ! band energies in eV
  open( unit=pdadat, file='pdadat', status=u7 ) ! transition matrix elments
  rewind cpbd
  rewind pdadat
  do i = 1, nk
     write(6,*) i
     !
     ! read and change format of band energies
     call enkread( enku, i, .false., ilftl, irgth, wtmp )
     wv = wtmp( ilftl : ilfth )
     wc = wtmp( irgtl : irgth )
!     call enkread( enku, i, .false., irgtl, irgth, wc )
     write ( cpbd, '(1f14.8)' ) eryd * wv( ilftl : ilfth )
     write ( cpbd, '(1f14.8)' ) eryd * wc( irgtl : irgth )
     !
     ! determine occupancies
!     do j = 1, nbv
!        band = j - 1 + ilftl
!        if ( wv( band ) .lt. efermi ) occ( i, j ) = 1
!     end do
!     do j = 1, nbc
!        band = j - 1 + irgtl
!        if ( wc( band ) .lt. efermi ) occ( i, j + nbv ) = 1
!     end do
     !
     ! keep track of highest band energy sampled at all k-points
     if ( i .eq. 1 ) then
        minccap = wc( irgth )
        new = .true.
     else
        new = ( wc( irgth ) .lt. minccap )
        minccap = min( minccap, wc( irgth ) )
     end if
     if ( new ) write ( 6, '(2x,2i5,1f10.5)' ) i, irgth, wc( irgth )
    ! 
    ! read and change format of matrix elements
     call melread( melu, i, ilftl, ilfth, irgtl, irgth, pdota )

     ! Normalizing pdota by inv_qlength
!     pdota = inv_qlength*pdota
!     write(6,*) " "
!     write(6,*) "Finished normalizing pdota. Now writing to file pdadat."
!     write(6,*) " "
!     do j = 0, 3
        do il= ilftl, ilfth
           do ir = irgtl, irgth
!              write ( pdadat, '(2(1x,1e22.15))' ) pdota( j, il, ir )
              write ( pdadat, '(2(1x,1e22.15))' ) pdota( il, ir )
           end do
        end do
!     end do
  end do  !nk loop
  close( unit=cpbd )
  close( unit=pdadat )
  !
  open( unit=99, file='ldaclips', form=f9, status=u7 )
  rewind 99
  write ( 99, '(2x,2f20.15)' ) eryd * efermi, eryd * minccap
  close( unit=99 )
  !
  open( unit=99, file='occ.dat', form='unformatted', status=u7 )
  rewind 99
  write ( 99 ) occ
  close( unit=99 )
  !
  deallocate( wc, wv, pdota, occ )
  !
end program setup1
