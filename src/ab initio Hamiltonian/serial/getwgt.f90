subroutine getwgt( nkf, nbv, nbc, psi )
  implicit none
  !
  integer, parameter :: pdadat = 34
  !
  integer nkf, nbv, nbc
  double complex psi( nbv, nbc, nkf )
  !
  integer ik, ic, i, j, nptot, ntot, riter, nkpts, numbands, atno, ncore, lc, nc, mc, k
  double complex pol( 0 : 3 )
  double precision sum
  logical dipok
  character * 3 meltyp
  !
  complex( kind = kind( 1.0d0 ) ), allocatable :: wgt( :,:, : ) ! no polariz array
  complex(kind=kind(1.0d0)) :: subtot, rm1
  real(kind=kind(1.0d0)), allocatable :: pcoefr(:,:), pcoefi(:,:), rex(:,:)
  real(kind=kind(1.0d0)) :: tau(3), br, bi
  character * 10 :: str
  !
  integer :: nedges, edge_iter, icml, icms, ivms
  double complex, allocatable :: tmp_psi( :, :, : )
  character*10 :: echamp_file
  rm1 = -1.0d0
  rm1 = sqrt(rm1)
  write ( 6, * ) 'entering getwgt'
  write(6,*) " "
  write(6,*) " Enter a case: dip, cor, or default (nothing)"
  write(6,*) " "

  allocate( wgt( 0 : 3, nbv, nbc ) )
!  allocate( wgt( nbv, nbc ) )
  !
  ! HERE we are with the polarization input via file ...
  !
  read ( 5, '(1a3)' ) meltyp
!  write(6,*) " Only doing the default valence case, with no separate polarization vector"
  !
  pol = 0
  !
  !
  select case( meltyp )
  case( 'dip')
!  if ( meltyp .eq. 'dip' ) then
     open( unit=99, file='dipok', form='formatted', status='unknown' )
     rewind 99
     read ( 99, * ) dipok
     close( unit=99 )
     if ( .not. dipok ) stop 'large Q: dipole calc forbidden'
     open( unit=99, file='polariz.h', form='formatted',status='unknown' )
     rewind 99
     do i = 1, 3
        read ( 99, * ) sum
        pol( i ) = sum
     end do
     close( unit=99 )
     sum = 0.0d0
     do i = 1, 3
        sum = sum + dble( pol( i ) * conjg( pol( i ) ) )
     end do
     sum = 1.0d0 / sqrt( sum )
     do i = 1, 3
        pol( i ) = sum * pol( i )
     end do
!  else
!     pol( 0 ) = 1
!  end if
!10 continue
!     allocate( wgt( 0 : 3, nbv, nbc ) )
  write ( 6, '(1a3)' ) meltyp
  write ( 6, '(4(3x,2f7.4))' ) pol
  !
  open( unit=pdadat, file='pdadat', status='unknown' )
  rewind pdadat
  do ik = 1, nkf
     do ic = 0, 3
        do i = 1, nbv
           do j = 1, nbc
              read ( pdadat, '(2(1x,1e22.15))' ) wgt( ic, i, j )
           end do
        end do
     end do
     do j = 1, nbc
        do i = 1, nbv
           psi( i, j, ik ) = dot_product( wgt(:,i,j),pol(:) )
        end do
     end do
  end do
  close( unit=pdadat )
  !
  deallocate( wgt )
  !
  case( 'cor')
    write(6,*) 'cor'
    ! First setp, read in |y(w)>, which is to say the exciton wavefunction amplitudes
    open( unit=99, file='wvfcninfo', form='unformatted', status='old')
    rewind( 99 )
    read( 99 ) numbands, nkpts
    close( 99 )
    !
    open( unit=99, file='ZNL', form='formatted', status='old')
    read( 99, * ) atno, ncore, lc
    close( 99 )
    !
    nc = 4 * ( 2 * lc + 1 )
    if ( numbands .gt. nbc ) stop 'Not enough bands...'
    if ( nkpts .ne. nkf ) stop 'K-point counts disagree'
    !
    allocate(rex( numbands*nkpts*nc, 2) )


    ! If there is more than one site per cell then we will need to 
    !  coherently sum the contributions from each of them.
    ! For now we will do this by assuming nedges is the correct input
    open( unit=99, file='nedges', form='formatted', status='old' )
    read( 99, * ) nedges
    close( 99 )
    allocate( tmp_psi( nbv, nbc, nkf ) )
    allocate(pcoefr(nkf,nbv), pcoefi(nkf,nbv) )

    psi(:,:,:) = 0.0d0
    do edge_iter = 1, nedges
      write(6,*) edge_iter
 
      tmp_psi(:,:,:) = 0.d0
      write( echamp_file, '(A7,I1)') 'echamp.', edge_iter
      write(6,*) echamp_file
!      open( unit=99, file='echamp', form='unformatted', status='old')
      open( unit=99, file=echamp_file, form='unformatted', status='old')
      rewind( 99 )
      read( 99 ) rex
      close( 99 )
    
!      riter = 1
!      do ic = 1, nc
!        do ik = 1, nkf
!          do j = 1, nbc
!            do i = 1, nbv
!              tmp_psi(i, j, ik ) = tmp_psi(i, j, ik ) + rex(riter,1) + rm1 * rex(riter,2)
!            enddo
!            riter = riter + 1 
!          enddo
!        enddo
!      enddo
      ic = 0
      do icms = -1, 1, 2
        do icml = -lc, lc
          do ivms = -1, 1, 2
            ic = ic + 1
            if( .true. ) then
!            if ( icms .eq. ivms ) then
              ! Not doing spin stuff
              write( str, '(1a4,1i2.2,a1,i1)' ) 'beff', 1 + icml + lc, '.', edge_iter
              write(6,*) str
              open( unit=98, file=str, form='unformatted', status='old' )
              rewind( 98 )
              read( 98 ) tau( : )
              do ik = 1, nkf
                do i = 1, nbv
                  read( 98 ) br, bi
                  do j = 1, nbc
                    riter = ( ic - 1 ) * nkf * nbc + ( ik - 1 ) * nbc + j
! beff gives the XES matrix elements which should be conjugated
                    tmp_psi(i, j, ik ) = tmp_psi(i, j, ik ) &
                                       + rex(riter,1) * br + rex(riter,2) * bi  &
                                       - rm1 * rex(riter,1) * bi + rm1 * rex(riter,2) * br
                  end do
                end do
              end do
              close( 98 )
            end if
          end do
        end do
     end do

!    do ik = 1, nkf
!      do j = 1, nbc
!        subtot = 0.0d0
!        do i=1,nc
!          subtot = subtot + rex(riter,1) + rm1 * rex(riter,2)
!          riter = riter + 1
!        enddo
!        do i = 1, nbv
!          psi(j, i, ik) = subtot
!        enddo
!      enddo
!    enddo
    !
      ! Step two is read in the transition matrix elements between the valence and core hole
      open( unit=99, file='wvfvainfo', form='unformatted', status='old')
      rewind( 99 )
      read( 99 ) numbands, nkpts
      close( 99 )
      if ( numbands .ne. nbv) stop 'Valence band mismatch'
      if ( nkpts .ne. nkf ) stop 'K-point counts disagree for valence'
      !
!      pcoefr(:,:) = 0.0d0
!      pcoefi(:,:) = 0.0d0
!      do mc = -lc, lc
!        write( str, '(1a4,1i2.2,a1,i1)' ) 'beff', 1 + mc + lc, '.', edge_iter
!        open( unit=98, file=str, form='unformatted', status='old' )
!        rewind( 98 )
!        read( 98 ) tau( : )
!        do ik = 1, nkf
!          do i=1, nbv
!            read( 98 ) br, bi
!            pcoefr( ik, i ) = pcoefr( ik, i ) + br
!            pcoefi( ik, i ) = pcoefi( ik, i ) + bi
!          enddo
!        enddo
!        close( 98 )
!      enddo
!      do j= 1, nbc
!        do ik = 1, nkf
!          do i = 1, nbv
!! switched - <-> j
!!!!!JTV!            tmp_psi(i, j, ik) = tmp_psi(i, j, ik ) * CMPLX(pcoefr( ik, i ), pcoefi( ik, i ), kind(1.0d0) ) 
!          enddo
!        enddo
!      enddo
      psi( :, :, : ) = psi( :, :, : ) + tmp_psi( :, :, : )
    enddo ! edge_iter
    deallocate(pcoefr, pcoefi)    
    deallocate( rex )
!    open(unit=99,file='vacoremels',form='unformatted',status='old')
!    read(99) nptot, ntot
!    read(99) tau(:)
!    if (ntot .ne. nkf*nbv) stop 'nbv not the same as nbuse.ipt'
!    allocate(pcoefr(nptot,nkf,nbv), pcoefi(nptot,nkf,nbv) )
!    read(99) pcoefr
!    read(99) pcoefi
!    close(99)
!    do ik = 1, nkf
!      do i = 1, nbv
!        subtot = 1.d0
!        do j = 1, nptot
!          subtot = pcoefr(j,ik,i) + rm1 * pcoefi(j,ik,i)
!        enddo
!        do j = 1, nbc
!          psi(j, i, ik) = psi(j, i, ik) * subtot
!        enddo
!      enddo
!    enddo
!    deallocate(pcoefr, pcoefi)
!
!
  case default

  write(6,*) " Doing the default case : just regular valence calc"
 ! assume we still want normal valence behavior
    pol(0) = 1
!    goto 10
!     allocate( wgt( 0 : 3, nbv, nbc ) )
  write ( 6, '(1a3)' ) meltyp
  write ( 6, '(4(3x,2f7.4))' ) pol
  !
  open( unit=pdadat, file='pdadat', status='unknown' )
  rewind pdadat
  do ik = 1, nkf
     do ic = 0, 3
        do i = 1, nbv
           do j = 1, nbc
              read ( pdadat, '(2(1x,1e22.15))' ) wgt( ic, i, j )
!              read ( pdadat, '(2(1x,1e22.15))' ) wgt( i, j )
           end do
        end do
     end do
     do j = 1, nbc
        do i = 1, nbv
           psi( i, j, ik ) = dot_product( wgt(:,i,j),pol(:) )
!           psi( i, j, ik ) = wgt(i,j)
        end do
     end do
  end do
  close( unit=pdadat )
  !
  deallocate( wgt )

  end select
  write ( 6, * ) 'done with getwgt'
  open( unit=99,file='psi_test',form='formatted')
  do i=1,nkf
    do j=1,nbc
      do k=1,nbv
        write(99,'(E12.5,X)') psi(k,j,i)
      enddo
    enddo
  enddo
  close( 99 )

  !
  return
end subroutine getwgt

subroutine getwgt_new( nkf, nbv, nbc, ilftl, ilfth, irgtl, irgth, n, psi )
  implicit none
  !
  integer, parameter :: pdadat = 34
  !
  integer nkf, nbv, nbc, n,ilftl,ilfth,irgtl,irgth,il,ir
  double complex psi( nbv, nbc, nkf ), ptemp(nbv,nbc)
  double precision :: rptemp(nbv,nbc), iptemp(nbv,nbc)
  !
  integer ik, ic, i, j
  double complex pol( 0 : 3 )
  double precision sum
  logical dipok
  character * 3 meltyp
  !
!  complex( kind = kind( 1.0d0 ) ), allocatable :: wgt( :,:, : ) ! no polariz array
  complex( kind = kind( 1.0d0 ) ), allocatable :: wgt( :, : ) ! to follow OCEAN format
  complex(kind=kind(1.0d0)) :: subtot, rm1

  rm1 = -1.0d0
  rm1 = sqrt(rm1)
  write ( 6, * ) 'entering getwgt'
  write(6,*) " "
!  write(6,*) " Enter a case: dip, cor, or default (nothing)"
  write(6,*) " "

!  allocate( wgt( 0 : 3, nbv, nbc ) )
  allocate( wgt( nbv, nbc ) )
  !
  ! HERE we are with the polarization input via file ...
  !
!  write(6,*) " Enter polarization type: 'dip' or default/press enter"
!  read ( 5, '(1a3)' ) meltyp
  write(6,*) " Polarization type: 'dip' or default"
  meltyp = ' '

!  write(6,*) " Only doing the default valence case, with no separate polarization vector"
  !
  pol = 0
  !
  !
  if ( meltyp .eq. 'dip' ) then
     open( unit=99, file='dipok', form='formatted', status='unknown' )
     rewind 99
     read ( 99, * ) dipok
     close( unit=99 )
     if ( .not. dipok ) stop 'large Q: dipole calc forbidden'
     open( unit=99, file='polariz.h', form='formatted',status='unknown' )
     rewind 99
     do i = 1, 3
        read ( 99, * ) sum
        pol( i ) = sum
     end do
     close( unit=99 )
     sum = 0.0d0
     do i = 1, 3
        sum = sum + dble( pol( i ) * conjg( pol( i ) ) )
     end do
     sum = 1.0d0 / sqrt( sum )
     do i = 1, 3
        pol( i ) = sum * pol( i )
     end do
  else
     pol( 0 ) = 1
  end if
!10 continue
!     allocate( wgt( 0 : 3, nbv, nbc ) )
  write ( 6, '(1a3)' ) meltyp
  write ( 6, '(4(3x,2f7.4))' ) pol

!  ! Directly reading ocean data in ptmels.dat
!  open( unit=150, file='ptmels.dat', form='unformatted', status='old',access='stream')
!  rewind 150
!
!  psi=0.0
!  do ik=1,nkf
!     ptemp = 0.0
!     read(150) ptemp
!     do j = 1, nbc
!        do i = 1, nbv
!           psi(i,j,ik) = ptemp(i,j)
!        end do
!     end do
!  end do
!  close(150)

  ! Directly reading ocean data in tmels
  open( unit=150, file='tmels', form='formatted', status='unknown')
  rewind 150

  psi=0.0
  do ik=1,nkf
     ptemp = 0.0
     do j = 1, nbc
        do i = 1, nbv
           read(150,*) rptemp(i,j),iptemp(i,j)
        end do
     end do
     do i = 1, nbv
        do j = 1, nbc
           psi(i,j,ik) = cmplx(rptemp(i,j),-iptemp(i,j))
        end do
     end do
  end do
  close(150)
  write(6,*) 'done with getwgt'

  return
end subroutine getwgt_new
