program bridgegw
  implicit none
  !
  integer :: i, j, nk, nb, bwflg
! integer :: nbc, nbv
  !
  real( kind = kind( 1.0d0 ) ) :: lumo, homo, egap, eshift, newgap, efermi, sef
  real( kind = kind( 1.0d0 ) ) :: egw, elda, vs, cs
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: e( :, : ), edft( :, : )
  !
  open( unit=99, file='gwipt', form='formatted', status='unknown' )
  rewind 99
  read( 99, * ) egw, elda, vs, cs
  write(6,*) " Real band gap =",egw,"DFT band gap =",elda,"T_V =",vs,"T_S =",cs
  write(6,*) " "
  close( unit=99 )
  !
  open( unit=99, file='efermiinrydberg.ipt', form='formatted', status='unknown' )
  rewind 99
  read( 99, * ) efermi
  close( unit=99 )
  efermi = efermi * 13.6057d0
  !
  write(6,*) " Fermi energy = ",efermi, " eV "
  write(6,*) " "

  open( unit=99, file='kandb.h', form='formatted', status='unknown' )
  call igetval( nk, 'nk   ' )
  call igetval( nb, 'nb   ' )
! call igetval( nbc, 'nbc  ' )
! call igetval( nbv, 'nbv  ' )

  bwflg = 0
  if (bwflg .eq. 1) nb = 2*nb

  allocate( e( nb, nk ), edft( nb, nk ) )
  !
  write(6,*) " Opening and reading in the cpbd file with band energies in eV"
  write(6,*) " "
  open( unit=99, file='cpbd', form='formatted', status='unknown' )
  rewind 99
  do i = 1, nk
     do j = 1, nb
        read ( 99,  *  ) e( j, i )
     end do
  end do
  close( unit=99 )
  edft = e
  !
  ! OLD 1 fwbw April 7 2006
  ! get homo
! homo = e( nbv, 1 )
! do i = 1, nk
!    do j = 1, nbv
!       if ( e( j, i ) .gt. homo ) homo = e( j, i )
!    end do
! end do
  !
  ! get lumo
! write ( 6, * ) 'homo=', homo
! lumo = e( nbv + 1, 1 )
! do i = 1, nk
!    do j = nbv + 1, nb
!       if ( e( j, i ) .lt. lumo ) lumo = e( j, i )
!    end do
! end do
! write ( 6, * ) 'lumo=', lumo
  !
  ! NEW 1 fwbw April 7 2006 begin: this change allows for forward plus backward pairs
  homo = e( 1, 1 )
  lumo = e( nb, 1 )
  do i = 1, nk
     do j = 1, nb
        if ( e( j, i ) .lt. efermi ) then
           homo = max( homo, e( j, i ) )
!CD        write ( 87, * ) efermi, j, e( j, i )
        else
           lumo = min( lumo, e( j, i ) )
!CD        write ( 88, * ) efermi, j, e( j, i )
        end if
     end do
  end do
  ! NEW 1 fwbw April 7 2006 end
  !
  !
  ! get gap
  egap = lumo - homo
  !
  ! adjust to vbm
! do i = 1, nk
!    e( 1 : nb, i ) = e( 1 : nb, i ) - homo
! end do
  e( :, : ) = e( :, : ) - homo
  sef = efermi - homo
  !
  ! determine shift
  eshift = egw - elda
  newgap = egap + eshift
  write(6,*) " sef = efermi - homo =",sef
  write ( 6, * ) 'eshift = ', eshift, " eV"
  write ( 6, * ) 'newgap = ', newgap, " eV"
  !     
  ! OLD 2 fwbw April 7 2006
! ! shift
! do i = 1, nk
!    do j = nbv + 1, nb
!       e( j, i ) = e( j, i ) + eshift
!    end do
! end do
  !
  ! stretch
! do i = 1, nk
!    do j = 1, nbv
!       e( j, i ) = ( 1.d0 + vs ) * e( j, i )
!    end do
! end do
! do i = 1, nk
!    do j = nbv + 1, nb
!       e( j, i ) = newgap + ( 1.d0 + cs ) * ( e( j, i ) - newgap )
!    end do
! end do
  ! NEW 2 fwbw April 7 2006 begin: do it based on fermi energy criterion for fw+bw
  ! shift and stretch
  do i = 1, nk
     do j = 1, nb
        if ( e( j, i ) .lt. sef ) then
           e( j, i ) = ( 1.d0 + vs ) * e( j, i )
           write ( 97, * ) i, j
        else
           e( j, i ) = e( j, i ) + eshift
           e( j, i ) = newgap + ( 1.d0 + cs ) * ( e( j, i ) - newgap )
           write ( 98, * ) i, j
        end if
     end do
  end do
  ! NEW 2 fwbw April 7 2006 end
  !
  !
  open( unit=36, file='ebdat', form='formatted', status='unknown' )
  open( unit=40, file='bs.dat', form='formatted', status='unknown' )
  rewind 36
  rewind 40
  do j = 1, nb
     do i = 1, nk
!        write ( 36, '(2(1x,2f12.6))' ) e( j, i ), edft( j, i )
        write ( 36, * ) e( j, i ), edft( j, i )
        write ( 40, * ) i, e( j, i )
     end do
  end do
  close( unit=40 )
  close( unit=36 )
  !
  write ( 6, * ) 'gap before moving:', egap
  write ( 6, * ) 'gap after stretching:', newgap
  !
  stop
end program bridgegw
