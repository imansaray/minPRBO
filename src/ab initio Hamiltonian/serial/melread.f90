subroutine melread_edit( melu, i, ilftl, ilfth, irgtl, irgth, pdota )
  implicit none
  !
  integer :: melu, i, ilftl, ilfth, irgtl, irgth
  complex( kind = kind( 1.0d0 ) ) :: pdota( ilftl : ilfth, irgtl : irgth )
  !
  if ( i .eq. 1 ) then
     open( unit=melu, file='tmels', form='formatted', status='unknown' )
     rewind melu
  end if
!  read( melu, * ) pdota( :, : )
  read( melu, '(8(1x,1e22.15))' ) pdota( :, : )
  !
  return
end subroutine melread_edit

subroutine melread( melu, i, ilftl, ilfth, irgtl, irgth, pdota )
  implicit none
  !
  integer :: melu, i, ilftl, ilfth, irgtl, irgth
!  complex( kind = kind( 1.0d0 ) ) :: pdota( 0 : 3, ilftl : ilfth, irgtl : irgth )
  complex( kind = kind( 1.0d0 ) ) :: pdota( ilftl : ilfth, irgtl : irgth )
  complex( kind = kind( 1.0d0 ) ) :: ptemp( ilftl : ilfth, irgtl : irgth ) ! to be able to read OCEAN ptmels.dat bin file
  !
  if ( i .eq. 1 ) then
!     open( unit=melu, file='tmels', form='formatted', status='unknown' )
     open( unit=melu, file='ptmels.dat', form='unformatted', status='old',access='stream' )
     rewind melu
  end if
  ptemp = 0.0
  read( melu ) ptemp
!  write(6,*) shape(ptemp)
!  write(6,*) " "
!  write(6,*) ptemp
!  write(6,*) " "
  pdota(:,:) = ptemp
  !pdota(0,:,:) = ptemp
  !pdota(1,:,:) = 0.0
  !pdota(2,:,:) = 0.0
  !pdota(3,:,:) = 0.0
!  read( melu, '(8(1x,1e22.15))' ) pdota( :, :, : )
  !
  return
end subroutine melread
