subroutine enkread( enku, i, val, il, ih, w )
  implicit none
  !
  integer :: enku, i, il, ih
  real( kind = kind( 1.0d0 ) ) :: w( il : ih )
  logical :: val
  !
!  write(6,*) " Writing out the band energies in Rydberg from enkfile "
  open( unit=enku, file='enkfile', form='formatted', status='unknown' )
  if ( ( i .eq. 1 ) .and. ( val ) ) then
!     open( unit=enku, file='enkfile', form='formatted', status='unknown' )
     rewind enku
  end if
  read ( enku, * ) w( : )
!  write(6,*) " "
!  write(6,*) w
  !
  return
end subroutine enkread
