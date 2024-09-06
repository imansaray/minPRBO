subroutine sizereport( siz, var )
  implicit none
  !
  integer :: siz
  character * 10 :: var
  !
  integer :: runtot, oldsiz, maxtot, i, maxsiz
  character * 16 :: fnam
  !
  write ( fnam, '(1a6,1a10)' ) '.memu_', var
  !
  ! try to read old memory use ... an error sets it to zero
  open( unit=99, file=fnam, form='formatted', iostat=i, status='old' )
  rewind 99
  if ( i .eq. 0 ) then
     read ( 99, * ) oldsiz, maxsiz
  else
     oldsiz = 0
     maxsiz = 0
  end if
  close( unit=99 )
  !
  ! try to read cumulative memory use ... error implies no previous use
  open( unit=99, file='cumlog', form='formatted', iostat=i, status='old' )
  rewind 99
  if ( i .eq. 0 ) then
     read ( 99, * ) runtot
     read ( 99, * ) maxtot
  else
     runtot = 0
     maxtot = 0
  end if
  close( unit=99 )
  !
  ! output cumulative usage
  open( unit=99, file='cumlog', form='formatted', status='unknown' )
  rewind 99
  runtot = runtot + siz - oldsiz
  write ( 99, '(1i12,5x,1a6)' ) runtot, 'runtot'
  write ( 99, '(1i12,5x,1a6)' ) max( maxtot, runtot ), 'maxtot'
  close( unit=99 )
  !
  ! output particular usage
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(2i12,5x,1a10)' ) siz, max( maxsiz, siz ), var
  close( unit=99 )
  !
  return
end subroutine sizereport
