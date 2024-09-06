subroutine adump( chan, nam, pnorm, nk, n, a )
  implicit none
  !
  integer :: chan, nk, n
  character * 5 :: nam
  double precision :: pnorm 
 ! complex( kind = kind( 1.0d0 ) ):: a( 0 : n )
  real( kind = kind( 1.0d0 ) ):: a( 0 : n )
  !
  integer i
  !
  open( unit=chan, file=nam, form='formatted', status='unknown' )
  rewind chan
  write ( chan, * ) pnorm, nk
  do i = 0, n
     write ( chan, * ) i, a( i )
  end do
  close( unit=chan )
  return
end subroutine adump

subroutine cadump( chan, nam, pnorm, nk, n, ca )
  implicit none
  !
  integer :: chan, nk, n
  character * 5 :: nam
  double precision :: pnorm 
  complex( kind = kind( 1.0d0 ) ):: ca( 0 : n )
  real( kind = kind( 1.0d0 ) ):: a( 0 : n )
  !
  integer i
  !
  open( unit=chan, file=nam, form='formatted', status='unknown' )
  rewind chan
  write ( chan, * ) pnorm, nk
  do i = 0, n
     write ( chan, * ) i, ca( i )
  end do
  close( unit=chan )
  return
end subroutine cadump

subroutine fdump( chan, nam, pnorm, nk, n, a )
  implicit none
  !
  integer :: chan, nk, n
  character * 5 :: nam
  double precision :: pnorm 
 ! complex( kind = kind( 1.0d0 ) ):: a( 0 : n )
  real( kind = kind( 1.0d0 ) ):: a( 0 : n )
  !
  integer i
  !
  open( unit=chan, file=nam, form='formatted', status='unknown' )
  rewind chan
!  write ( chan, * ) pnorm, nk
  do i = 1, n
     write(6,*)"i =",i,"a(i) =",a(i)
     write(6,*)" "
     write ( chan, * ) i, a( i )
  end do
  close( unit=chan )
  return
end subroutine fdump
