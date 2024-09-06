subroutine getband( nb, nk, e, elda, qpflg )
  implicit none
  !
  integer :: nb, nk, qpflg
  double precision :: e( nb, nk ), elda( nb, nk )
  !
  integer :: ib, ik
  !
  open( unit=99, file='ebdat', status='unknown' )
  rewind 99
  write(6,*) " nb = ",nb, "nk =",nk
  write(6,*) " "
  do ib = 1, nb
     do ik = 1, nk
        read ( 99, * ) e( ib, ik ), elda( ib, ik )
!        read ( 99, '(2(1x,2f12.6))' ) e( ib, ik ), elda( ib, ik )
     end do
  end do
  close( unit=99 )
  if ( qpflg .eq. 0 ) e( :, : ) = elda( :, : )
  !
  return
end subroutine getband
