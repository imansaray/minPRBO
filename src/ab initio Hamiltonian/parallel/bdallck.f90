! this now allows backwards + forwards pairs
subroutine bdallck( nbv, nbc, nk, allow, sfact, ebd, eblda, qpflg, bwflg )
  implicit none
  !
  integer :: nk, qpflg, bwflg, nbv, nbc
  real( kind = kind( 1.0d0 ) ) :: allow( nbv, nbc, nk )
  real( kind = kind( 1.0d0 ) ) :: sfact( nbv, nbc, nk )
  real( kind = kind( 1.0d0 ) ) :: ebd( nbv + nbc, nk )
  real( kind = kind( 1.0d0 ) ) :: eblda( nbv + nbc, nk )
  !
  integer :: i, j, jj, k
  real( kind = kind( 1.0d0 ) ) :: ev, ec, clipl, cliph, efermi
  !
  ! get e( k ) on the regular mesh of k-points, allow an electron-hole pair
  ! for ev < clipl < ec < cliph, based on LDA band topology
  allow( :, :, : ) = 0.0d0
  sfact( :, :, : ) = 0.0d0
  open( unit=99, file='ldaclips', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) clipl, cliph
  close( unit=99 )
  open( unit=99, file='efermiinrydberg.ipt', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) efermi
  close( unit=99 )
  efermi = efermi * 13.6057d0
  !
  write(6,*) " "
  write(6,*) " qpflg = ",qpflg, "n = ",nk*nbv*nbc
  write(6,*) " "
  call getband( nbv+nbc, nk, ebd, eblda, qpflg )
  do k = 1, nk
     do j = 1, nbc
        jj = j + nbv
        ec = eblda( jj, k ) !, original
!        ec = ebd( jj, k )
        do i = 1, nbv
           ev = eblda( i, k ) !, original
!           ev = ebd( i, k )
           if ( max( ev, ec ) .lt. cliph ) then
              if ( ( ev .lt. efermi ) .and. ( ec .gt. efermi ) ) then
                 allow( i, j, k ) = 1.0d0
                 sfact( i, j, k ) = 1.0d0
              end if
              if ( ( ev .gt. efermi ) .and. ( ec .lt. efermi ) .and. ( bwflg .eq. 1 ) ) then
                 allow( i, j, k ) = 1.0d0
                 sfact( i, j, k ) = -1.0d0
              end if
           end if
        end do
     end do
  end do
  !
  return
end subroutine bdallck
