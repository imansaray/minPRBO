subroutine smkxc( c, u, nk, ng, ngx, ngy, ngz, nb, nbc, nbv )
  implicit none
  !
  integer :: nk, ng, ngx, ngy, ngz, nb, nbc, nbv
  complex( kind = kind( 1.0d0 ) ) :: u( nb, nk, ng ), c( nbv, nbc, nk )
  !
  integer :: ix, ik, j, jj, i, igx, igy, igz
  complex( kind = kind( 1.0d0 ) ) :: rm1, s1, s2
  real( kind = kind( 1.0d0 ) ) :: pi, omega, dum
  real( kind = kind( 1.0d0 ) ), allocatable :: dnr( :, :, : ), dni( :, :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: kxc( :, :, : )
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  pi = 4.0d0 * atan( 1.0d0 )
  !
  allocate( dnr( ngz, ngy, ngx ), dni( ngz, ngy, ngx ), kxc( ngz, ngy, ngx ) )
  !
  dnr = 0
  dni = 0
  ix = 0
  do igx = 1, ngx
     do igy = 1, ngy
        do igz = 1, ngz
           ix = ix + 1
           s2 = 0
           do ik = 1, nk
              do j = 1, nbc
                 jj = nbv + j
                 s1 = 0
                 do i = 1, nbv
                    s1 = s1 + conjg( u( i, ik, ix ) ) * c( i, j, ik )
                 end do
                 s2 = s2 + s1 * u( jj, ik, ix )
              end do
           end do
           dnr( igz, igy, igx ) = s2 * ng
           dni( igz, igy, igx ) = -s2 * ng * rm1
        end do
     end do
  end do
  !
  open( unit=99, file='rhokxc.txt', form='formatted', status='unknown' )
  rewind 99
  do igx = 1, ngx
     do igy = 1, ngy
        do igz = 1, ngz
           read ( 99, * ) dum, kxc( igz, igy, igx )
        end do
     end do
  end do
  close( unit=99 )
  open( unit=99, file='omega.h', form='formatted', status='unknown' )
  rewind 99
  call rgetval( omega, 'volbo' )
  close( unit=99 )
  dnr( :, :, : ) = dnr( :, :, : ) * kxc( :, :, : ) * 27.2114d0 / ( dble( nk ) * omega )
  dni( :, :, : ) = dni( :, :, : ) * kxc( :, :, : ) * 27.2114d0 / ( dble( nk ) * omega )
  !
  c = 0
  ix = 0
  do igx = 1, ngx
     do igy = 1, ngy
        do igz = 1, ngz
           ix = ix + 1
           s1 = dnr( igz, igy, igx ) + rm1 * dni( igz, igy, igx )
           do ik = 1, nk
              do j = 1, nbc
                 jj = nbv + j
                 s2 = s1 * conjg( u( jj, ik, ix ) )
                 do i = 1, nbv
                    c( i, j, ik ) = c( i, j, ik ) + s2 * u( i, ik, ix )
                 end do
              end do
           end do
        end do
     end do
  end do
  !
  deallocate( dnr, dni, kxc )
  !
  return
end subroutine smkxc
