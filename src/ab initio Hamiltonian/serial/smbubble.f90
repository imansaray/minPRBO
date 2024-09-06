subroutine smbubble( c, u, nk, ng, ngx, ngy, ngz, nb, nbc, nbv )
  implicit none
  !
  integer :: nk, ng, ngx, ngy, ngz, nb, nbc, nbv
  complex( kind = kind( 1.0d0 ) ) :: u( nb, nk, ng ), c( nbv, nbc, nk )
  !
  integer :: ix, ik, j, jj, i, idwrk, igx, igy, igz, mode, nn1
  complex( kind = kind( 1.0d0 ) ) :: rm1, s1, s2
  real( kind = kind( 1.0d0 ) ), allocatable :: dnr( :, :, : ), dni( :, :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: wrk( : )
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  !
  allocate( dnr( ngz, ngy, ngx ), dni( ngz, ngy, ngx ) )
  idwrk = 2 * max( ngx ** 2 + ngx, ngy ** 2 + ngy, ngz ** 2 + ngz )
  allocate( wrk( idwrk ) )
  !
  dnr = 0
  dni = 0
  ix = 0
!  print *, "u( :, 1, 1 ) = ",u( :,1,1 )
!  print *
!  print *, "c( :, 1, 1 ) = ",c( :,1,1 )
!  print *
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
  mode = -1
  nn1 = ngz
  call cfft( dnr, dni, nn1, ngz, ngy, ngx, mode, wrk, idwrk )
  call coulmult( ngz, ngy, ngx, nk, dnr, dni )
  mode = +1
  call cfft( dnr, dni, nn1, ngz, ngy, ngx, mode, wrk, idwrk )
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
  deallocate( dnr, dni, wrk )
  !
  return
end subroutine smbubble
