subroutine ctrdfft( ar, ai, nxl, nxh, nyl, nyh, nzl, nzh, mode, wrk, idwrk )
  implicit none
  !
  integer :: nxl, nxh, nyl, nyh, nzl, nzh, mode, idwrk
  !
  real( kind = kind( 1.0d0 ) ) :: ar( nxl : nxh, nyl : nyh, nzl : nzh )
  real( kind = kind( 1.0d0 ) ) :: ai( nxl : nxh, nyl : nyh, nzl : nzh )
  real( kind = kind( 1.0d0 ) ) :: wrk( idwrk )
  !
  integer :: nx, ny, nz, ix, iy, iz, jx, jy, jz
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: tr( : , : , : ), ti( :, :, : )
  !
  nx = 1 + nxh - nxl
  ny = 1 + nyh - nyl
  nz = 1 + nzh - nzl
  allocate( tr( nx, ny, nz ), ti( nx, ny, nz ) )
  !
  do ix = nxl, nxh
     do iy = nyl, nyh
        do iz = nzl, nzh
           jx = 1 + ix
           if ( jx .le. 0 ) jx = jx + nx
           jy = 1 + iy
           if ( jy .le. 0 ) jy = jy + ny
           jz = 1 + iz
           if ( jz .le. 0 ) jz = jz + nz
           tr( jx, jy, jz ) = ar( ix, iy, iz )
           ti( jx, jy, jz ) = ai( ix, iy, iz )
        end do
     end do
  end do
  !
  call cfft( tr, ti, nx, nx, ny, nz, mode, wrk, idwrk )      

  do jx = 1, nx
     do jy = 1, ny
        do jz = 1, nz
           ix = jx - 1
           if ( ix .gt. nxh ) ix = ix - nx
           iy = jy - 1
           if ( iy .gt. nyh ) iy = iy - ny
           iz = jz - 1
           if ( iz .gt. nzh ) iz = iz - nz
           ar( ix, iy, iz ) = tr( jx, jy, jz )
           ai( ix, iy, iz ) = ti( jx, jy, jz )
        end do
     end do
  end do
  !
  deallocate( tr, ti )
  !
  return
end subroutine ctrdfft
