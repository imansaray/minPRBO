subroutine coulmult( ngz, ngy, ngx, nk, dnr, dni )
  implicit none
  !
  integer :: ngx, ngy, ngz, nk
  real( kind = kind( 1.0d0 ) ) :: dnr( ngz, ngy, ngx ), dni( ngz, ngy, ngx )
  !
  integer :: star, mstar, ix( 3, 48 ), i, ixx, iyy, izz, idum, j, k
  real( kind = kind( 1.0d0 ) ), allocatable :: enr( :, :, : ), eni( :, :, : )
  real( kind = kind( 1.0d0 ) ) :: gsqd, mul, omega, pi
  real( kind = kind( 1.0d0 ) ) :: qt( 3 ), qq( 3 ), bmet( 3, 3 )
  logical :: fault
  character * 60 :: fstr
  !
  pi = 4.0d0 * atan( 1.0d0 )
  allocate( enr( ngz, ngy, ngx ), eni( ngz, ngy, ngx ) )
  enr = dnr
  eni = dni
  dnr = 0
  dni = 0
  !
  open( unit=99, file='omega.h', form='formatted', status='unknown' )
  rewind 99
  call rgetval( omega, 'volbo' )
  call rgetval( qt( 1 ), 'qtra1' )
  call rgetval( qt( 2 ), 'qtra2' )
  call rgetval( qt( 3 ), 'qtra3' )
  close( unit=99 )
  open( unit=99, file='vecy', form='formatted', status='unknown' )
  rewind 99
  do i = 1, 3
     do j = 1, 3
        read ( 99, * ) bmet( i, j )
     end do
  end do 
  !
  open( unit=99, file='gvectors', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) idum ! skip G=0
  read ( 99, * ) idum
  fault = .false.
  do while ( .not. fault )
     read ( 99, * ) star, mstar
     i = 1
     do while ( ( i .le. mstar ) .and. ( .not. fault ) )
        read ( 99, * ) idum, idum, ix( :, i )
        if ( 2 * abs( ix( 1, i ) ) .ge. ngx ) fault = .true.
        if ( 2 * abs( ix( 2, i ) ) .ge. ngy ) fault = .true.
        if ( 2 * abs( ix( 3, i ) ) .ge. ngz ) fault = .true.
        i = i + 1
     end do
     if ( .not. fault ) then
        fstr = '(1i8,2(1x,1e15.8))'
        do i = 1, mstar
           do j = 1, 3
              qq( j ) = qt( j ) + ix( j, i )
           end do
           gsqd = 0
           do j = 1, 3
              do k = 1, 3
                 gsqd = gsqd + qq( j ) * qq( k ) * bmet( j, k )
              end do
           end do 
           mul = 4.0d0 * pi * 27.2114d0 / ( nk * omega * gsqd )
           ixx = 1 + ix( 1, i )
           iyy = 1 + ix( 2, i )
           izz = 1 + ix( 3, i )
           if ( ixx .le. 0 ) ixx = ixx + ngx
           if ( iyy .le. 0 ) iyy = iyy + ngy
           if ( izz .le. 0 ) izz = izz + ngz
           dnr( izz, iyy, ixx ) = enr( izz, iyy, ixx ) * mul
           dni( izz, iyy, ixx ) = eni( izz, iyy, ixx ) * mul
        end do
     end if
  end do
  close( unit=99 )
  !
  deallocate( enr, eni ) 
  !
  return
end subroutine coulmult
