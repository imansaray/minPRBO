program tid
  implicit none
  integer n, i
  !---------------------------------------------------v
  character * 3 req, bs, as
  character * 5 ev
  character * 9 ct
  integer i1, i2, n2, need, iwrk
  double precision :: f( 2 )
  complex( kind = kind( 1.0d0 ) ) :: rm1
  double complex, allocatable :: v1( : ), v2( : ), cwrk( : )
  !---------------------------------------------------^
  double complex, allocatable :: x( : ), b( : )
  double complex, allocatable :: a( : ) 
  double complex, allocatable :: e( : ), be( : )
  double complex, allocatable :: xtrue( : )
  double precision el, eh, v, step, im, err
  double complex, external :: hilbcd
  read ( 5, * ) n, n2, el, eh, im, f( 1 )
  rm1 = -1
  rm1 = sqrt( rm1 )
  allocate( x( n ), b( n ), a( n ) )
  step = ( eh - el ) / dble( n - 1 )
  v = el
  do i = 1, n
     a( i ) = cmplx( v, im )
     v = v + step
  end do
  b = 1
  allocate( xtrue( n ), e( n ), be(n ) )
  xtrue = b / a
  !---------------------------------------------------v
  iwrk = 1
  allocate( cwrk( iwrk ), v1( n ), v2( n ) )
  ev = 'zerox'
  ct = 'beginning'
  req = '---'
  do while ( req .ne. 'end' )
     call invdrv( x, b, n, i1, i2, n2, need, iwrk, cwrk, v1, v2, bs, as, req, ct, ev, f )
     select case( req )
     case( 'all' )
        deallocate( cwrk )
        iwrk = need
        allocate( cwrk( need ) )
        !  in these two options, act and prc, must leave v1 untouched.
     case( 'act' )
        v2 = a * v1        ! HOW I'M ACTING
     case( 'prc' )
        v2( : ) = v1( : ) / ( a( : ) + 5.0d0 * rm1 )            ! HOW I'M PRECONDITIONING (NOT)
        e = xtrue - x
        be = a * e
        err = hilbcd( n, be, e )
        write ( 6, '(2x,3i5,3e16.8)' ) i1, i2, n2, f( 2 ), f( 1 ), err
     end select
  end do
  deallocate( cwrk, v1, v2 )
  !---------------------------------------------------^
  deallocate( x, b, a )
end program tid
