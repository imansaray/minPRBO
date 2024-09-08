program mainact
  implicit none
  !
  integer, parameter :: adat = 45
  integer, parameter :: psitmp = 55
  integer, parameter :: eav = 80
  integer, parameter :: rcd = 46
  !
  integer, allocatable :: kk( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: ebd(:,:), eblda(:,:)
  real( kind = kind( 1.0d0 ) ), allocatable :: xk(:,:)
  real( kind = kind( 1.0d0 ) ), allocatable :: a(:), b(:)
  real( kind = kind( 1.0d0 ) ), allocatable :: mr(:,:,:)
  complex( kind = kind( 1.0d0 ) ), allocatable :: u(:,:,:)
  !
  ! variables with the shape of a wave function ... used to be ( nbv, nbc, nk ), now just ( n )
  real( kind = kind( 1.0d0 ) ), allocatable :: allow( : ), sfact( : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: psi_p( : ), psi( : ), psi_m( : ), tmp( : )
  !
  integer :: i, j, n
  real( kind = kind( 1.0d0 ) ) :: pnorm, vol,inv_qlength,qinb(3),bvec(3,3)
  !
  integer :: ilftl, ilfth, irgtl, irgth
  integer :: nkx, nky, nkz, nk, iter
  integer :: ngx, ngy, ngz, ng, nbv, nbc, nb
  integer :: niter, bflag, lflag, backf, aldaf, qpflg, bwflg, bande
  complex( kind = kind( 1.0d0 ) ), external :: hilbcd
  !
  complex( kind = kind( 1.0d0 ) ), external :: wvfdot
  !
  character * 3 :: cmd, meth
  logical :: loud
  !
  integer :: nkret
  integer, allocatable :: kret( : )
  !
  integer :: nloop, iwrk, i1, i2, need
  real( kind = kind( 1.0d0 ) ) :: f( 2 ), gres, gprc, ener, pi, fact, de, emax,epsre,epsim,rterm
  complex( kind = kind( 1.0d0 ) ) :: rm1,zterm
  complex( kind = kind( 1.0d0 ) ), allocatable :: v1( : ), v2( : ), cwrk( : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: x( : ), xprev( : ), rhs( : ), pcdiv( : ), bra( : ), xnew( : )
  character * 3 :: req, bs, as
  character * 5 :: eval
  character * 9 :: ct
  character * 7, parameter :: u7='unknown'
  character * 9, parameter :: f9='formatted'
  !
  pi = 4.0d0 * atan( 1.0d0 )
  rm1 = -1
  rm1 = sqrt( rm1 )
  !
  open( unit=99, file='kandb.h', form='formatted', status='unknown' )
  call igetval( nkx, 'nkx  ' ); call igetval( nky, 'nky  ' ); call igetval( nkz, 'nkz  ' ); nk = nkx * nky * nkz
  call igetval( ngx, 'ngx  ' ); call igetval( ngy, 'ngy  ' ); call igetval( ngz, 'ngz  ' ); ng = ngx * ngy * ngz
  call igetval( nbv, 'nbv  ' ); call igetval( nbc, 'nbc  ' )
  call igetval( ilftl, 'ilftl' ); call igetval( ilfth, 'ilfth' ); call igetval( irgtl, 'irgtl' ); call igetval( irgth, 'irgth' )
! nb = irgth - ilftl + 1
! if ( nb .ne. nbv + nbc ) stop 'metal option not yet reinstated'
  if ( nbv .ne. 1 + ilfth - ilftl ) stop 'bad nbv'
  if ( nbc .ne. 1 + irgth - irgtl ) stop 'bad nbc'
  close( unit=99 )
  nb = nbv + nbc
  n = nbv * nbc * nk
  !
  open( unit=99, file='niter.h', form='formatted', status='unknown' )
  call igetval( niter, 'niter' )
  call igetval( bflag, 'bflag' )
  call igetval( lflag, 'lflag' )
  call igetval( backf, 'backf' )
  call igetval( aldaf, 'aldaf' )
  call igetval( qpflg, 'qpflg' )
  call igetval( bwflg, 'bwflg' )
  call igetval( bande, 'bande' )
  close( unit=99 )

!  if ( (lflag .ne. 0) .and. ( (backf .ne. 0)  .or. (aldaf .ne. 0) ) ) stop 'ladder incompatible with Kxc or backwards pairs'
  !
  allocate( ebd( nb, nk ), eblda( nb, nk ) )
  allocate( u( nb, nk, ng ) )
  allocate( kk( nk, 3 ), xk( nk, 3 ) )
  allocate( a( 0 : niter ), b( 0: niter + 1 ) )
  !
  ! allocate( allow( nbv, nbc, nk ) )
  ! allocate( psi_p( nbv, nbc, nk ), psi( nbv, nbc, nk ) )
  ! allocate( psi_m( nbv, nbc, nk ) )
  ! allocate( tmp( nbv, nbc, nk ) )
  allocate( allow( n ), sfact( n ), psi( n ), psi_p( n ), psi_m( n ), tmp( n ) )
  !
  open( unit=99, file='ladder.dat', form='unformatted', status='unknown' )
  rewind 99
  read ( 99 ) nkret
  allocate( kret( nkret ) )
  allocate( mr( nkret, ng, ng ) )
  read ( 99 ) kret
  do i = 1, ng
     open( unit=98, file='prog', form='formatted', status='unknown' )
     rewind 98
     write ( 98, '(1a11,3x,2i6)' ) 'ladder read', i, ng
     close( unit=98 )
     do j = 1, ng
        read ( 99 ) mr( :, j, i )
     end do
  end do
  close( unit=99 )
  !

  ! Reading in the q vectors (in units of b) and the bvecs
  open( unit=99, file='qinunitsofbvectors.ipt', form=f9,status=u7)
  rewind 99
  read(99,*) qinb(:)
  print *
  print *, "q in b",qinb
  print *
  close(unit=99)

  open(unit=99, file='bvecs', form=f9,status=u7)
  rewind 99
  read(99,*) bvec(:,:)
  print *," bvecs"
  print *,bvec
  print *
  close(unit=99)
  
  ! Calculating the inverse q length to normalize ptmels.dat
  inv_qlength = (qinb(1)*bvec(1,1)+qinb(2)*bvec(1,2)+qinb(3)*bvec(1,3))**2 &
              + (qinb(1)*bvec(2,1)+qinb(2)*bvec(2,2)+qinb(3)*bvec(2,3))**2 &
              + (qinb(1)*bvec(3,1)+qinb(2)*bvec(3,2)+qinb(3)*bvec(3,3))**2

  inv_qlength = 1.0d0/sqrt(inv_qlength)
  print *
  print *," inverse q length = ",inv_qlength
  print *
  print *," inverse q length from OCEAN = 10029.4937358645"
  print *  

  ! get k-points
  call kmapr( nk, nkx, nky, nkz, xk, kk )
  !
  ! get initial state
!  call getwgt(nk,nbv,nbc,psi) ! original 
  call getwgt_new(nk,nbv,nbc,ilftl,ilfth,irgtl,irgth,n,psi) !for OCEAN data
  psi = inv_qlength*psi 
  !
  ! get u( k ) functions at real-space grid points, x
!  call getumat( nb, nk, ng, u ) ! original
  call getumat( nb,ilftl,ilfth,irgtl,irgth, nk, ng, u ) ! original

  !
  call bdallck( nbv, nbc, nk, allow, sfact, ebd, eblda, qpflg, bwflg )
  if ( backf .eq. 1 ) call bfnorm( psi, ebd, nk, nb, nbc, nbv )
  !
  call uphase( u, nb, nk, ng, ngx, ngy, ngz, xk )
  !
!  print *, ' Now enter a method : hay or inv '
!  read ( 5, '(1a3)' ) meth
  meth='hay'
  print *
  print *," Calculation type : haydock/Hermitian Lanczos('hay') or GMRES('inv')"
  print *," Doing a  '",meth, "'  type of calculation"
  print *
  select case( meth )
  case ( 'hay' )
     cmd = '---'
     loud = .true.
     do while ( cmd .ne. 'end' )
        call haydock( n, psi, psi_p, psi_m, allow, pnorm, niter, iter, a, b, loud, cmd, tmp )
        if ( cmd .eq. 'act' ) then
           write(6,*) 'call act'
           call act( psi_p, u, ebd, kk, nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv, xk, nkret, kret, mr, &
                     lflag, bflag, backf, aldaf, allow, sfact, bande )
        end if
        if ( cmd .eq. 'opt' ) then
           call adump( adat, 'a.dat', pnorm, nk, iter, a )
           call adump( adat, 'b.dat', pnorm, nk, iter + 1, b )
           open( unit=99, file='ndone.h', form='formatted', status='unknown' )
           rewind 99
           call iputval( iter, 'ndone' )
           close( unit=99 )
        end if
     end do
  case ( 'inv' )
     open( unit=100, file='opconsi_gmres.dat', form='formatted', status='unknown' )
     rewind 100
     write(100,*) "# energy,  e1,   e2,   ELF "

     write ( 98, '(1a11,3x,2i6)' ) 'ladder read', i, ng
     open( unit=99, file='omega.h', form='formatted', status='unknown' )
     call rgetval( vol, 'volbo' )
     close( unit=99 )
     print *, " Enter the following parameter values:"
     print *,"nloop,  gres, gprc, tol, energy, de, emax"
     print *
     open(unit=99, file='gmres_in', form='formatted', status='unknown')
     rewind 99
     read(99,*)
     read ( 99, * ) nloop, gres, gprc, f( 1 ), ener, de, emax
     close(unit=99)
     print *, nloop, gres, gprc, f( 1 ), ener, de, emax
     print *
     eval = 'zerox'
     allocate( x( n ), xprev( n ), xnew( n ) )
     do while ( ener .lt. emax )
     print *
     print *,"------- Energy = ",ener,"----------"
     print *
     xnew = x + x - xprev
     xprev = x
     x = xnew
     ! nloop = depth of GMRES, gres = resolution gamma, gprc = precondition gamma, f( 1 ) = tolerance, ener = E in eV
     fact = 8.0d0 * pi * 27.2114d0 / ( vol * dble( nk ) )
     write ( 6, '(1p,1a7,1e15.8)' ) 'fact = ', fact
     fact = sqrt( fact )
     write ( 6, '(1p,1a28,1e15.8)' ) 'state-vector scaling fact = ', fact
     iwrk = 1
     allocate( rhs( n ), bra( n ), v1( n ), v2( n ), cwrk( iwrk ), pcdiv( n ) )
     rhs( : ) = fact * psi( : )
     rhs( : ) = rhs( : ) * sfact( : )
     bra( : ) = rhs( : )
     v1( : ) = 1
     v1( : ) = v1( : ) * sfact( : )
      write(6,*) 'call act'
     call act( v1, u, ebd, kk, nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv, xk, nkret, kret, mr, &
               0, 0, backf, 0, allow, sfact, bande )
     do i = 1, n
        if ( abs( v1( i ) ) .lt. 0.00000001d0 ) v1( i ) = 1.0d0
        pcdiv( i ) = sfact( i ) * ( ener - v1( i ) - rm1 * gprc ) / ( ( ener - v1( i ) ) ** 2 + gprc ** 2 )
     end do  
     ct = 'beginning'
     req = '---'
     do while ( req .ne. 'end' )
        call invdrv( x, rhs, n, i1, i2, nloop, need, iwrk, cwrk, v1, v2, bs, as, req, ct, eval, f )
        select case( req )
        case( 'all' ) ! meaning, allocate an array as shown
           deallocate( cwrk )
           iwrk = need
           allocate( cwrk( need ) )
        case( 'act' ) ! meaning, multiply by SE-S(SH0+H1)S = S ( S ( E - H0 ) - H1 ) S ... in what follows, v1 must be untouched
           v2( : ) = v1( : )
           v2( : ) = v2( : ) * sfact( : )
           write(6,*) 'call act'
           call act( v2, u, ebd, kk, nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv, xk, nkret, kret, mr, &
                     lflag, bflag, backf, aldaf, allow, sfact, bande )
           v2( : ) = ( ener + rm1 * gres ) * v1( : ) - v2( : )
           v2( : ) = v2( : ) * sfact( : )
        case( 'prc' ) ! meaning, divide by S(E-H0) ... in what follows, v1 must be untouched
           v2( : ) = v1( : ) * pcdiv( : )
           write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) i1, i2, nloop, f( 2 ), f( 1 ), ener, 1.0d0 - dot_product( bra, x )
           write ( 66, '(1p,2x,3i5,5(1x,1e15.8))' ) i1, i2, nloop, f( 2 ), f( 1 ), ener, 1.0d0 - dot_product( bra, x )
        end select
     end do
     write ( 66, * )
     zterm = 1.0d0 - dot_product(bra,x)
     epsre = real(zterm)
     epsim = aimag(zterm)
     rterm = epsim/(epsre**2 + epsim**2)

!     write ( 76, '(1p,1i5,3(1x,1e15.8))' ) i1, ener, 1.0d0 - dot_product( bra, x )
!     write ( 76, '(1p,1i5,4(1x,1e15.8))' ) i1, ener, 1.0d0 - dot_product( bra, x ),rterm
     write ( 100, '(1p,4(1x,1e15.8))' ) ener, 1.0d0 - dot_product( bra, x ),rterm
     print *
     print *, '---- spectra -----'
     write ( 6, '(1p,1i5,4(1x,1e15.8))' ) i1, ener, 1.0d0 - dot_product( bra, x ),rterm
     print *
     deallocate( bra, rhs, v1, v2, cwrk, pcdiv )
     ener = ener + de
     if ( eval .eq. 'zerox' ) xprev = x
     eval = 'havex'
     end do 
     deallocate( x, xprev, xnew )
     close(100)
  end select
  !
  deallocate( ebd, u, kk, xk, a, b, tmp, allow, psi_p, psi, psi_m )
  !
end program mainact
