program mainact
! This is re written to closely match the format of the lanczos_non_herm routine.
! It does away with the 'cmd' , act, end etc, and the GMRES
  implicit none
  !
  integer, parameter :: adat = 45, fdat = 47,orth=48
  integer, parameter :: psitmp = 55
  integer, parameter :: eav = 80
  integer, parameter :: rcd = 46
  !
  integer, allocatable :: kk( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: ebd(:,:), eblda(:,:)
  real( kind = kind( 1.0d0 ) ), allocatable :: xk(:,:)
!  real( kind = kind( 1.0d0 ) ), allocatable :: a(:), b(:)
  complex( kind = kind( 1.0d0 ) ), allocatable :: a(:), b(:), c(:)
  real( kind = kind( 1.0d0 ) ), allocatable :: mr(:,:,:)
  complex( kind = kind( 1.0d0 ) ), allocatable :: u(:,:,:)

 
  ! variables with the shape of a wave function ... used to be ( nbv, nbc, nk ), now just ( n )
  real( kind = kind( 1.0d0 ) ), allocatable :: allow( : ), sfact( : ),tmparray(:)
  complex( kind = kind( 1.0d0 ) ), allocatable :: psi_p( : ), psi( : ), psi_m( : ), tmp( : ),ztmparray(:)
  complex( kind = kind( 1.0d0 ) ), allocatable :: backpsi_p( : ), backpsi( : ), backpsi_m( : ), backtmp( : ),&
                                                  tpsi(:),tbackpsi(:),pp(:),qq(:)
  !
  ! Arrays, parameters for doing partial reorthogonalization
  complex( kind = kind( 1.0d0 ) ), allocatable :: Planc( :, : ), Qlanc( :, : ),Pj(:,:),Qj(:,:)
  real( kind = kind( 1.0d0 ) ), allocatable :: fnorm(:),orth_lev(:)
  real( kind = kind( 1.0d0 ) ), parameter :: eps = 1*1.11E-16, eps1 = sqrt(eps)
  real( kind = kind( 1.0d0 ) ) :: frob,term,backterm,max_orth
!  logical, parameter :: full = .false. , biortho = .true. 
  logical :: full, biortho,sec_orth 
  integer :: orthtot, oterm, orthfreq,iterm1,iterm2,F_j,prbo2

  integer :: i, j, n
  real( kind = kind( 1.0d0 ) ) :: pnorm, vol, backpnorm,inv_qlength, qinb(3),bvec(3,3)
  !
  ! Lapack, BLAS functions
  integer :: vecdim
  real( kind = kind( 1.0d0 ) ) :: DZNRM2, DLANGE, DDOT,rterm1,rterm2
  complex( kind = kind( 1.0d0 ) ) :: ZDOTC,zterm,zterm1,zterm2,lnorm1,lnorm2


  integer :: ilftl, ilfth, irgtl, irgth
  integer :: nkx, nky, nkz, nk, iter
  integer :: ngx, ngy, ngz, ng, nbv, nbc, nb
  integer :: niter, bflag, lflag, backf, aldaf, qpflg, bwflg, bande, nonherm
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
  real( kind = kind( 1.0d0 ) ) :: f( 2 ), gres, gprc, ener, pi, fact, de, emax
  complex( kind = kind( 1.0d0 ) ) :: rm1
  complex( kind = kind( 1.0d0 ) ), allocatable :: v1( : ), v2( : ), cwrk( : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: x( : ), xprev( : ), rhs( : ), pcdiv( : ), bra( : ), xnew( : )
  character * 3 :: req, bs, as
  character * 5 :: eval
  character * 9 :: ct
  character * 7, parameter :: u7 = 'unknown'
  character * 9, parameter :: f9 = 'formatted'
  !
  pi = 4.0d0 * atan( 1.0d0 )
  rm1 = -1
  rm1 = sqrt( rm1 )
  !
  open( unit=99, file='kandb.h', form='formatted', status='unknown' )
  rewind 99
  call igetval( nkx, 'nkx  ' ); call igetval( nky, 'nky  ' ); call igetval( nkz, 'nkz  ' ); nk = nkx * nky * nkz
  call igetval( ngx, 'ngx  ' ); call igetval( ngy, 'ngy  ' ); call igetval( ngz, 'ngz  ' ); ng = ngx * ngy * ngz
  call igetval( nbv, 'nbv  ' ); call igetval( nbc, 'nbc  ' )
  call igetval( ilftl, 'ilftl' ); call igetval( ilfth, 'ilfth' ); call igetval( irgtl, 'irgtl' ); call igetval( irgth, 'irgth' )
! nb = irgth - ilftl + 1
! if ( nb .ne. nbv + nbc ) stop 'metal option not yet reinstated'
  if ( nbv .ne. 1 + ilfth - ilftl ) stop 'bad nbv'
  if ( nbc .ne. 1 + irgth - irgtl ) stop 'bad nbc'
  close( unit=99 )

  open( unit=99, file='niter.h', form='formatted', status='unknown' )
  call igetval( niter, 'niter' )
  call igetval( bflag, 'bflag' )
  call igetval( lflag, 'lflag' )
  call igetval( backf, 'backf' )
  call igetval( aldaf, 'aldaf' )
  call igetval( qpflg, 'qpflg' )
  call igetval( bwflg, 'bwflg' )
  call igetval( bande, 'bande' )
  call igetval( nonherm, 'nonhe' )
  close( unit=99 )


  ! reorthogonalization flags and parameters
  open( unit=99, file='ortho_flag', form='formatted', status='unknown' )
  rewind 99
  read(99,*) 
  read(99,*) iterm1,iterm2
  if(iterm1 .eq. 1)then
    biortho = .true.
  else
    biortho = .false.
  endif

  if(iterm2 .eq. 1)then
    full = .true.
  else
    full = .false.
  endif
  print *
  print *, " biortho =",biortho," full =",full
  print *
  close( unit=99 )

  ! Use this to do the second rebiorthog if a rebiorthog is done in the previous iteration
  prbo2 = 2


  ! Setting the number of bands
  if (bwflg .eq. 0) then
     nb = nbv + nbc
     write(6,*) " "
     write(6,*) "Doing TDA BSE, nb = ",nb
  else
     nb = nbv + nbc
     write(6,*) " "
     write(6,*) "Doing full BSE, nb = ",nb
     print *, "bande =",bande
     print *
     print *, "bflag =",bflag
     print *
     print *, "lflag =",lflag
     print *
     print *, "bwflag =",bwflg
     print *
     print *, "nonherm =",nonherm
     print *
  endif

  n = nbv * nbc * nk
  write(6,*) " "
  write(6,*) " nk = ",nk," nbv = ", nbv, " nbc = ",nbc, " n = ",n
  write(6,*) " "
  !

!  nonherm = 1 ! I will add this later to niter.h
 

  if ( (lflag .ne. 0) .and. ( (backf .ne. 0)  .or. (aldaf .ne. 0) ) ) stop 'ladder incompatible with Kxc or backwards pairs'
  !
  allocate( ebd( nb, nk ), eblda( nb, nk ) )
  allocate( u( nb, nk, ng ) )
  allocate( kk( nk, 3 ), xk( nk, 3 ) )
  allocate( a( 0 : niter ), b( 0: niter + 1 ), c( 0: niter + 1 ) )
  !
  ! allocate( allow( nbv, nbc, nk ) )
  ! allocate( psi_p( nbv, nbc, nk ), psi( nbv, nbc, nk ) )
  ! allocate( psi_m( nbv, nbc, nk ) )
  ! allocate( tmp( nbv, nbc, nk ) )
  allocate( allow( n ), sfact( n ), psi( n ), psi_p( n ), psi_m( n ))
  allocate( backpsi( n ), backpsi_p( n ), backpsi_m( n ))
  allocate(tpsi(n),tbackpsi(n),pp(n),qq(n))
  !
  ! Allocate arrays for saving Lanczos vectors
  ! Planc for the psi's, Qlanc for the backpsi's

  if (nonherm .ne. 0) then
     allocate( Planc( n, niter ), Qlanc( n,niter ), fnorm(niter),orth_lev(niter))
     Planc = 0.0
     Qlanc = 0.0
     fnorm = 0.0
     orth_lev = 0.0
  end if

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

  ! Reading in the q vector (in units of b) and the bvecs
  open( unit=99, file='qinunitsofbvectors.ipt', form=f9,status=u7)
  rewind 99
  read(99,*) qinb(:)
  write(6,*) " q in b",qinb
  write(6,*) " "
  close(unit=99)

  open( unit=99, file='bvecs', form=f9,status=u7)
  rewind 99
  read(99,*) bvec(:,:)
  write(6,*) " bvecs"
  write(6,*)bvec
  write(6,*) " "
  close(unit=99)

  ! Calculating the inverse q length to normalize ptmels.dat
  inv_qlength = (qinb(1)*bvec(1,1)+qinb(2)*bvec(1,2)+qinb(3)*bvec(1,3))**2 &
              + (qinb(1)*bvec(2,1)+qinb(2)*bvec(2,2)+qinb(3)*bvec(2,3))**2 &
              + (qinb(1)*bvec(3,1)+qinb(2)*bvec(3,2)+qinb(3)*bvec(3,3))**2

  inv_qlength = 1.0/sqrt(inv_qlength)
  write(6,*) " "
  write(6,*) " inverse q length = ",inv_qlength
  write(6,*) " "
  write(6,*) "inverse q length from OCEAN = 10029.4937358645"
  write(6,*) " "

  write(6,*) " Size of lanczos matrix of orthog vectors = ",n,"by",niter
  write(6,*) " "
  inv_qlength = 10029.493735864515

  ! get k-points
  call kmapr( nk, nkx, nky, nkz, xk, kk )
  !
  ! get initial state
!  call getwgt( nk, nbv, nbc, psi ) ! original
  call getwgt_new( nk, nbv, nbc,ilftl,ilfth,irgtl,irgth,n, psi ) ! for OCEAN data
  psi = inv_qlength*psi


  if (nonherm .ne. 0)backpsi=psi
  !
  ! get u( k ) functions at real-space grid points, x
!  call getumat( nb, nk, ng, u ) ! original
  call getumat( nb,ilftl,ilfth,irgtl,irgth, nk, ng, u )
  !
  call bdallck( nbv, nbc, nk, allow, sfact, ebd, eblda, qpflg, bwflg )
  write(6,*) " nbv = ",nbv,"nbc =",nbc,"nk = ",nk
!  write(6,*) " Now printing the allow matrix of shape: ",shape(allow)
!  write(6,*) " "
!  write(6,*) allow
!  write(6,*) " "
!  write(6,*) " sfact matrix:"
!  write(6,*) " "
!  write(6,*) sfact
  write(6,*) " "


  if ( backf .eq. 1 ) call bfnorm( psi, ebd, nk, nb, nbc, nbv )
  !
  call uphase( u, nb, nk, ng, ngx, ngy, ngz, xk )
  !
     
  loud = .true.

! Normalize psi and backpsi, psi == p, backpsi == q in my lanczos_non_herm routine
!   allow = 1.0

   pnorm = dble( wvfdot( psi, psi, allow, n ) )
   pnorm = sqrt( pnorm )

     ! back psi
   backpnorm = dble( wvfdot( backpsi, backpsi, allow, n ) )
   backpnorm = sqrt( backpnorm )

!     if ( loud ) write ( 6, '(1a6,1x,1e15.8)' ) 'pnorm=', pnorm
   if ( loud ) write ( 6, * ) 'pnorm=', pnorm, 'backpnorm=', backpnorm
   term = 1.0d0 / pnorm
   backterm = 1.0d0 / backpnorm

   ! ------- Starting the Lanczos iteration ---------
   psi( : ) = psi( : ) * term

     ! back psi
   backpsi( : ) = backpsi( : ) * backterm

   ! Add efect of allow
   psi = psi * allow
   backpsi = psi * allow

   ! Add the first normalized vectors
   Planc(:,1) = psi
   print *
   print *, " Starting the Lanczos algorithm "
   print *
   Qlanc(:,1) = backpsi


  ! Save info about reorthog, first # is orthfreq, second orthtot
   open( orth, file='orthinfo', form='formatted', status='unknown' )
   rewind(orth) 
   write(orth,*) "# iter,   oterm,   2*oterm,   orthfreq,   orthtot "
   orthtot = 0
   orthfreq = 0
   F_j = 0
   do iter=1,niter

      if (iter==1) then
         b(iter) = 0.0
         c(iter) = 0.0
      endif

      if (iter .ne. 1) then
        psi_m = Planc(:,iter-1)
        backpsi_m = Qlanc(:,iter-1)
      else
        psi_m = 0.0
        backpsi_m = 0.0
      endif

      write(6,*) " "
      write(6,*) "************ Lanczos iteration =",iter,"**********"

      ! Checking biorthogonality with frobenius norm
      allocate(Pj(n,iter),Qj(n,iter))
      Pj = 0.0
      Qj = 0.0
      Pj = Planc(:,1:iter)
      Qj = Qlanc(:,1:iter)
      call zcheck_orthog_frob(Pj,Qj,n,iter,frob)
      write(6,*) " "
      write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^"
      write(6,*) "frob =",frob
      write(6,*) " "
      write(6,*) "machine prec. =",eps
      write(6,*) " "
      write(6,*) "tolerance =",eps1
      write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^"
      write(6,*) " "
      fnorm(iter) = frob

      if (iter .eq. 1) then
        orth_lev(iter) = frob
      else
        orth_lev(iter) = max_orth
      endif
      deallocate(Pj,Qj)

   ! Act on psi and backpsi, getting tpsi, tbackpsi, according to lanczos_non_herm routine in lanczos_8.f90
      ! psi
      tpsi = psi
      write(6,*) " Acting on psi with H^dag, p in my formulation"
      write(6,*) " "
!      call act( tpsi, u, ebd, kk, nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv, &
!             xk, nkret, kret, mr, lflag, bflag, backf, aldaf, allow, sfact, bande )
      call act_cmplx_dag( tpsi, u, ebd, kk, nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv, &
             xk, nkret, kret, mr, lflag, bflag, backf, aldaf, allow, sfact, bande )

      !backpsi
      tbackpsi = backpsi
      write(6,*) " Acting on backpsi with H, q in my formulation"
      write(6,*) " "
!      call act( tbackpsi, u, ebd, kk, nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv, &
!             xk, nkret, kret, mr, lflag, bflag, backf, aldaf, allow, sfact, bande )
      call act_cmplx( tbackpsi, u, ebd, kk, nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv, &
             xk, nkret, kret, mr, lflag, bflag, backf, aldaf, allow, sfact, bande )
      zterm = wvfdot( psi, tbackpsi, allow, n ) 
      write(6,*) " "
      write(6,*) "===== Before Lanczos, overlap ======="
      write(6,*) " <psi|H|backpsi>", zterm

    ! First part of three term recurrence
    ! pp = tp - conjg(beta)*p_m-1, tp = H^dag*p
    ! qq = tq - gamma*q_(m-1) , tq = H*q, beta = b, gamma = c, alpha = a
      psi_p = tpsi - conjg(b(iter))*psi_m
      backpsi_p = tbackpsi - c(iter)*backpsi_m

    ! Calculating a(iter)
      zterm = wvfdot( psi, backpsi_p, allow, n ) 
      a(iter) = zterm
      write(6,*) " "
      write(6,*) "a(iter) =",zterm
      write(6,*) " "

    ! Other term of recurrence
      psi_p = psi_p - conjg(a(iter))*psi
      backpsi_p = backpsi_p - a(iter)*backpsi

      ! Doing a rebiorthog of psi_p and backpsi_p, can do no re orthog,full reorthog,and 
      ! partial reorthog
       
!    ! Second re-biorthog, dependent on whether in the iter-1 there occured an orthog,
!    ! It is also dependent on prbo2.
      if (prbo2 .eq. 2) then
              if (iter .eq. F_j) then
                 print *
                 print *, " Doing a second full re-biorthogonalization "
                 print *
                 print *," iter = ",iter," F_j =", F_j
                 print *
                 call zgs_blas_sing_two_side_x(Planc(:,1:iter),Qlanc(:,1:iter),n,iter,eps1,&
                                         psi_p,backpsi_p,pp,qq,max_orth) 
                 psi_p = pp
                 backpsi_p = qq
              endif
      endif
      if (biortho) then
         if (full) then
            write(6,*) "########  Doing a full re-biorthogonalization ########"
            ! The three different routines do the matrix*matrix*vector product differently, with
            ! the most robust but memory intensive being twoside_GS_sing_x, followed by zgs_blas_sing_two_side_x, and the least robust being
            ! twoside_GS_naive_sing_x

!            call twoside_GS_sing_x(Planc(:,1:iter),Qlanc(:,1:iter),n,iter,eps1,&
!                                 psi_p,backpsi_p,pp,qq) 
!            call twoside_GS_naive_sing_x(Planc(:,1:iter),Qlanc(:,1:iter),n,iter,eps1,&
!                                 psi_p,backpsi_p,pp,qq) 
            call zgs_blas_sing_two_side_x(Planc(:,1:iter),Qlanc(:,1:iter),n,iter,eps1,&
                                 psi_p,backpsi_p,pp,qq,max_orth) 
         else
            write(6,*) "########  Doing a partial re-biorthogonalization ########"
            write(6,*) " tolerance = ",eps1
            write(6,*) " "
            ! Three different ways of doing the matrix-matrix-vector product needed for rebiorthog
            ! the most robust but memory intensive for different kinds of nonhermitian hamiltonians is twoside_GS_sing_minprbo,
            ! then zgs_blas_sing_two_side_minprbo_x and then zgs_blas_sing_two_side_prbo_x does the PRBO-1 and PRBO-2

!            call twoside_GS_sing_minprbo(Planc,Qlanc,n,iter,eps1,&
!                                 psi_p,backpsi_p,pp,qq,max_orth)  ! better for challenging non hermitian calc, but needs more memory 
!            call zgs_blas_sing_two_side_minprbo_x(Planc,Qlanc,n,iter,eps1,&
!                                 psi_p,backpsi_p,pp,qq,oterm,max_orth)
            call zgs_blas_sing_two_side_prbo_x(Planc,Qlanc,n,iter,eps1,&
                                 psi_p,backpsi_p,pp,qq,oterm,max_orth,F_j)
!            print *, "oterm = ",oterm
            if (oterm .gt. 0) then
              orthfreq = orthfreq + 1
            endif 

         endif 
         psi_p = pp
         backpsi_p = qq
      else
        write(6,*) " "
        write(6,*) "####### No reorthogonalization is carried out #########"
        oterm = 0
        write(6,*) " "
      endif


    ! Calculating b(iter+1), c)iter+1)
      zterm = wvfdot( backpsi_p, psi_p, allow, n ) 
      lnorm1 = sqrt(abs(zterm))
      lnorm2 = conjg(zterm)/lnorm1
      write(6,*) "b(iter+1) = ",lnorm1
      write(6,*) "c(iter+1) = ",lnorm2
      write(6,*) " "

!      if (iter .ne. niter) then
      b(iter+1) = lnorm1
      c(iter+1) = lnorm2
!      endif

    ! Updating the new psi's and backpsi's for next iteration
      if (iter .ne. niter) then
         psi = (1.0/conjg(lnorm2))*psi_p
         backpsi = (1.0/lnorm1)*backpsi_p
     
    ! Add efect of allow
         psi = psi * allow
         backpsi = psi * allow

!         pnorm = dble( wvfdot( psi, psi, allow, n ) )
!         pnorm = sqrt( pnorm )
!
!     ! back psi
!         backpnorm = dble( wvfdot( backpsi, backpsi, allow, n ) )
!         backpnorm = sqrt( backpnorm )
!         print *
!         print *," After reorthog, pnorm=",pnorm,"backpnorm=",backpnorm
!         print *
     ! Add new Lanczos vectors
         Planc(:,iter+1) = psi
         Qlanc(:,iter+1) = backpsi  
!         zterm = wvfdot( psi, Qlanc(:,iter+1), allow, n ) 
!         zterm1 = wvfdot( Planc(:,iter+1), backpsi, allow, n ) 
!         print *,"Overlap of new vectors "
!         print *,"<p|q> =",zterm
!         print *
!         print *,"<q|p> =",zterm1
!         print *
      endif
!      print *
!      print *," # re orthogs =",2*oterm
      orthtot = orthtot + 2*oterm

      ! Update orthogs info
      write(orth,*)iter,oterm,2*oterm,orthfreq,orthtot
   enddo ! end of niter loop

!! Save info about reorthog, first # is orthfreq, second orthtot
!   open( unit=99, file='orthinfo', form='formatted', status='unknown' )
!   rewind 99
!   write(99,*)orthfreq,orthtot
   close( orth )

! Saving the Lanczos coefficients to file   
!           call adump( adat, 'a.dat', pnorm, nk, iter, a )
!           call adump( adat, 'b.dat', pnorm, nk, iter + 1, b )
   call cadump( adat, 'a.dat', pnorm, nk, iter, a )
   call cadump( adat, 'b.dat', pnorm, nk, iter + 1, b )
   call cadump( adat, 'c.dat', pnorm, nk, iter + 1, c )
   open( unit=99, file='ndone.h', form='formatted', status='unknown' )
   rewind 99
   call iputval( iter, 'ndone' )
   close( unit=99 )

! write the file fnorm which saves the level of orthogonality 
   if (nonherm .ne. 0) then
      open( unit=98, file='orth_lev.dat', form='formatted', status='unknown' )
      open( unit=99, file='f.dat', form='formatted', status='unknown' )
      open( unit=100, file='im_a.dat', form='formatted', status='unknown' )
      open( unit=101, file='orthog_sum', form='formatted', status='unknown' )
      rewind 98
      rewind 99
      rewind 100
      rewind 101
      do i=1,niter
         write(98,*)i,orth_lev(i)
         write(99,*)i,fnorm(i)
         write(100,*)i,abs(aimag(a(i)))
         write(101,*)i,int(i*(i+1)/2),i*(i+1)
      enddo
      close( unit=98 )
      close( unit=99 )
      close( unit=100 )
      close( unit=101 )
   endif
   deallocate( ebd, u, kk, xk, a, b, c, allow, psi_p, psi, psi_m,tpsi,tbackpsi,pp,qq )
   if (nonherm .ne. 0) deallocate( Planc, Qlanc, fnorm,orth_lev )
  !
end program mainact

