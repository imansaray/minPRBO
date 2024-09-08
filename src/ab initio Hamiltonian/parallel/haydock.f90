! Original Hermitian Lanczos

subroutine haydock( n, psi, psi_p, psi_m, allow, pnorm, niter, iter, a, b, loud, cmd, tmp )
  !
  integer :: n, niter, iter
  real( kind = kind( 1.0d0 ) ) :: pnorm, allow( n )
  real( kind = kind( 1.0d0 ) ) :: a( 0 : niter ), b( 0 : niter + 1 )
  complex( kind = kind( 1.0d0 ) ) :: psi( n ), psi_p( n ), psi_m( n ), tmp( n )
  logical :: loud
  character * 3 :: cmd
  !
  real( kind = kind( 1.0d0 ) ) :: pstot, x
  complex( kind = kind( 1.0d0 ) ), external :: wvfdot
  !
  select case( cmd )
  case ( '---' )
     pnorm = dble( wvfdot( psi, psi, allow, n ) )
     do i = 1, n
        write ( 95, '(1i8,2(1x,1e15.8))' ) i, psi( i )
        write ( 955, '(1i8,(1x,1e15.8))' ) i, allow( i )
     end do
     pnorm = sqrt( pnorm )
     if ( loud ) write ( 6, '(1a6,1x,1e15.8)' ) 'pnorm=', pnorm
     x = 1.0d0 / pnorm
     psi( : ) = psi( : ) * x
     psi_p( : ) = allow( : ) * psi( : )
     psi_m( : ) = 0
     b = 0
     iter = 0
     cmd = 'clk'
  case ( 'clk' )
     if ( loud ) write ( 6, * ) ' iter. no. ', iter
     cmd = 'act'
  case ( 'act' )
     a( iter ) = dble( wvfdot( psi, psi_p, allow, n ) )
     psi_p( : ) = psi_p( : ) - a( iter ) * psi( : ) - b( iter ) * psi_m( : )
     b( iter + 1 ) = dble( wvfdot( psi_p, psi_p, allow, n ) )
     b( iter + 1 ) = sqrt( b( iter + 1 ) )
     x = 1.0d0 / b( iter + 1 )   
     psi_p( : ) = psi_p( : ) * x
     pstot = wvfdot( psi_p, psi_p, allow, n )
     if ( loud ) write ( 6, '(1a6,1x,1f12.8)' ) 'pstot=', pstot
     psi_m( : ) = allow( : ) * psi( : )
     psi( : ) = allow( : ) * psi_p( : )
     if ( loud ) write ( 6, '(2x,1a2,1f10.5,3x,1a2,1f10.5)' ) 'a=', a( iter ), 'b=', b( iter + 1 )
     cmd = 'opt'
  case ( 'opt' )
     iter = iter + 1
     if ( iter .le. niter ) then
        cmd = 'clk'
     else
        cmd = 'end'
     end if
  end select
  !
  return
end subroutine haydock



subroutine nonherm_haydock( n, psi, psi_p, psi_m, allow, pnorm, &
                            backpsi, backpsi_p, backpsi_m, backpnorm, &
                            niter, iter, a, b, c, loud, cmd, tmp )
  !
  integer :: n, niter, iter
  real( kind = kind( 1.0d0 ) ) :: pnorm, backpnorm, allow( n )
!  real( kind = kind( 1.0d0 ) ) :: a( 0 : niter ), b( 0 : niter + 1 )
  complex( kind = kind( 1.0d0 ) ) :: a( 0 : niter ), b( 0 : niter + 1 ), c( 0 : niter + 1 )
  complex( kind = kind( 1.0d0 ) ) :: psi( n ), psi_p( n ), psi_m( n ), tmp( n ), &
                                     backpsi( n ), backpsi_p( n ), backpsi_m( n ), backtmp( n )
  complex( kind = kind( 1.0d0 ) ) :: zterm, zterm2, zx, backzx

  logical :: loud
  character * 3 :: cmd
  !
  real( kind = kind( 1.0d0 ) ) :: pstot, x, backpstot, backx, term
  complex( kind = kind( 1.0d0 ) ), external :: wvfdot

  ! Lapack routines
  complex( kind = kind( 1.0d0 ) ):: ZDOTC
  !
  select case( cmd )
  case ( '---' )
     pnorm = dble( wvfdot( psi, psi, allow, n ) )
     pnorm = sqrt( pnorm )

     ! back psi
     backpnorm = dble( wvfdot( backpsi, backpsi, allow, n ) )
     backpnorm = sqrt( backpnorm )

!     if ( loud ) write ( 6, '(1a6,1x,1e15.8)' ) 'pnorm=', pnorm
     if ( loud ) write ( 6, * ) 'pnorm=', pnorm, 'backpnorm=', backpnorm
     x = 1.0d0 / pnorm
     backx = 1.0d0 / backpnorm

     psi( : ) = psi( : ) * x
     psi_p( : ) = allow( : ) * psi( : )
     psi_m( : ) = 0
     b = 0

     ! back psi
     backpsi( : ) = backpsi( : ) * backx
     backpsi_p( : ) = allow( : ) * backpsi( : )
     backpsi_m( : ) = 0
     c = 0
     iter = 0
     cmd = 'clk'
  case ( 'clk' )
     if ( loud ) write ( 6, * ) ' iter. no. ', iter
     cmd = 'act'
  case ( 'act' )
     ! First part of three term recurrence     
     psi_p( : ) = psi_p( : ) - b( iter ) * psi_m( : )
     backpsi_p( : ) = backpsi_p( : ) - conjg(c( iter )) * backpsi_m( : )
 
     ! Calculating alpha/a
     a( iter ) = dble( wvfdot( backpsi, psi_p, allow, n ) )

     ! Second part of three term recurrence
     psi_p( : ) = psi_p( : ) - a( iter ) * psi( : )
     backpsi_p( : ) = backpsi_p( : ) - conjg(a( iter )) * backpsi( : )

     ! Calculating w_j = <psi_p,backpsi_p>
     zterm = dble( wvfdot( psi_p, backpsi_p, allow, n ) )
     if ( loud ) write ( 6, * ) 'w_j, <psi_p|backpsi_p=', zterm
     
     ! updating coefficients and new vectors
     if (iter .ne. niter) then
         zterm2 = sqrt(abs(zterm))
         c( iter + 1 ) = zterm2
         b( iter + 1 ) = conjg(zterm)/zterm2

         zx = 1.0d0 / b( iter + 1 )   
         backzx = 1.0d0 / c( iter + 1 )  

      ! new vectors  
         psi_p( : ) = psi_p( : ) * backzx
         backpsi_p( : ) = backpsi_p( : ) * conjg(zx)

      ! new norms   
         pstot = wvfdot( psi_p, psi_p, allow, n )
         backpstot = wvfdot( backpsi_p, backpsi_p, allow, n )

!         if ( loud ) write ( 6, '(1a6,1x,1f12.8)' ) 'pstot=', pstot
         if ( loud ) write ( 6, * ) 'pstot=', pstot, 'backpstot=', backpstot
    
      !  p_(j-1), q_(j-1)   
         psi_m( : ) = allow( : ) * psi( : )
         psi( : ) = allow( : ) * psi_p( : )

         backpsi_m( : ) = allow( : ) * backpsi( : )
         backpsi( : ) = allow( : ) * backpsi_p( : )

        !     if ( loud ) write ( 6, '(2x,1a2,1f10.5,3x,1a2,1f10.5)' ) 'a=', a( iter ), 'b=', b( iter + 1 )
         if ( loud ) write ( 6, * ) 'a=', a( iter ), 'b=', b( iter + 1 ),'c=',c(iter + 1)
         cmd = 'opt'
      else
         cmd = 'end'
      end if
  case ( 'opt' )
     iter = iter + 1
     if ( iter .le. niter ) then
        cmd = 'clk'
     else
        cmd = 'end'
     end if
  end select
  !
  return
end subroutine nonherm_haydock


subroutine nonherm_haydock_reorth( n, psi, psi_p, psi_m, allow, pnorm, &
                            backpsi, backpsi_p, backpsi_m, backpnorm, &
                            niter, iter, a, b, c, loud, cmd, tmp, &
                            Planc,Qlanc,fnorm, full, eps )

  ! Non hermitian lanczos with Partial Rebiorthogonalization,
  ! following the implementation in my lanc_nonherm_pro_1 routine in
  ! lanczos_8.f90

  integer :: n, niter, iter,i,j,k,l,iterm,jterm,f
  real( kind = kind( 1.0d0 ) ) :: pnorm, backpnorm, allow( n )
!  real( kind = kind( 1.0d0 ) ) :: a( 0 : niter ), b( 0 : niter + 1 )
  complex( kind = kind( 1.0d0 ) ) :: a( 0 : niter ), b( 0 : niter + 1 ), c( 0 : niter + 1 )
  complex( kind = kind( 1.0d0 ) ) :: psi( n ), psi_p( n ), psi_m( n ), tmp( n ), &
                                     backpsi( n ), backpsi_p( n ), backpsi_m( n ), backtmp( n )
  complex( kind = kind( 1.0d0 ) ) :: zterm, zterm2, zx, backzx

  logical :: loud, full
  character * 3 :: cmd
  !
  real( kind = kind( 1.0d0 ) ) :: pstot, x, backpstot, backx, term
  complex( kind = kind( 1.0d0 ) ), external :: wvfdot
  !
  ! Arrays and things for partial reorthog
  complex( kind = kind( 1.0d0 ) ) :: Planc(n,niter),Qlanc(n,niter)
  real( kind = kind( 1.0d0 ) ) :: fnorm(niter),eps,eps1,term_re,term_im,frob,fterm,tol
  integer :: orthtot,poterm,qoterm,chk,oterm
  complex( kind = kind( 1.0d0 ) ) :: term1,lnorm1,lnorm2,vp(n),vq(n),wp(n),wq(n),pp(n),qq(n)
  integer,allocatable :: porthopt(:),qorthopt(:),porthindx(:),qorthindx(:)
  complex( kind = kind( 1.0d0 ) ),allocatable :: iden(:,:),Pj(:,:),Qj(:,:),iden1(:,:),&
                                                 iden2(:,:),pw_jk(:),qw_jk(:),Clanc(:,:),Dlanc(:,:)
  ! Lapack routines
  complex( kind = kind( 1.0d0 ) ):: ZDOTC

! Setting the precision
  eps1 = eps
  select case( cmd )
  case ( '---' )
     pnorm = dble( wvfdot( psi, psi, allow, n ) )
     pnorm = sqrt( pnorm )

     ! back psi
     backpnorm = dble( wvfdot( backpsi, backpsi, allow, n ) )
     backpnorm = sqrt( backpnorm )

!     if ( loud ) write ( 6, '(1a6,1x,1e15.8)' ) 'pnorm=', pnorm
     if ( loud ) write ( 6, * ) 'pnorm=', pnorm, 'backpnorm=', backpnorm
     x = 1.0d0 / pnorm
     backx = 1.0d0 / backpnorm

     psi( : ) = psi( : ) * x
     psi_p( : ) = allow( : ) * psi( : )
!     psi_p( : ) = psi( : )
     psi_m( : ) = 0
     b = 0

     ! back psi
     backpsi( : ) = backpsi( : ) * backx
     backpsi_p( : ) = allow( : ) * backpsi( : )
!     backpsi_p( : ) = backpsi( : )
     backpsi_m( : ) = 0
     c = 0

     ! Add the initial vectors to Planc , Qlanc
     ! psi_p and backpsi_p are like my p and q
     Planc(:,1) = psi_p
     Qlanc(:,1) = backpsi_p
!     do i=1,n
!       write(6,*) "Planc(i,1) =",Planc(i,1)
!       write(6,*) " "
!     end do

     ! initialize total num of orthogonalizations
!     orthtot = 0

     iter = 0
     cmd = 'clk'
  case ( 'clk' )
     if ( loud ) write ( 6, * ) ' iter. no. ', iter
     cmd = 'act'
  case ( 'act' )

     ! Check level of orthogonality of vectors
     if (iter .ne. 0) then
        allocate(Pj(n,iter),Qj(n,iter))
        Pj = 0.0
        Qj = 0.0
        Pj = Planc(:,1:iter)
        Qj = Qlanc(:,1:iter)
        call zcheck_orthog_frob(Pj,Qj,n,iter,frob)
!        call zcheck_orthog_frob(Planc(:,1:iter),Qlanc(:,1:iter),n,iter,frob)
        write(6,*) " "
        write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
        write(6,*) "iteration = ",iter,"frob =",frob,"prec. =",eps1
        write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
        fnorm(iter) = frob
!        write(6,*) "fnorm(iter) =",fnorm(iter)
!        write(6,*) " "
        deallocate(Pj,Qj)
     end if    

     ! Now using Planc and Qlanc to set psi_m,backpsi_m.
     ! Need to double check this
!     if (iter .ne. 1) then
!        psi_m = Planc(:,iter-1)
!        backpsi_m = Qlanc(:,iter-1)
!     else
!        psi_m = 0.0
!        backpsi_m = 0.0
!     end if 

     ! First part of three term recurrence     
     psi_p( : ) = psi_p( : ) - b( iter ) * psi_m( : )
     backpsi_p( : ) = backpsi_p( : ) - conjg(c( iter )) * backpsi_m( : )
 
     ! Calculating alpha/a
     a( iter ) = wvfdot( backpsi, psi_p, allow, n ) 

     ! Second part of three term recurrence
     psi_p( : ) = psi_p( : ) - a( iter ) * psi( : )
     backpsi_p( : ) = backpsi_p( : ) - conjg(a( iter )) * backpsi( : )

     ! Calculating w_j = <psi_p,backpsi_p>
     zterm = wvfdot( psi_p, backpsi_p, allow, n ) 
     if ( loud ) write ( 6, * ) 'w_j, <psi_p|backpsi_p=', zterm
     write(6,*) " "

     ! updating coefficients and new vectors
     if (iter .ne. niter) then

!      ! Start checking for loss of orthogonality and doing partial reorthog
!         if(iter .ne. 0) then
!
!
!!! Do a Gram Schmidt orthog of new vector with previous basis  
!!!        subroutine twoside_GS_sing(A,B,n,iter,eps1,p,q,pp,qq)
!!            call twoside_GS_sing_pro(Planc(:,1:iter),Qlanc(:,1:iter),n,iter,eps1,&
!!                                 psi_p,backpsi_p,pp,qq)
!!            call twoside_GS_sing_x(Planc(:,1:iter),Qlanc(:,1:iter),n,iter,eps1,&
!!                                 psi_p,backpsi_p,pp,qq)
!            call twoside_GS_sing_x(Planc(:,1:iter),Qlanc(:,1:iter),n,iter,eps1,&
!                                 psi_p,backpsi_p,pp,qq)
!            psi_p = pp
!            backpsi_p = qq
!
!        end if !iter .ne. 0
!         
!         write(6,*) " "
!         write(6,*) "Done with partial re-biorthogonalization"
       
! ! Adding the new vectors to Planc, Qlanc
!
!         psi( : ) = allow( : ) * psi_p( : )
!
!!         backpsi_m( : ) = allow( : ) * backpsi( : )!, now set at beginning
!         backpsi( : ) = allow( : ) * backpsi_p( : )
!
!         Planc(:,iter+1) = psi
!         Qlanc(:,iter+1)=backpsi
!
!! Reorthogonalize the new iter+1 vectors        
!        allocate(Clanc(n,iter+1),Dlanc(n,iter+1))
!        Clanc = 0.0
!        Dlanc = 0.0
!        write(6,*) " "
!!              write(6,*) " Starting the first full biorthog"
!        write(6,*) " Starting Gram Schmidt with reorthog"
!        write(6,*) " "
!!        call zgs_class_full_two_side_x(Planc(:,1:iter+1),Qlanc(:,1:iter+1),n,iter+1,&
!!                                    Clanc,Dlanc)
!        call zgs_blas_full_two_side_x(Planc(:,1:iter+1),Qlanc(:,1:iter+1),n,iter+1,&
!                                    Clanc,Dlanc)
!!              call zgs_mod_full_two_side_naive_xx(Planc(:,1:iter+1),Qlanc(:,1:iter+1),n,iter+1,&
!!                                    Clanc,Dlanc)
!! Update the new Lanczos vectors
!        Planc(:,1:iter+1) = Clanc
!        Qlanc(:,1:iter+1) = Dlanc
!        deallocate(Clanc,Dlanc)
!
!! Check the level of orthog
!        call zcheck_orthog_frob(Planc(:,1:iter+1),Qlanc(:,1:iter+1),n,iter+1,frob)
!        write(6,*) " "
!        write(6,*) "**************************** " 
!        write(6,*) "After GS with reorthog for iter+1 "
!        write(6,*) " "
!        write(6,*) "^^^^^^^^^^^^^^^^^^"    
!        write(6,*) "# column vecs = ",iter+1,"frob=",frob 
!        write(6,*) "^^^^^^^^^^^^^^^^^^"
!        write(6,*) " "
!
    ! calculated updated dot prod for coeffs
        zterm = wvfdot(backpsi_p,psi_p,allow,n)
        lnorm1 = sqrt(abs(zterm))
        lnorm2 = conjg(zterm)/lnorm1
        write(6,*) " "
        write(6,*) "<backpsi_p,psi_p> =",zterm
        write(6,*) " lnorm =",lnorm1, "lnorm2 =",lnorm2
        write(6,*) " "

     ! Update coeffs    
        c( iter + 1 ) = lnorm1
        b( iter + 1 ) = lnorm2

        zx = 1.0d0 / b( iter + 1 )   
        backzx = 1.0d0 / c( iter + 1 )  

      ! new vectors  
        psi_p( : ) = psi_p( : ) * backzx
        backpsi_p( : ) = backpsi_p( : ) * conjg(zx)

      ! new norms   
        pstot = wvfdot( psi_p, psi_p, allow, n )
        backpstot = wvfdot( backpsi_p, backpsi_p, allow, n )

!         if ( loud ) write ( 6, '(1a6,1x,1f12.8)' ) 'pstot=', pstot
        if ( loud ) write ( 6, * ) 'pstot=', pstot, 'backpstot=', backpstot
    
      !  p_(j-1), q_(j-1)   
        psi_m( : ) = allow( : ) * psi( : )!, now set at beginnning
        psi( : ) = allow( : ) * psi_p( : )

        backpsi_m( : ) = allow( : ) * backpsi( : )!, now set at beginning
        backpsi( : ) = allow( : ) * backpsi_p( : )

 ! Adding the new vectors to Planc, Qlanc
        Planc(:,iter+1) = psi
        Qlanc(:,iter+1)=backpsi

! Reorthogonalize the new iter+1 vectors        
        allocate(Clanc(n,iter+1),Dlanc(n,iter+1))
        Clanc = 0.0
        Dlanc = 0.0
        write(6,*) " "
        write(6,*) " Starting full basis Gram Schmidt, no reorthog"
        write(6,*) " "
        call zgs_class_full_two_side_x(Planc(:,1:iter+1),Qlanc(:,1:iter+1),n,iter+1,&
                                    Clanc,Dlanc)
!        call zgs_blas_full_two_side_x(Planc(:,1:iter+1),Qlanc(:,1:iter+1),n,iter+1,&
!                                    Clanc,Dlanc)
!              call zgs_mod_full_two_side_naive_xx(Planc(:,1:iter+1),Qlanc(:,1:iter+1),n,iter+1,&
!                                    Clanc,Dlanc)
! Update the new Lanczos vectors
        Planc(:,1:iter+1) = Clanc
        Qlanc(:,1:iter+1) = Dlanc
        deallocate(Clanc,Dlanc)

! Check the level of orthog
        call zcheck_orthog_frob(Planc(:,1:iter+1),Qlanc(:,1:iter+1),n,iter+1,frob)
        write(6,*) " "
        write(6,*) "**************************** " 
        write(6,*) "After GS with reorthog for iter+1 "
        write(6,*) " "
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "# column vecs = ",iter+1,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

        !     if ( loud ) write ( 6, '(2x,1a2,1f10.5,3x,1a2,1f10.5)' ) 'a=', a( iter ), 'b=', b( iter + 1 )
         if ( loud ) write ( 6, * ) 'a=', a( iter ), 'b=', b( iter + 1 ),'c=',c(iter + 1)
         cmd = 'opt'
      else
!         fnorm(iter+1) = frob
         cmd = 'end'
      end if ! iter .ne. niter
  case ( 'opt' )
     iter = iter + 1
     if ( iter .le. niter ) then
        cmd = 'clk'
     else
        cmd = 'end'
     end if
  end select
  !
  return
end subroutine nonherm_haydock_reorth


        subroutine zcheck_orthog_frob(A,B,m,k,frob)
!  Calculates the level of orthogonality between two dualspaces of complex matrices manually.
        implicit none
        integer :: i,j,k,m,f
        real(kind=8) :: frob,fterm
        complex(kind=8) :: A(m,k),B(m,k) 
        complex(kind=8),allocatable :: iden(:,:) 

        allocate(iden(k,k))
        iden = 0.0
!        write(6,*) "shape Pj=",shape(A)
!        write(6,*) " "  

!        do i=1,m
!          write(6,*) "A(i,1) =",A(i,1)
!          write(6,*) " "
!        end do
!  Set the identity matrix
        do f = 1,k
          iden(f,f) = 1.0d0
        enddo    
        iden = iden-matmul(conjg(transpose(A)),B)

!  Start calc of frobenius norm    
        fterm = 0.0
        do f=1,k
          do j=1,k
             fterm = fterm + (abs(iden(f,j)))**2
!                 write(6,*) "i=",f,"j=",j,"fterm_ij=",fterm
!                 write(6,*) " "
          enddo
        enddo 
        frob=sqrt(fterm) 
!        write(6,*) " " 
!        write(6,*) "^^^^^^^^^^^^^^^^^^"    
!        write(6,*) "frob=",frob
!        write(6,*) "^^^^^^^^^^^^^^^^^^"
!        write(6,*) " "
!!          fnorm(k) = fterm 
        deallocate(iden)
        return
        end subroutine zcheck_orthog_frob



        subroutine twoside_GS_sing_minprbo(A,B,n,iter,eps1,p,q,pp,qq,max_orth)
!  Does a 2-sided, single step Gram Schmidt for vectors p and q on basis
! A and B, and then returns pp and qq vectors, not normalized.This does minimal partial rebiorthogonalization(minPRBO).
! It checks for the level of biorthogonality first before deciding to orthogonalize against a given vector.
        implicit none
        integer :: i,j,k,m,f,n,iter,iterm,jterm
        complex(kind=kind(1.0d0)) :: A(n,iter),B(n,iter),p(n),q(n),pp(n),qq(n) 
        real( kind = kind( 1.0d0 ) ) :: eps,eps1,term_re,term_im,frob,fterm,tol,max_orth
        integer :: orthtot,poterm,qoterm,chk,oterm
        complex( kind = kind( 1.0d0 ) ) :: zterm,term1,lnorm1,lnorm2,vp(n),vq(n),wp(n),wq(n)
        integer,allocatable :: porthopt(:),qorthopt(:),porthindx(:),qorthindx(:)
        complex( kind = kind( 1.0d0 ) ),allocatable :: iden(:,:),Pj(:,:),Qj(:,:),iden1(:,:),&
                                                 iden2(:,:),pw_jk(:),qw_jk(:)
!            write(6,*)" "
!            write(6,*) "iter =",iter
!            write(6,*) " "
            pp = 0.0
            qq = 0.0
            allocate(Pj(n,iter),Qj(n,iter),pw_jk(iter),qw_jk(iter))
            Pj = 0.0
            Qj = 0.0
            pw_jk = 0.0
            qw_jk = 0.0
            Pj = A(:,1:iter)
            Qj = B(:,1:iter)
!            write(6,*) "shape of Pj/Qj =",shape(Qj),"shape Planc =",shape(A)
!            write(6,*)" "
            
        ! Calculating local orthogonality, first normalizing p and q
            vp = 0.0
            vq = 0.0    
            term1 = dot_product(q,p)
!            term1 = wvfdot(backpsi_p,psi_p,allow,n)
            lnorm1 = sqrt(abs(term1))
            lnorm2 = conjg(term1)/lnorm1
            vp = (1.0/conjg(lnorm2))*p
            vq = (1.0/lnorm1)*q
            pw_jk = matmul(conjg(transpose(Qj)),vp)
            qw_jk = matmul(conjg(transpose(Pj)),vq)

            ! estimate of level of biorthogonality
            max_orth = max(maxval(abs(pw_jk)),maxval(abs(qw_jk)))

!            write(6,*) " "
!            write(6,*)"orthog,<p_j+1,q_k>, max",maxval(abs(pw_jk))
!            write(6,*) " "
!            write(6,*)pw_jk
!            write(6,*) " "
!            write(6,*)"orthog,<q_j+1,p_k>, max",maxval(abs(qw_jk))
!            write(6,*) " "
!            write(6,*) qw_jk
!            write(6,*)" "

            deallocate(Pj,Qj)

         ! Check for partial orthogonality
           allocate(porthopt(iter),qorthopt(iter))
           porthopt = 0
           qorthopt = 0
           poterm = 0
           qoterm = 0
         
         ! P's
           do i=1,iter
             zterm = pw_jk(i)
             term_re = abs(real(zterm))
             term_im = abs(aimag(zterm))
             if ((term_re .ge. eps1).or.(term_im .ge. eps1)) then
!               write(6,*) " "
!               write(6,*) "vector",i,"lost orth, pw_jk=",zterm,eps1
!               write(6,*) " "
               porthopt(i) = 1
             endif
           enddo
           poterm = sum(porthopt)
         ! Q's
           do i=1,iter
             zterm = qw_jk(i)
             term_re = abs(real(zterm))
             term_im = abs(aimag(zterm))
             if ((term_re .ge. eps1).or.(term_im .ge. eps1)) then
!               write(6,*) " "
!               write(6,*) "vector",i,"lost orth, qw_jk=",zterm,eps1
!               write(6,*) " "
               qorthopt(i) = 1
             endif
           enddo
           qoterm = sum(qorthopt)
           write(6,*)"---# vecs that lost biortho",poterm,qoterm,"---"
           write(6,*) " "

        ! Get the index of these vectors
           if (poterm .ne. qoterm) then
              oterm = max(poterm,qoterm)

           ! Make the indices the same if the # or orthog is diff
              if(poterm .ge. qoterm) then
                 qorthopt=qorthopt+(porthopt-qorthopt)
              else
                 porthopt=porthopt+(qorthopt-porthopt)
              endif
              poterm=oterm
              qoterm=oterm
           endif 
           
           allocate(porthindx(poterm),qorthindx(qoterm))
           porthindx = 0
           qorthindx = 0

    ! P's
           iterm = 1
           do i=1,iter
             if(porthopt(i) .eq. 1) then
               porthindx(iterm) = i
               iterm = iterm + 1
             endif
           enddo  

    ! Q's
           iterm = 1
           do i=1,iter
             if(qorthopt(i) .eq. 1) then
               qorthindx(iterm) = i
               iterm = iterm + 1
             endif
           enddo  
               
         deallocate(porthopt,qorthopt)  
!
    ! Construct the Pj and Qj
         if ((poterm .ge. 1) .or. (qoterm .ge. 1))then
!           write(6,*) "doing partial 2-sided GS as biorthog failed"
!           write(6,*) " poterm =",poterm,"qoterm =",qoterm
!           write(6,*) " "
           allocate(Pj(n,poterm),Qj(n,qoterm),iden1(n,n),iden2(n,n))
!           write(6,*) "allocated Pj, Qj"
!           write(6,*) "allocated iden1,iden2"
           iden1 = 0.0
           iden2 = 0.0
           Pj = 0.0
           Qj = 0.0
           do i=1,poterm
             iterm = porthindx(i)
             jterm = qorthindx(i)
!             write(6,*) "iterm =",iterm,"jterm =",jterm
!             write(6,*) " "
             Pj(:,i) = A(:,iterm)
             Qj(:,i) = B(:,jterm)
           enddo
    ! Set the identity matrix
           do f=1,n
              iden1(f,f) = 1.0d0
              iden2(f,f) = 1.0d0
           enddo
!           write(6,*) " Identity matrix is set"
!           write(6,*) " "
           
           iden1 = iden1-matmul(Qj,conjg(transpose(Pj)))
           iden2 = iden2-matmul(Pj,conjg(transpose(Qj)))

    ! re biorthogonalizing
           qq = matmul(iden1,q) 
           pp = matmul(iden2,p)
         !  write(6,*) " "
         !  write(6,*) "pp =",pp
         !  write(6,*) " "
         !  write(6,*) "----- Comparing the level of orthog after rebiortho ----"
         !  write(6,*) " "

         !   pw_jk = matmul(conjg(transpose(A(:,1:iter))),pp)
         !   qw_jk = matmul(conjg(transpose(B(:,1:iter))),qq)

         !   write(6,*) " "
         !   write(6,*)"orthog,<p_j+1,q_k>, max",maxval(abs(pw_jk))
         !   write(6,*) " "
         !   write(6,*)pw_jk
         !   write(6,*) " "
         !   write(6,*)"orthog,<q_j+1,p_k>, max",maxval(abs(qw_jk))
         !   write(6,*) " "
         !   write(6,*) qw_jk
         !   write(6,*)" "
           deallocate(iden1,iden2,Pj,Qj)  
         else
!           write(6,*) " "
!           write(6,*) " No orthog was needed "
           pp = p
           qq = q 
         endif
         deallocate(porthindx,qorthindx) 
         deallocate(pw_jk,qw_jk)
         
!         write(6,*) " "
!         write(6,*) "Done with re-biorthogonalization"
         return
         end subroutine twoside_GS_sing_minprbo

        subroutine zgs_class_full_two_side_x(A,B,m,n,C,D)
  ! Two sided full, classical Gram-Schmidt for complex matrices
        implicit none
        integer :: g,i,j,k,l,m,n,f
        real(kind=8) :: temp1,frob,lterm,term
        complex(kind=8) :: zterm, lnorm1,lnorm2
        complex(kind=8) :: A(m,n),B(m,n),C(m,n),D(m,n)
        complex(kind=8) :: a_j(m),b_j(m),q_j(m),p_j(m)
        complex(kind=8),allocatable ::iden1(:,:),iden2(:,:),&
                                      Pj(:,:),Qj(:,:)
        
        allocate(iden1(m,m),iden2(m,m))
        C = 0.0
        D = 0.0
        a_j = 0.0
        b_j = 0.0
        q_j = 0.0
        p_j = 0.0
         
        write(6,*) " "
        write(6,*) "Shape of A/B is",shape(A)
        write(6,*) " "
        write(6,*) " Two sided full Gram-Schmidt, one step"

! Checking the level of orthog with Frobenius norm before two-sided GS
        call zcheck_orthog_frob(A,B,m,n,frob)
        write(6,*) " "
        write(6,*) "!!!!!!!!!!!!!!!!!!!! " 
        write(6,*) "Before Gram Schmidt "
        write(6,*) " "
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "# column vecs = ",n,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

!------- Starting the two sided Gram-Schmidt  -------
        do g=1,n
!           write(6,*) " "
!           write(6,*) "********* Gram-Schmidt iteration=",g,"*********"
!           write(6,*) " "
           a_j = A(:,g)  
           b_j = B(:,g)

 
!  Set the identity matrix
           iden1 = 0.0
           iden2 = 0.0
           do f = 1,m
             iden1(f,f) = 1.0d0
             iden2(f,f) = 1.0d0
           enddo  
           if (g .ne. 1) then

              allocate(Pj(m,g-1),Qj(m,g-1))
!              write(6,*) " "
!              write(6,*) "g = ",g
!              write(6,*) " Done allocating Pj and Qj"
              Pj = 0.0
              Qj = 0.0
              Pj = C(:,1:g-1) !originally A, B
              Qj = D(:,1:g-1)
! For Hermitian case, works well  
              iden1 = iden1-matmul(Qj,conjg(transpose(Pj)))
              iden2 = iden2-matmul(Pj,conjg(transpose(Qj)))

! Gram Schmidt step
              q_j = matmul(iden1,b_j)
              p_j = matmul(iden2,a_j)
              deallocate(Pj,Qj)
            else
!              write(6,*) " "
!              write(6,*) " No need for re-biorthog "
              q_j = b_j
              p_j = a_j
            endif
 
!!  New norm
            zterm = dot_product(q_j,p_j) ! following netlib Alg 7.13
            lnorm1 = sqrt(abs(zterm))  ! following yambo paper
            lnorm2 = conjg(zterm)/lnorm1
!            write(6,*) "lnorm2 =",lnorm2
!            write(6,*) ""
!            write(6,*) "g =",g,"zterm = ",zterm
           

! Normalizing and updating Lanczos vectors.
! Should keep the last vector added un normalized, since this will need to be normalized
! to calculate the coefficients. Or does it not matter if it is also normalized?
            q_j = (1.0/lnorm1)*q_j
            p_j = (1.0/conjg(lnorm2))*p_j


            C(:,g) = p_j
            D(:,g) = q_j
        enddo ! loop over Lanczos columns for two-sided CGS
        deallocate(iden1,iden2)

! Checking the level of orthog with Frobenius norm after two-sided GS
        call zcheck_orthog_frob(C,D,m,n,frob)
        write(6,*) " " 
        write(6,*) "After Gram Schmidt "
        write(6,*) " "
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "# column vecs = ",n,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

        return
        end subroutine zgs_class_full_two_side_x


        subroutine zgs_class_two_side_pro(A,B,k,m,n,op,oq,p,q,pp,qq,&
                                          indxp,indxq,frob)
! Two sided Gram-Schmidt for partial reorthog, for complex matrices
        implicit none
        integer :: g,i,j,k,l,m,n,f,op,oq,indxp(op),indxq(oq),&
                   omin,iterm,jterm
        real(kind=8) :: temp1,frob,lterm,term
        complex(kind=8) :: zterm, lnorm1,lnorm2
        complex(kind=8) :: A(m,n),B(m,n),C(m,n),D(m,n)
        complex(kind=8) ::a_j(m),b_j(m),q_j(m),p_j(m),p(m),q(m),&
                          pp(m),qq(m)
        complex(kind=8),allocatable ::iden1(:,:),iden2(:,:),&
                                      Pj(:,:),Qj(:,:)
        
        allocate(iden1(m,m),iden2(m,m))
        qq = 0.0
        pp = 0.0
    
!------- Starting the two sided Gram-Schmidt  -------
 
!  Set the identity matrix
        iden1 = 0.0
        iden2 = 0.0
        do f = 1,m
          iden1(f,f) = 1.0d0
          iden2(f,f) = 1.0d0
        enddo  

!  Make sure they arrays match     
        write(6,*) " "
        write(6,*) "op =",op,"oq =",oq
        write(6,*) " "      
        if (op .ne. oq) then 
          omin=min(op,oq)
        else
          omin = op
        endif
        allocate(Pj(m,omin),Qj(m,omin))
        Pj = 0.0
        Qj = 0.0

! Form the matrix that holds the selected vectors
 
!  P's, Q's
        do i=1,omin
          iterm = indxp(i)
          jterm = indxq(i)
          Pj(:,i) = A(:,iterm)
          Qj(:,i) = B(:,jterm)
        enddo

! For Hermitian case, works well  
        iden1 = iden1-matmul(Qj,conjg(transpose(Pj)))
        iden2 = iden2-matmul(Pj,conjg(transpose(Qj)))

! re biorthogonalizing
        qq = matmul(iden1,q)
        pp = matmul(iden2,p)

        deallocate(iden1,iden2,Pj,Qj)
          
!!  New norm
!        zterm = dot_product(qq,pp) ! following netlib Alg 7.13
!        lnorm1 = sqrt(abs(zterm))  ! following yambo paper
!        lnorm2 = conjg(zterm)/lnorm1
!        write(6,*) "lnorm2 =",lnorm2
!        write(6,*) ""
!
!! Normalizing and updating Lanczos vectors
!        qq = (1.0/lnorm1)*qq
!        pp = (1.0/conjg(lnorm2))*pp
!
!! Adding the new vecs to the Lanczos basis
!        A(:,k+1) = pp
!        B(:,k+1) = qq

!! Checking the level of orthog with Frobenius norm after two-sided GS
!        call zcheck_orthog_frob(A(:,1:k+1),B(:,1:k+1),m,k+1,frob)
!        write(6,*) " " 
!        write(6,*) "^^^^^^^^^^^^^^^^^^"    
!        write(6,*) "frob=",frob 
!        write(6,*) "^^^^^^^^^^^^^^^^^^"
!        write(6,*) " "

        return
        end subroutine zgs_class_two_side_pro


        subroutine twoside_GS_sing_x(A,B,n,iter,eps1,p,q,pp,qq)
!  Does a 2-sided, single step Gram Schmidt for vectors p and q on basis
! A and B, and then returns r and s vectors, not normalized.This does a full reorthogonalization.

        implicit none
        integer :: i,j,k,m,f,n,iter,iterm,jterm
        complex(kind=kind(1.0d0)) :: A(n,iter),B(n,iter),p(n),q(n),pp(n),qq(n) 
        real( kind = kind( 1.0d0 ) ) :: eps,eps1,term_re,term_im,frob,fterm,tol
        integer :: orthtot,poterm,qoterm,chk,oterm
        complex( kind = kind( 1.0d0 ) ) :: zterm,term1,lnorm1,lnorm2,vp(n),vq(n),wp(n),wq(n)
        integer,allocatable :: porthopt(:),qorthopt(:),porthindx(:),qorthindx(:)
        complex( kind = kind( 1.0d0 ) ),allocatable :: iden(:,:),Pj(:,:),Qj(:,:),iden1(:,:),&
                                                 iden2(:,:),pw_jk(:),qw_jk(:)
            write(6,*)" "
            write(6,*) "iter =",iter
            write(6,*) " "
            pp = 0.0
            qq = 0.0
            write(6,*) "shape of Planc =",shape(A),"shape Qlanc =",shape(B)
            write(6,*)" "
            

    ! Set the identity matrix
           allocate(iden1(n,n),iden2(n,n))
           do f=1,n
              iden1(f,f) = 1.0d0
              iden2(f,f) = 1.0d0
           enddo
           write(6,*) " Identity matrix is set"
           write(6,*) " "
           
           iden1 = iden1-matmul(B,conjg(transpose(A)))
           iden2 = iden2-matmul(A,conjg(transpose(B)))

    ! re biorthogonalizing
           qq = matmul(iden1,q) 
           pp = matmul(iden2,p)
           deallocate(iden1,iden2)  
         
         return
         end subroutine twoside_GS_sing_x


        subroutine twoside_GS_naive_sing_x(A,B,n,iter,eps1,p,q,pp,qq)
!  Does a 2-sided, single step Gram Schmidt for vectors p and q on basis
! A and B, and then returns r and s vectors, not normalized.
! It does it without matrix multiplication, the naive way.

        implicit none
        integer :: i,j,k,m,f,n,iter,iterm,jterm
        complex(kind=kind(1.0d0)) :: A(n,iter),B(n,iter),p(n),q(n),pp(n),qq(n) 
        real( kind = kind( 1.0d0 ) ) :: eps,eps1,term_re,term_im,frob,fterm,tol
        integer :: orthtot,poterm,qoterm,chk,oterm
        complex( kind = kind( 1.0d0 ) ) :: zterm,pterm,qterm,term1,lnorm1,lnorm2,vp(n),vq(n),wp(n),wq(n),&
                                           vp_j(n),vq_j(n),p_k(n),q_k(n)
        integer,allocatable :: porthopt(:),qorthopt(:),porthindx(:),qorthindx(:)
        complex( kind = kind( 1.0d0 ) ),allocatable :: iden(:,:),Pj(:,:),Qj(:,:),iden1(:,:),&
                                                 iden2(:,:),pw_jk(:),qw_jk(:)
  ! Lapack routines
        complex( kind = kind( 1.0d0 ) ):: ZDOTC
            write(6,*)" "
            write(6,*) "iter =",iter
            write(6,*) " "
            pp = 0.0
            qq = 0.0
            vp_j = 0.0
            vq_j = 0.0
            p_k = 0.0
            q_k = 0.0
            write(6,*) "shape of Planc =",shape(A),"shape Qlanc =",shape(B)
            write(6,*)" "
            


           if (iter .ne. 1) then
! 1st orthogonalization
                   vp_j = p
                   vq_j = q

                   do j=1,iter
                     p_k = A(:,j) ! initially A, B
                     q_k = B(:,j)
!                     pterm = dot_product(vp_j,q_k)
!                     qterm = dot_product(vq_j,p_k)
                     pterm = ZDOTC(n,vp_j,1,q_k,1)
                     qterm = ZDOTC(n,vq_j,1,p_k,1)
                     vp_j = vp_j - pterm*p_k
                     vq_j = vq_j - conjg(qterm)*q_k
                   enddo  
           else
             vp_j = p
             vq_j = q
           endif
           pp = vp_j
           qq = vq_j
         return
         end subroutine twoside_GS_naive_sing_x

        subroutine zgs_class_full_two_side_xx(A,B,m,n,C,D)
  ! Two sided full, classical Gram-Schmidt for complex matrices,
  ! with reorthogonalization.

        implicit none
        integer :: g,i,j,k,l,m,n,f
        real(kind=8) :: temp1,frob,lterm,term
        complex(kind=8) :: zterm, lnorm1,lnorm2
        complex(kind=8) :: A(m,n),B(m,n),C(m,n),D(m,n)
        complex(kind=8) :: a_j(m),b_j(m),q_j(m),p_j(m),qq_j(m),pp_j(m)
        complex(kind=8),allocatable ::iden1(:,:),iden2(:,:),&
                                      Pj(:,:),Qj(:,:)
        
        allocate(iden1(m,m),iden2(m,m))
        C = 0.0
        D = 0.0
        a_j = 0.0
        b_j = 0.0
        q_j = 0.0
        p_j = 0.0
         
        write(6,*) " "
        write(6,*) "Shape of A/B is",shape(A)
        write(6,*) " "
        write(6,*) " Two sided full Gram-Schmidt, with reorthog"

! Checking the level of orthog with Frobenius norm before two-sided GS
        call zcheck_orthog_frob(A,B,m,n,frob)
        write(6,*) " "
        write(6,*) "!!!!!!!!!!!!!!!!!!!! " 
        write(6,*) "Before Gram Schmidt "
        write(6,*) " "
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "# column vecs = ",n,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

!------- Starting the two sided Gram-Schmidt  -------
        do g=1,n
!           write(6,*) " "
!           write(6,*) "********* Gram-Schmidt iteration=",g,"*********"
!           write(6,*) " "
           a_j = A(:,g)
           b_j = B(:,g)

 
!  Set the identity matrix
           iden1 = 0.0
           iden2 = 0.0
           do f = 1,m
             iden1(f,f) = 1.0d0
             iden2(f,f) = 1.0d0
           enddo  
           if (g .ne. 1) then

              allocate(Pj(m,g-1),Qj(m,g-1))
!              write(6,*) " "
!              write(6,*) "g = ",g
!              write(6,*) " Done allocating Pj and Qj"
              Pj = 0.0
              Qj = 0.0
              Pj = A(:,1:g-1)
              Qj = B(:,1:g-1)
! For Hermitian case, works well  
              iden1 = iden1-matmul(Qj,conjg(transpose(Pj)))
              iden2 = iden2-matmul(Pj,conjg(transpose(Qj)))

! Gram Schmidt step
              q_j = matmul(iden1,b_j)
              p_j = matmul(iden2,a_j)

! Reorthogonalization step
              qq_j = q_j
              pp_j = p_j
              q_j = matmul(iden1,qq_j)
              p_j = matmul(iden2,pp_j)
              deallocate(Pj,Qj)
            else
!              write(6,*) " "
!              write(6,*) " No need for re-biorthog "
              q_j = b_j
              p_j = a_j
            endif
 
!!  New norm
            zterm = dot_product(q_j,p_j) ! following netlib Alg 7.13
            lnorm1 = sqrt(abs(zterm))  ! following yambo paper
            lnorm2 = conjg(zterm)/lnorm1
!            write(6,*) "lnorm2 =",lnorm2
!            write(6,*) ""

! Normalizing and updating Lanczos vectors.
! Should keep the last vector added un normalized, since this will need to be normalized
! to calculate the coefficients. Or does it not matter if it is also normalized?
!            if (g .ne. n) then
            q_j = (1.0/lnorm1)*q_j
            p_j = (1.0/conjg(lnorm2))*p_j
!            endif

            C(:,g) = p_j
            D(:,g) = q_j
        enddo ! loop over Lanczos columns for two-sided CGS
        deallocate(iden1,iden2)

! Checking the level of orthog with Frobenius norm after two-sided GS
        call zcheck_orthog_frob(C,D,m,n,frob)
        write(6,*) " " 
        write(6,*) "After Gram Schmidt with reorthog"
        write(6,*) " "
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "# column vecs = ",n,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

        return
        end subroutine zgs_class_full_two_side_xx


        subroutine twoside_GS_sing_xx(A,B,n,iter,eps1,p,q,pp,qq)
!  Does a 2-sided, single step Gram Schmidt for vectors p and q on basis
! A and B, and then returns r and s vectors, not normalized.This does a full reorthogonalization.
! And also does a second step of Gram Schmidt reorthog
        implicit none
        integer :: i,j,k,m,f,n,iter,iterm,jterm
        complex(kind=kind(1.0d0)) :: A(n,iter),B(n,iter),p(n),q(n),pp(n),qq(n),p2(n),q2(n) 
        real( kind = kind( 1.0d0 ) ) :: eps,eps1,term_re,term_im,frob,fterm,tol
        integer :: orthtot,poterm,qoterm,chk,oterm
        complex( kind = kind( 1.0d0 ) ) :: zterm,term1,lnorm1,lnorm2,vp(n),vq(n),wp(n),wq(n)
        integer,allocatable :: porthopt(:),qorthopt(:),porthindx(:),qorthindx(:)
        complex( kind = kind( 1.0d0 ) ),allocatable :: iden(:,:),Pj(:,:),Qj(:,:),iden1(:,:),&
                                                 iden2(:,:),pw_jk(:),qw_jk(:)
            write(6,*)" "
            write(6,*) "iter =",iter
            write(6,*) " "
            pp = 0.0
            qq = 0.0
            p2 = 0.0
            q2 = 0.0
            write(6,*) "shape of Planc =",shape(A),"shape Qlanc =",shape(B)
            write(6,*)" "
            
        ! First normalizing p and q
            vp = 0.0
            vq = 0.0    
            term1 = dot_product(q,p)
!            term1 = wvfdot(backpsi_p,psi_p,allow,n)
            lnorm1 = sqrt(abs(term1))
            lnorm2 = conjg(term1)/lnorm1
            vp = (1.0/conjg(lnorm2))*p
            vq = (1.0/lnorm1)*q

    ! Set the identity matrix
           allocate(iden1(n,n),iden2(n,n))
           do f=1,n
              iden1(f,f) = 1.0d0
              iden2(f,f) = 1.0d0
           enddo
           write(6,*) " Identity matrix is set"
           write(6,*) " "
           
           iden1 = iden1-matmul(B,conjg(transpose(A)))
           iden2 = iden2-matmul(A,conjg(transpose(B)))

    ! Gram Schmidt step
           qq = matmul(iden1,q) 
           pp = matmul(iden2,p)

    ! Doing reorthog of Gram Schmidt
           q2 = qq
           p2 = pp
           qq = matmul(iden1,q2) 
           pp = matmul(iden2,p2)

           deallocate(iden1,iden2)  
         
         return
         end subroutine twoside_GS_sing_xx
         
         subroutine gram_schmidt_class_reorth(A,m,fnorm)
! Classical Gram-Schmidt from paper: "The loss of orthogonality in the
! Gram-Schmidt orthogonalization process"
!        Use lapack_interfaces
!        Use lapack_interfaces, Only: DGEEV
!        Use lapack_precision, Only: dp
        implicit none
        integer :: i,j,k,l,m,f
        real(kind=8) :: b0,temp1,tol,fnorm(m),fterm,lterm,term
        real(kind=8) :: u0,lnorm,DLANGE
        real(kind=8) :: rand_num(m)
        real,parameter :: pi=3.141592654d0
        complex,parameter :: ii=complex(0.d0,1.d0)
        
!        real(kind=8) :: alpha(mlanc),beta(mlanc)
        real(kind=8):: vm1(m),a_j(m),s(m),q_j(m),v_j(m),w_j(m),q_k(m),&
                       t_j(m)
        
        real(kind=8) :: A(m,m)

!  Matrices for calculating biorthogonality
        real(kind=8) :: Planc(m,m) ! since doing orthogonalization
!       from 1 to j-1
        real(kind=8), allocatable :: iden(:,:),Pj(:,:)
        integer,parameter :: sedsz=33
        integer :: seed(sedsz) 

! Lapack parameters for getting Eigenvalues (DGEEV example program in
! Fortran, www.intel.com

!        integer, parameter :: LDA=m, N=m, LDVL=m, LDVR=m,LWMAX=1000
        integer :: LDA, N, LDVL, LDVR,LWMAX,INFO,LWORK
        integer, parameter :: nb = 64
!        real(kind=8) :: VL(LDVL,m), VR(LDVR,m),WR(m), WI(m),WORK(LMAX),&
!                        B(LDA,m)
        real(kind=8),allocatable :: VL(:,:), VR(:,:),WR(:), WI(:),&
                                    WORK(:),B(:,:),C(:,:)
        real(kind=8) :: dummy(1,1)

        
!  Printing out matrix elements
!        do i=1,m
!          do j=1,m
!             write(6,*) ""
!             write(6,*) "--- i = ",i,"j = ",j,"A(i,j) = ",A(i,j)
!             write(6,*) ""
!          enddo
!        enddo
        write(6,*) "||A||_F =",norm2(A)
        write(6,*) " "
        write(6,*) " Doing Gram-Schmidt with reorthogonalization" 

!   Initializing arrays to 0
!        alpha=0.0
!        beta=0.0
!        c=0.0
        vm1=0.0
        a_j=0.0
        q_k=0.0
        q_j=0.0
        v_j=0.0
        w_j=0.0
        t_j=0.0
        Planc = 0.0

!  Generating the test vector , v_1

        call random_seed(put=seed(1:sedsz))
        call random_number(rand_num)
        
        a_j=A(:,1)
        u0=dot_product(a_j,a_j)
        u0=sqrt(u0)
        write(6,*) "--- The Euclidean norm =",u0
        q_j=a_j/u0
!        do i=1,m
!          write(6,*)""
!          write(6,*) "q_j(i) = ",q_j(i)
!          write(6,*) ""
!        enddo
        Planc(:,1) = q_j


!----------  Starting the Lanczos iteration ------------------------
        do k=1,m
          write(6,*) " "
          write(6,*) "******** Gram-Schmidt iteration =",k,"***********"

! Traditional formulation : "Loss of orthogonality in the Gram-Schmidt
! Orthogonalization process"

!  1st orthogonlization
          a_j = A(:,k)
          v_j = a_j

          do j=1,k-1
             q_k = Planc(:,j)
             term = dot_product(a_j,q_k)
!          write(6,*) ""
!          write(6,*) "j =",j,"alpha(k,j) =",term
!          write(6,*) ""
             do i=1,m
!            write(6,*) ""
!            write(6,*) " --- i = ",i,"------"
!            write(6,*) "before updating, t(i)=",t(i)
               v_j(i)=v_j(i)-term*q_k(i) 
!            write(6,*) ""
             enddo
           enddo
!!!  New norm
            lnorm = dot_product(v_j,v_j)
            lnorm = sqrt(lnorm)
            !            write(6,*) "norm w_j=",lnorm
!            write(6,*) " "  
        
            do i=1,m
              q_j(i)=(1/lnorm)*w_j(i)
            enddo

!!  Add new vectors 
            Planc(:,k) = q_j


!  Checking biorthogonality with the frobenius norm and spectral norm
            allocate(iden(k,k),Pj(m,k))
            iden = 0.0
            Pj = 0.0
            Pj = Planc(:,1:k)
!            write(6,*) "shape Pj=",shape(Pj)
!            write(6,*) " "  
!  Set the identity matrix
            do f = 1,k
              iden(f,f) = 1.0d0
            enddo    
            iden = iden-matmul(transpose(Pj),Pj)

!  Start calc of frobenius norm    
            fterm = 0.0
            do f=1,k
              do j=1,k
                 fterm = fterm + (abs(iden(f,j)))**2
!                 write(6,*) "i=",f,"j=",j,"fterm_ij=",fterm
!                 write(6,*) " "
              enddo
            enddo 
            fterm=sqrt(fterm) 


            write(6,*) " " 
            write(6,*) "^^^^^^^^^^^^^^^^^^"    
            write(6,*) "frob=",norm2(real(iden,kind=8)) !,"L2 nrm=",term
            write(6,*) "^^^^^^^^^^^^^^^^^^"
            write(6,*) " "
            fnorm(k) = fterm 

! deallocate necessary arrays
            deallocate(iden,Pj)
!            deallocate(VL,VR,WR,WI,WORK,B)

        enddo ! END OF MGS LOOP
        return
        end subroutine gram_schmidt_class_reorth

        subroutine zgs_mod_full_two_side_naive_xx(A,B,m,n,C,D)
  ! Two sided full, modified Gram-Schmidt for complex matrices,
  ! with reorthogonalization. Does not use matrix multiplication to carry out the Gram-Schmidt,
  ! hence dont need to create these large matrices.

        implicit none
        integer :: g,i,j,k,l,m,n,f
        real(kind=8) :: temp1,frob,lterm,term
        complex(kind=8) :: zterm,zterm1,zterm2,pterm,qterm,lnorm1,lnorm2
        complex(kind=8) :: A(m,n),B(m,n),C(m,n),D(m,n)
        complex(kind=8) :: a_j(m),b_j(m),q_j(m),p_j(m),qq_j(m),pp_j(m),&
                          q_k(m),p_k(m),vp_j(m),vq_j(m),wp_j(m),wq_j(m),&
                          tp_j(m),tq_j(m),vpp_j(m),vqq_j(m)
        complex(kind=8),allocatable ::iden1(:,:),iden2(:,:),&
                                      Pj(:,:),Qj(:,:)
        
  ! Lapack routines
        complex( kind = kind( 1.0d0 ) ):: ZDOTC

        allocate(iden1(m,m),iden2(m,m))
        C = 0.0
        D = 0.0
        a_j = 0.0
        b_j = 0.0
        q_j = 0.0
        p_j = 0.0
        q_k = 0.0
        p_k = 0.0
        vp_j = 0.0
        vq_j = 0.0
        wp_j = 0.0
        wq_j = 0.0
        tp_j = 0.0
        tq_j = 0.0

        write(6,*) " "
        write(6,*) "Shape of A/B is",shape(A)
        write(6,*) " "
        write(6,*) " Two sided full Gram-Schmidt, with reorthog"

! Checking the level of orthog with Frobenius norm before two-sided GS
        call zcheck_orthog_frob(A,B,m,n,frob)
        write(6,*) " "
        write(6,*) "!!!!!!!!!!!!!!!!!!!! " 
        write(6,*) "Before Gram Schmidt "
        write(6,*) " "
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "# column vecs = ",n,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

!------- Starting the two sided Gram-Schmidt  -------
        do g=1,n
!           write(6,*) " "
!           write(6,*) "********* Gram-Schmidt iteration=",g,"*********"
!           write(6,*) " "
           a_j = A(:,g)
           b_j = B(:,g)
           
           if (g .ne. 1) then
! 1st orthogonalization
                   vp_j = a_j
                   vq_j = b_j

                   do j=1,g-1
                     p_k = C(:,j)
                     q_k = D(:,j)
!                     pterm = dot_product(vp_j,q_k)
!                     qterm = dot_product(vq_j,p_k)
                     pterm = ZDOTC(m,vp_j,1,q_k,1)
                     qterm = ZDOTC(m,vq_j,1,p_k,1)
                     vp_j = vp_j - pterm*p_k
                     vq_j = vq_j - qterm*q_k
                   enddo 

! 2nd orthogonalization       
                   vpp_j = vp_j
                   vqq_j = vq_j
                   do j=1,g-1
                      p_k = C(:,j)
                      q_k = D(:,j)
!                      pterm = dot_product(vpp_j,q_k)
!                      qterm = dot_product(vqq_j,p_k)
                      pterm = ZDOTC(m,vpp_j,1,q_k,1)
                      qterm = ZDOTC(m,vqq_j,1,p_k,1)
                      vpp_j = vpp_j - pterm*p_k
                      vqq_j = vqq_j - qterm*q_k
                   enddo
                   vp_j = vpp_j
                   vq_j = vqq_j
           
           else
             vp_j = a_j
             vq_j = b_j
           endif
!!  New norm
!            zterm = dot_product(vq_j,vp_j) ! following netlib Alg 7.13
            zterm = ZDOTC(m,vq_j,1,vp_j,1) ! following netlib Alg 7.13
            lnorm1 = sqrt(abs(zterm))  ! following yambo paper
            lnorm2 = conjg(zterm)/lnorm1
!            write(6,*) "lnorm2 =",lnorm2
!            write(6,*) ""

! Normalizing and updating Lanczos vectors
            q_j = (1.0/lnorm1)*vq_j
            p_j = (1.0/conjg(lnorm2))*vp_j


            C(:,g) = p_j
            D(:,g) = q_j
        enddo ! loop over Lanczos columns for two-sided CGS

! Checking the level of orthog with Frobenius norm after two-sided GS
        call zcheck_orthog_frob(C,D,m,n,frob)
        write(6,*) " " 
        write(6,*) "After Gram Schmidt with reorthog"
        write(6,*) " "
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "# column vecs = ",n,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

        return
        end subroutine zgs_mod_full_two_side_naive_xx


        subroutine zgs_mod_full_two_side_naive_x(A,B,m,n,C,D)
  ! Two sided full, modified Gram-Schmidt for complex matrices,
  ! without reorthogonalization. Does not use matrix multiplication to carry out the Gram-Schmidt,
  ! hence dont need to create these large matrices.

        implicit none
        integer :: g,i,j,k,l,m,n,f
        real(kind=8) :: temp1,frob,lterm,term
        complex(kind=8) :: zterm,zterm1,zterm2,pterm,qterm,lnorm1,lnorm2
        complex(kind=8) :: A(m,n),B(m,n),C(m,n),D(m,n)
        complex(kind=8) :: a_j(m),b_j(m),q_j(m),p_j(m),qq_j(m),pp_j(m),&
                          q_k(m),p_k(m),vp_j(m),vq_j(m),wp_j(m),wq_j(m),&
                          tp_j(m),tq_j(m),vpp_j(m),vqq_j(m)
        complex(kind=8),allocatable ::iden1(:,:),iden2(:,:),&
                                      Pj(:,:),Qj(:,:)
        
  ! Lapack routines
        complex( kind = kind( 1.0d0 ) ):: ZDOTC

        allocate(iden1(m,m),iden2(m,m))
        C = 0.0
        D = 0.0
        a_j = 0.0
        b_j = 0.0
        q_j = 0.0
        p_j = 0.0
        q_k = 0.0
        p_k = 0.0
        vp_j = 0.0
        vq_j = 0.0
        wp_j = 0.0
        wq_j = 0.0
        tp_j = 0.0
        tq_j = 0.0

        write(6,*) " "
        write(6,*) "Shape of A/B is",shape(A)
        write(6,*) " "
        write(6,*) " Two sided full Gram-Schmidt, with reorthog"

! Checking the level of orthog with Frobenius norm before two-sided GS
        call zcheck_orthog_frob(A,B,m,n,frob)
        write(6,*) " "
        write(6,*) "!!!!!!!!!!!!!!!!!!!! " 
        write(6,*) "Before Gram Schmidt "
        write(6,*) " "
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "# column vecs = ",n,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

!------- Starting the two sided Gram-Schmidt  -------
        do g=1,n
!           write(6,*) " "
!           write(6,*) "********* Gram-Schmidt iteration=",g,"*********"
!           write(6,*) " "
           a_j = A(:,g)
           b_j = B(:,g)
           
           if (g .ne. 1) then
! 1st orthogonalization
                   vp_j = a_j
                   vq_j = b_j

                   do j=1,g-1
                   p_k = C(:,j) ! initially A, B
                     q_k = D(:,j)
!                     pterm = dot_product(vp_j,q_k)
!                     qterm = dot_product(vq_j,p_k)
                     pterm = ZDOTC(m,vp_j,1,q_k,1)
                     qterm = ZDOTC(m,vq_j,1,p_k,1)
                     vp_j = vp_j - pterm*p_k
                     vq_j = vq_j - conjg(qterm)*q_k
                   enddo  
           else
             vp_j = a_j
             vq_j = b_j
           endif
!!  New norm
!            zterm = dot_product(vq_j,vp_j) ! following netlib Alg 7.13
            zterm = ZDOTC(m,vq_j,1,vp_j,1) ! following netlib Alg 7.13
            lnorm1 = sqrt(abs(zterm))  ! following yambo paper
            lnorm2 = conjg(zterm)/lnorm1
!            write(6,*) "lnorm2 =",lnorm2
!            write(6,*) ""

! Normalizing and updating Lanczos vectors
            q_j = (1.0/lnorm1)*vq_j
            p_j = (1.0/conjg(lnorm2))*vp_j


            C(:,g) = p_j
            D(:,g) = q_j
        enddo ! loop over Lanczos columns for two-sided CGS

! Checking the level of orthog with Frobenius norm after two-sided GS
        call zcheck_orthog_frob(C,D,m,n,frob)
        write(6,*) " " 
        write(6,*) "After Gram Schmidt with reorthog"
        write(6,*) " "
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "# column vecs = ",n,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

        return
        end subroutine zgs_mod_full_two_side_naive_x


        subroutine zgs_blas_full_two_side_x(A,B,m,n,C,D)
  ! Two sided full, classical Gram-Schmidt for complex matrices.
  ! This avoids the construction of the intermediate matrix needed for carrying out the G.S step.
  ! Instead it uses the BLAS complex dot products to get the final vector from Q*P^H*v

        implicit none
        integer :: g,i,j,k,l,m,n,f
        real(kind=8) :: temp1,frob,lterm,term
        complex(kind=8) :: zterm, lnorm1,lnorm2
        complex(kind=8) :: A(m,n),B(m,n),C(m,n),D(m,n)
        complex(kind=8) :: a_j(m),b_j(m),q_j(m),p_j(m),qq_j(m),pp_j(m)
        complex(kind=8),allocatable ::iden1(:,:),iden2(:,:),&
                                      Pj(:,:),Qj(:,:),atmp(:)
  ! Lapack routines
        complex( kind = kind( 1.0d0 ) ):: ZDOTC
        
      
        C = 0.0
        D = 0.0
        a_j = 0.0
        b_j = 0.0
        q_j = 0.0
        p_j = 0.0
        qq_j = 0.0
        pp_j = 0.0
         
        write(6,*) " "
        write(6,*) "Shape of A/B is",shape(A)
        write(6,*) " "
        write(6,*) " Two sided full Gram-Schmidt, one step"

! Checking the level of orthog with Frobenius norm before two-sided GS
        call zcheck_orthog_frob(A,B,m,n,frob)
        write(6,*) " "
        write(6,*) "!!!!!!!!!!!!!!!!!!!! " 
        write(6,*) "Before Gram Schmidt "
        write(6,*) " "
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "# column vecs = ",n,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

!------- Starting the two sided Gram-Schmidt  -------
        do g=1,n
!           write(6,*) " "
!           write(6,*) "********* Gram-Schmidt iteration=",g,"*********"
!           write(6,*) " "
           a_j = A(:,g)  
           b_j = B(:,g)

 
!!  Set the identity matrix
!           iden1 = 0.0
!           iden2 = 0.0
!           do f = 1,m
!             iden1(f,f) = 1.0d0
!             iden2(f,f) = 1.0d0
!           enddo  
           if (g .ne. 1) then

              allocate(Pj(m,g-1),Qj(m,g-1),atmp(g-1))
!              write(6,*) " "
!              write(6,*) "g = ",g
!              write(6,*) " Done allocating Pj and Qj"
              Pj = 0.0
              Qj = 0.0
              atmp = 0.0
              Pj = C(:,1:g-1) !originally A, B
              Qj = D(:,1:g-1)

! Using subroutine AB_dagger_v(A,B,m,n,v,d)
!              call AB_dagger_v(Pj,Qj,m,g-1,a_j,pp_j)  ! Doing matmul(Pj,conjg(transpose(Qj)))*a_j
!              call AB_dagger_v(Qj,Pj,m,g-1,b_j,qq_j)  ! Doing matmul(Pj,conjg(transpose(Qj)))*a_j

! Alternative formulation using associativity. Pj*Qj^dag*a_j = Pj*(Qj^dag*a_j)   
              atmp = matmul(conjg(transpose(Qj)),a_j)
              pp_j = matmul(Pj,atmp)

              p_j = a_j - pp_j
              q_j = b_j - qq_j

!! For Hermitian case, works well  
!              iden1 = iden1-matmul(Qj,conjg(transpose(Pj)))
!              iden2 = iden2-matmul(Pj,conjg(transpose(Qj)))
!
!! Gram Schmidt step
!              q_j = matmul(iden1,b_j)
!              p_j = matmul(iden2,a_j)
              deallocate(Pj,Qj)
            else
!              write(6,*) " "
!              write(6,*) " No need for re-biorthog "
              q_j = b_j
              p_j = a_j
            endif
 
!!  New norm
!            zterm = dot_product(q_j,p_j) ! following netlib Alg 7.13
            zterm = ZDOTC(m,q_j,1,p_j,1) ! following netlib Alg 7.13
            lnorm1 = sqrt(abs(zterm))  ! following yambo paper
            lnorm2 = conjg(zterm)/lnorm1
!            write(6,*) "lnorm2 =",lnorm2
!            write(6,*) ""
!            write(6,*) "g =",g,"zterm = ",zterm
           

! Normalizing and updating Lanczos vectors.
! Should keep the last vector added un normalized, since this will need to be normalized
! to calculate the coefficients. Or does it not matter if it is also normalized?
            q_j = (1.0/lnorm1)*q_j
            p_j = (1.0/conjg(lnorm2))*p_j


            C(:,g) = p_j
            D(:,g) = q_j
        enddo ! loop over Lanczos columns for two-sided CGS
      

! Checking the level of orthog with Frobenius norm after two-sided GS
        call zcheck_orthog_frob(C,D,m,n,frob)
        write(6,*) " " 
        write(6,*) "After Gram Schmidt "
        write(6,*) " "
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "# column vecs = ",n,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

        return
        end subroutine zgs_blas_full_two_side_x


        subroutine zgs_blas_sing_two_side_x(A,B,n,iter,eps1,p,q,pp,qq,max_orth)
!  Does a 2-sided, single step Gram Schmidt for vectors p and q on basis
! A and B, and then returns pp and qq vectors, not normalized.
! It does it without explicit matrix multiplication.

        implicit none
        integer :: i,j,k,m,f,n,iter,iterm,jterm
        complex(kind=kind(1.0d0)) :: A(n,iter),B(n,iter),p(n),q(n),pp(n),qq(n) 
        real( kind = kind( 1.0d0 ) ) :: eps,eps1,term_re,term_im,frob,fterm,tol,max_orth
        integer :: orthtot,poterm,qoterm,chk,oterm
        complex( kind = kind( 1.0d0 ) ) :: zterm,pterm,qterm,term1,lnorm1,lnorm2,vp(n),vq(n),wp(n),wq(n),&
                                           vp_j(n),vq_j(n),p_k(n),q_k(n),pp_j(n),qq_j(n)
        integer,allocatable :: porthopt(:),qorthopt(:),porthindx(:),qorthindx(:)
        complex( kind = kind( 1.0d0 ) ),allocatable :: iden(:,:),Pj(:,:),Qj(:,:),iden1(:,:),&
                                                 iden2(:,:),pw_jk(:),qw_jk(:),ptmp(:),qtmp(:)
  ! Lapack routines
        complex( kind = kind( 1.0d0 ) ):: ZDOTC
            write(6,*)" "
            write(6,*) "iter =",iter
            write(6,*) " "
            pp = 0.0
            qq = 0.0
            write(6,*) "shape of Planc =",shape(A),"shape Qlanc =",shape(B)
            write(6,*)" "
            

           if (iter .ne. 1) then

              allocate(Pj(n,iter),Qj(n,iter),ptmp(iter),qtmp(iter))
!              write(6,*) " "
!              write(6,*) "g = ",g
!              write(6,*) " Done allocating Pj and Qj"
              Pj = 0.0
              Qj = 0.0
              ptmp = 0.0
              qtmp = 0.0
              Pj = A(:,1:iter) !originally A, B
              Qj = B(:,1:iter)

!! Using subroutine AB_dagger_v(A,B,m,n,v,d)
!              call AB_dagger_v(Pj,Qj,n,iter,p,pp_j)  ! Doing matmul(Pj,conjg(transpose(Qj)))*a_j
!              call AB_dagger_v(Qj,Pj,n,iter,q,qq_j)  ! Doing matmul(Pj,conjg(transpose(Qj)))*a_j

! Alternative formulation using associativity. Pj*Qj^dag*a_j = Pj*(Qj^dag*a_j)   
!              ptmp = matmul(conjg(transpose(Qj)),p)
!              qtmp = matmul(conjg(transpose(Pj)),q)
!              pp_j = matmul(Pj,ptmp)
!              qq_j = matmul(Qj,qtmp)

! Doing it all one shot              
              pp_j = matmul(Pj,matmul(conjg(transpose(Qj)),p) )
              qq_j = matmul(Qj,matmul(conjg(transpose(Pj)),q) )
! Updating pp and qq
              pp = p - pp_j
              qq = q - qq_j

              deallocate(Pj,Qj,ptmp,qtmp)
            else
!              write(6,*) " "
!              write(6,*) " No need for re-biorthog "
              qq = q
              pp = p
            endif
         return
         end subroutine zgs_blas_sing_two_side_x


        subroutine zgs_blas_sing_two_side_minprbo_x(A,B,n,iter,eps1,p,q,pp,qq,oterm,max_orth)
!  Does a 2-sided, single step Gram Schmidt for vectors p and q on basis
! A and B, and then returns pp and qq vectors, not normalized.This does minimal partial rebiorthogonalization.
! It checks for the level of biorthogonality first before deciding to orthogonalize against a given vector.
! Uses associativity to do the matrix-matrix-vector product.
        implicit none
        integer :: i,j,k,m,f,n,iter,iterm,jterm
        complex(kind=kind(1.0d0)) :: A(n,iter),B(n,iter),p(n),q(n),pp(n),qq(n),pp_j(n),qq_j(n) 
        real( kind = kind( 1.0d0 ) ) :: eps,eps1,term_re,term_im,frob,fterm,tol,max_orth
        integer :: orthtot,poterm,qoterm,chk,oterm
        complex( kind = kind( 1.0d0 ) ) :: zterm,term1,lnorm1,lnorm2,vp(n),vq(n),wp(n),wq(n)
        integer,allocatable :: porthopt(:),qorthopt(:),porthindx(:),qorthindx(:)
        complex( kind = kind( 1.0d0 ) ),allocatable :: iden(:,:),Pj(:,:),Qj(:,:),iden1(:,:),&
                                                 iden2(:,:),pw_jk(:),qw_jk(:),ptmp(:),qtmp(:)
  ! Lapack routines
        complex( kind = kind( 1.0d0 ) ):: ZDOTC
!            write(6,*)" "
!            write(6,*) "iter =",iter
!            write(6,*) " "
            pp = 0.0
            qq = 0.0
            allocate(Pj(n,iter),Qj(n,iter),pw_jk(iter),qw_jk(iter))
            Pj = 0.0
            Qj = 0.0
            pw_jk = 0.0
            qw_jk = 0.0
            Pj = A(:,1:iter)
            Qj = B(:,1:iter)
!            write(6,*) "shape of Pj/Qj =",shape(Qj),"shape Planc =",shape(A)
!            write(6,*)" "
            
        ! Calculating local orthogonality, first normalizing p and q
            vp = 0.0
            vq = 0.0    
            term1 = dot_product(q,p)
!            term1 = wvfdot(backpsi_p,psi_p,allow,n)
            lnorm1 = sqrt(abs(term1))
            lnorm2 = conjg(term1)/lnorm1
            vp = (1.0/conjg(lnorm2))*p
            vq = (1.0/lnorm1)*q
            pw_jk = matmul(conjg(transpose(Qj)),vp)
            qw_jk = matmul(conjg(transpose(Pj)),vq)

            ! estimate of level of biorthogonality
            max_orth = max(maxval(abs(pw_jk)),maxval(abs(qw_jk)))

!            write(6,*) " "
!            write(6,*)"orthog,<p_j+1,q_k>, max",maxval(abs(pw_jk))
!            write(6,*) " "
!            write(6,*)pw_jk
!            write(6,*) " "
!            write(6,*)"orthog,<q_j+1,p_k>, max",maxval(abs(qw_jk))
!            write(6,*) " "
!            write(6,*) qw_jk
!            write(6,*)" "

            deallocate(Pj,Qj)

         ! Check for partial orthogonality
           allocate(porthopt(iter),qorthopt(iter))
           porthopt = 0
           qorthopt = 0
           poterm = 0
           qoterm = 0
         
         ! P's
           do i=1,iter
             zterm = pw_jk(i)
             term_re = abs(real(zterm))
             term_im = abs(aimag(zterm))
             if ((term_re .ge. eps1).or.(term_im .ge. eps1)) then
!               write(6,*) " "
!               write(6,*) "vector",i,"lost orth, pw_jk=",zterm,eps1
!               write(6,*) " "
               porthopt(i) = 1
             endif
           enddo
           poterm = sum(porthopt)
         ! Q's
           do i=1,iter
             zterm = qw_jk(i)
             term_re = abs(real(zterm))
             term_im = abs(aimag(zterm))
             if ((term_re .ge. eps1).or.(term_im .ge. eps1)) then
!               write(6,*) " "
!               write(6,*) "vector",i,"lost orth, qw_jk=",zterm,eps1
!               write(6,*) " "
               qorthopt(i) = 1
             endif
           enddo
           qoterm = sum(qorthopt)
           write(6,*)"---# vecs that lost biortho",poterm,qoterm,"---"
           write(6,*) " "

        ! Get the index of these vectors
           if (poterm .ne. qoterm) then
              oterm = max(poterm,qoterm)

           ! Make the indices the same if the # or orthog is diff
              if(poterm .ge. qoterm) then
                 qorthopt=qorthopt+(porthopt-qorthopt)
              else
                 porthopt=porthopt+(qorthopt-porthopt)
              endif
              poterm=oterm
              qoterm=oterm
           endif 
           
           allocate(porthindx(poterm),qorthindx(qoterm))
           porthindx = 0
           qorthindx = 0

    ! P's
           iterm = 1
           do i=1,iter
             if(porthopt(i) .eq. 1) then
               porthindx(iterm) = i
               iterm = iterm + 1
             endif
           enddo  

    ! Q's
           iterm = 1
           do i=1,iter
             if(qorthopt(i) .eq. 1) then
               qorthindx(iterm) = i
               iterm = iterm + 1
             endif
           enddo  
               
         deallocate(porthopt,qorthopt)  
!
    ! Construct the Pj and Qj
         if ((poterm .ge. 1) .or. (qoterm .ge. 1))then
!           write(6,*) "doing partial 2-sided GS as biorthog failed"
!           write(6,*) " poterm =",poterm,"qoterm =",qoterm
!           write(6,*) " "
           allocate(Pj(n,poterm),Qj(n,qoterm),ptmp(poterm),qtmp(poterm))
!           write(6,*) "allocated Pj, Qj"
!           write(6,*) "allocated iden1,iden2"
           Pj = 0.0
           Qj = 0.0
           ptmp = 0.0
           qtmp = 0.0
           do i=1,poterm
             iterm = porthindx(i)
             jterm = qorthindx(i)
!             write(6,*) "iterm =",iterm,"jterm =",jterm
!             write(6,*) " "
             Pj(:,i) = A(:,iterm)
             Qj(:,i) = B(:,jterm)
           enddo

!! Using subroutine AB_dagger_v(A,B,m,n,v,d)
!           call AB_dagger_v(Pj,Qj,n,poterm,p,pp_j)  ! Doing matmul(Pj,conjg(transpose(Qj)))*a_j
!           call AB_dagger_v(Qj,Pj,n,qoterm,q,qq_j)  ! Doing matmul(Pj,conjg(transpose(Qj)))*a_j

!! Alternative formulation using associativity. Pj*Qj^dag*a_j = Pj*(Qj^dag*a_j)   
           ptmp = matmul(conjg(transpose(Qj)),p)
           qtmp = matmul(conjg(transpose(Pj)),q)
           pp_j = matmul(Pj,ptmp)
           qq_j = matmul(Qj,qtmp)

! Update pp and qq              
           pp = p - pp_j
           qq = q - qq_j
         !  write(6,*) " "
         !  write(6,*) "pp =",pp
         !  write(6,*) " "
         !  write(6,*) "----- Comparing the level of orthog after rebiortho ----"
         !  write(6,*) " "

         !   pw_jk = matmul(conjg(transpose(A(:,1:iter))),pp)
         !   qw_jk = matmul(conjg(transpose(B(:,1:iter))),qq)

         !   write(6,*) " "
         !   write(6,*)"orthog,<p_j+1,q_k>, max",maxval(abs(pw_jk))
         !   write(6,*) " "
         !   write(6,*)pw_jk
         !   write(6,*) " "
         !   write(6,*)"orthog,<q_j+1,p_k>, max",maxval(abs(qw_jk))
         !   write(6,*) " "
         !   write(6,*) qw_jk
         !   write(6,*)" "
           deallocate(Pj,Qj,ptmp,qtmp)  
         else
!           write(6,*) " "
!           write(6,*) " No orthog was needed "
           pp = p
           qq = q 
         endif
         deallocate(porthindx,qorthindx) 
         deallocate(pw_jk,qw_jk)
         
!         write(6,*) " "
!         write(6,*) "Done with re-biorthogonalization"
         return
         end subroutine zgs_blas_sing_two_side_minprbo_x


        subroutine zgs_blas_sing_two_side_prbo_x(A,B,n,iter,eps1,p,q,pp,qq,oterm,max_orth,F_j)
!  Does a 2-sided, single step Gram Schmidt for vectors p and q on basis
! A and B, and then returns pp and qq vectors, not normalized.This a full rebiorthogonalization when it is deemed that
! biorthogonality is lost..
! It checks for the level of biorthogonality first before deciding to orthogonalize against a given vector.
! Uses associativity to do the matrix-matrix-vector product.
        implicit none
        integer :: i,j,k,m,f,n,iter,iterm,jterm,F_j
        complex(kind=kind(1.0d0)) :: A(n,iter),B(n,iter),p(n),q(n),pp(n),qq(n),pp_j(n),qq_j(n) 
        real( kind = kind( 1.0d0 ) ) :: eps,eps1,term_re,term_im,frob,fterm,tol,max_orth
        integer :: orthtot,poterm,qoterm,chk,oterm
        complex( kind = kind( 1.0d0 ) ) :: zterm,term1,lnorm1,lnorm2,vp(n),vq(n),wp(n),wq(n)
        integer,allocatable :: porthopt(:),qorthopt(:),porthindx(:),qorthindx(:)
        complex( kind = kind( 1.0d0 ) ),allocatable :: iden(:,:),Pj(:,:),Qj(:,:),iden1(:,:),&
                                                 iden2(:,:),pw_jk(:),qw_jk(:),ptmp(:),qtmp(:)
  ! Lapack routines
        complex( kind = kind( 1.0d0 ) ):: ZDOTC
!            write(6,*)" "
!            write(6,*) "iter =",iter
!            write(6,*) " "
            pp = 0.0
            qq = 0.0
            allocate(Pj(n,iter),Qj(n,iter),pw_jk(iter),qw_jk(iter))
            Pj = 0.0
            Qj = 0.0
            pw_jk = 0.0
            qw_jk = 0.0
            Pj = A(:,1:iter)
            Qj = B(:,1:iter)
!            write(6,*) "shape of Pj/Qj =",shape(Qj),"shape Planc =",shape(A)
!            write(6,*)" "
            
        ! Calculating local orthogonality, first normalizing p and q
            vp = 0.0
            vq = 0.0    
            term1 = dot_product(q,p)
!            term1 = wvfdot(backpsi_p,psi_p,allow,n)
            lnorm1 = sqrt(abs(term1))
            lnorm2 = conjg(term1)/lnorm1
            vp = (1.0/conjg(lnorm2))*p
            vq = (1.0/lnorm1)*q
            pw_jk = matmul(conjg(transpose(Qj)),vp)
            qw_jk = matmul(conjg(transpose(Pj)),vq)

            ! estimate of level of biorthogonality
            max_orth = max(maxval(abs(pw_jk)),maxval(abs(qw_jk)))

!            write(6,*) " "
!            write(6,*)"orthog,<p_j+1,q_k>, max",maxval(abs(pw_jk))
!            write(6,*) " "
!            write(6,*)pw_jk
!            write(6,*) " "
!            write(6,*)"orthog,<q_j+1,p_k>, max",maxval(abs(qw_jk))
!            write(6,*) " "
!            write(6,*) qw_jk
!            write(6,*)" "

            deallocate(Pj,Qj)

         ! Check for partial orthogonality
           allocate(porthopt(iter),qorthopt(iter))
           porthopt = 0
           qorthopt = 0
           poterm = 0
           qoterm = 0
         
         ! P's
           do i=1,iter
             zterm = pw_jk(i)
             term_re = abs(real(zterm))
             term_im = abs(aimag(zterm))
             if ((term_re .ge. eps1).or.(term_im .ge. eps1)) then
!               write(6,*) " "
!               write(6,*) "vector",i,"lost orth, pw_jk=",zterm,eps1
!               write(6,*) " "
               porthopt(i) = 1
             endif
           enddo
           poterm = sum(porthopt)
         ! Q's
           do i=1,iter
             zterm = qw_jk(i)
             term_re = abs(real(zterm))
             term_im = abs(aimag(zterm))
             if ((term_re .ge. eps1).or.(term_im .ge. eps1)) then
!               write(6,*) " "
!               write(6,*) "vector",i,"lost orth, qw_jk=",zterm,eps1
!               write(6,*) " "
               qorthopt(i) = 1
             endif
           enddo
           qoterm = sum(qorthopt)
           write(6,*)"---# vecs that lost biortho",poterm,qoterm,"---"
           write(6,*) " "

        ! Set F_j for determining the second rebiorthogonalization in iter+1
           oterm = max(poterm,qoterm)
           if (oterm .ge. 1) then
                F_j = iter+1
                print *
                print *, " F_j = ",F_j
                print *
           endif 

        ! Get the index of these vectors
           if (poterm .ne. qoterm) then
              oterm = max(poterm,qoterm)

           ! Make the indices the same if the # or orthog is diff
              if(poterm .ge. qoterm) then
                 qorthopt=qorthopt+(porthopt-qorthopt)
              else
                 porthopt=porthopt+(qorthopt-porthopt)
              endif
              poterm=oterm
              qoterm=oterm
           endif 
           
           allocate(porthindx(poterm),qorthindx(qoterm))
           porthindx = 0
           qorthindx = 0

    ! P's
           iterm = 1
           do i=1,iter
             if(porthopt(i) .eq. 1) then
               porthindx(iterm) = i
               iterm = iterm + 1
             endif
           enddo  

    ! Q's
           iterm = 1
           do i=1,iter
             if(qorthopt(i) .eq. 1) then
               qorthindx(iterm) = i
               iterm = iterm + 1
             endif
           enddo  
               
         deallocate(porthopt,qorthopt) 

!
    ! Construct the Pj and Qj
         if ((poterm .ge. 1) .or. (qoterm .ge. 1))then
!           write(6,*) "doing partial 2-sided GS as biorthog failed"
!           write(6,*) " poterm =",poterm,"qoterm =",qoterm
!           write(6,*) " "
           allocate(Pj(n,iter),Qj(n,iter),ptmp(iter),qtmp(iter))
!           write(6,*) "allocated Pj, Qj"
!           write(6,*) "allocated iden1,iden2"
           Pj = 0.0
           Qj = 0.0
           ptmp = 0.0
           qtmp = 0.0
           do i=1,iter
!             iterm = porthindx(i)
!             jterm = qorthindx(i)
!             write(6,*) "iterm =",iterm,"jterm =",jterm
!             write(6,*) " "
             Pj(:,i) = A(:,i)
             Qj(:,i) = B(:,i)
           enddo

!! Using subroutine AB_dagger_v(A,B,m,n,v,d)
!           call AB_dagger_v(Pj,Qj,n,poterm,p,pp_j)  ! Doing matmul(Pj,conjg(transpose(Qj)))*a_j
!           call AB_dagger_v(Qj,Pj,n,qoterm,q,qq_j)  ! Doing matmul(Pj,conjg(transpose(Qj)))*a_j

!! Alternative formulation using associativity. Pj*Qj^dag*a_j = Pj*(Qj^dag*a_j)   
           ptmp = matmul(conjg(transpose(Qj)),p)
           qtmp = matmul(conjg(transpose(Pj)),q)
           pp_j = matmul(Pj,ptmp)
           qq_j = matmul(Qj,qtmp)

! Update pp and qq              
           pp = p - pp_j
           qq = q - qq_j
         !  write(6,*) " "
         !  write(6,*) "pp =",pp
         !  write(6,*) " "
         !  write(6,*) "----- Comparing the level of orthog after rebiortho ----"
         !  write(6,*) " "

         !   pw_jk = matmul(conjg(transpose(A(:,1:iter))),pp)
         !   qw_jk = matmul(conjg(transpose(B(:,1:iter))),qq)

         !   write(6,*) " "
         !   write(6,*)"orthog,<p_j+1,q_k>, max",maxval(abs(pw_jk))
         !   write(6,*) " "
         !   write(6,*)pw_jk
         !   write(6,*) " "
         !   write(6,*)"orthog,<q_j+1,p_k>, max",maxval(abs(qw_jk))
         !   write(6,*) " "
         !   write(6,*) qw_jk
         !   write(6,*)" "
           deallocate(Pj,Qj,ptmp,qtmp)  
         else
!           write(6,*) " "
!           write(6,*) " No orthog was needed "
           pp = p
           qq = q 
         endif
         deallocate(porthindx,qorthindx) 
         deallocate(pw_jk,qw_jk)
         
!         write(6,*) " "
!         write(6,*) "Done with re-biorthogonalization"
         return
         end subroutine zgs_blas_sing_two_side_prbo_x

        subroutine AB_dagger_v(A,B,m,n,v,d)
  ! This subroutine takes two matrices A(m,n), B(m,n), and a vector v(m,1),
  ! and then carries out the matrix*(matrix^H)*vector operation
  ! d = A*(B^H)*v, without saving the intermediate matrix C(m,m). This avoids the memory
  ! limitations that come with solid state problems. It will use the Lapack routine for complex
  ! dot product ZDOTC

        implicit none
        integer :: g,i,j,k,l,m,n,f
        real(kind=8) :: temp1,frob,lterm,term
        complex(kind=8) :: zterm,zterm1,zterm2
        complex(kind=8) :: A(m,n),B(m,n)
        complex(kind=8) :: v(m),d(m),c_ith_row(m)
        complex(kind=8),allocatable ::iden1(:,:)
        
  ! Lapack routines
        complex( kind = kind( 1.0d0 ) ):: ZDOTC
        
        d = 0.0
        do i=1,m

        ! Get ith row of C = A*B^H. C has dim (m,m), so m columns
          c_ith_row = 0.0
          do j=1,m
            c_ith_row(j) = ZDOTC(n,A(i,:),1,B(j,:),1)
          enddo

        ! Get ith element of d vector  
          d(i) = ZDOTC(m,c_ith_row,1,v,1)
        end do ! loop for constructing vector d(m) element by element
        return
        end subroutine AB_dagger_v
