program planczos_minprbo
        ! This program carries out the non hermitian lanczos algorithm on a given matrix
        ! hermitian or non hermitian, using the matrix * matrix * vector operation in parallel, 
        ! following the template from the parallel_lu_solve.  
        ! Only processor (0,0) will carry out the Lanczos process for now, and then
        ! when there's a need to orthogonalize, will do this mmv product in parallel, 
        ! and collect the new vector from the parallel mmv at processor (0,0).
        ! A minimal partial rebiorthogonalization is carried out in this program based on some
        ! set tolerance value eps1 = f(eps), some fraction of the machine precision eps.
        !

        implicit none
        integer, parameter :: wp = kind(1.d0)
        complex(kind=wp), allocatable, dimension(:,:) :: A
        complex(kind=wp), allocatable, dimension(:,:) :: b
        complex(kind=wp), allocatable, dimension(:,:) :: y
        complex(kind=wp),allocatable :: A_g(:,:)
        complex(kind=wp) :: alpha, beta
        integer, dimension(9) :: desc_A_g,desc_A, desc_b, desc_y
        integer :: sys_order_m,sys_order_n,block_size_m,block_size_k,&
                   block_size_n,block_size_thresh,block_size_opt,nproc,i,j
        integer :: nprow, npcol, icontxt, myrow, mycol,mypnum,numroc,indxl2g
        integer :: mloc_A, nloc_A, mloc_b, nloc_b,lld_g,lld,rsrc,csrc,&
                   info,iA,jA
        
        ! Lanczos variables, parameters
        integer, parameter :: m=1000,mlanc =50
        integer :: k,iterm,jterm
        real(kind=wp) :: b0,temp1,tol,eps1,term_re,term_im,md,&
                         mc,frob,fnorm(mlanc),rand_num(m)
        real(kind=wp) :: u0,DZNRM2,ZDOTC,lnorm,DLANGE ! fnorm, frobenius norm
        real,parameter :: pi=3.141592654_wp, eps=1.11E-16
        integer,parameter :: sedsz=9
        integer           :: seed(sedsz),herm,mi,mj,chknorm(mlanc),orthtot,&
                             poterm,qoterm,oterm
        complex(kind=wp),parameter :: ii=cmplx(0.0_wp,1.0_wp)
        !  Non allocatable        
        complex(kind=wp) ::zalpha(mlanc),zbeta(mlanc),zdelta(mlanc), &
                &           zgamma(mlanc),zborth(mlanc) 
        complex(kind=wp) :: hterm,term,term1,lnorm1,lnorm2,H(m,m)
        complex(kind=wp):: pm1(m),qm1(m),p(m),q(m),s(m),tq(m),tp(m), &
                      & pp(m),qp(m),vp(m),vq(m) ! u is now vp/ v^_(j+1) ,wp is

        !  Matrices for calculating biorthogonality
        integer :: orth_ord(1,1)
        integer, allocatable :: porthopt(:),qorthopt(:),porthindx(:),&
                                qorthindx(:)
        complex(kind=wp) :: Planc(m,mlanc), Qlanc(m,mlanc)
        complex(kind=wp), allocatable :: Pj(:,:),Qj(:,:),pw_jk(:),qw_jk(:)
        logical :: biortho
                                       
        ! ScaLAPACK variables
        integer, dimension(9) :: desc_Pj,desc_Qj,desc_Pj_loc,desc_Qj_loc,desc_T_Pj_loc,desc_T_Qj_loc,&
                                 desc_pp_g,desc_qp_g,desc_pp_orth,desc_qp_orth,desc_pp_loc,desc_qp_loc,&
                                 desc_pp_orth_loc,desc_qp_orth_loc

        integer :: mloc_Pj,nloc_Pj,mloc_pp,nloc_pp,i_loc,j_loc
        complex(kind=wp),allocatable :: Pj_loc(:,:),Qj_loc(:,:),T_Pj_loc(:,:),T_Qj_loc(:,:),&
                                        pp_loc(:,:),qp_loc(:,:),pp_orth_loc(:,:),qp_orth_loc(:,:),&
                                        pp_orth(:,:),qp_orth(:,:),pp_g(:,:),qp_g(:,:)

        ! Read inputs : get number of processors and block size

        call blacs_pinfo(mypnum,nproc)
        open(10,file='block.in',status='old')
        rewind(10)
        read(10,*) block_size_thresh,block_size_opt,block_size_m
        close(10)

        ! Initialize BLACS grid and array descriptor vectors
        ! Set grid topology
        nprow=2; npcol=nproc/nprow
        print *,nproc,block_size_m,npcol

        call blacs_get(0,0,icontxt) !blacs_get(icontxt,what,val), when what=0/2, icontxt is ignored
        call blacs_gridinit(icontxt,'R',nprow,npcol)
        call blacs_gridinfo(icontxt,nprow,npcol,myrow,mycol)

        ! Generate Hamiltonian: basic example, Non Hermitian, diag = mod(i,6), off diag =
        biortho = .true.
        H=0.0
        zalpha=0.0
        zbeta=0.0
        zdelta=0.0
        zgamma=0.0

! elements of the matrix, m1=2, m2=2 for i and j respectively
        mi=6
        mj=5
        md = 1.0
        mc = 0.001
!        rand_num = 0.0
!        call random_seed(put=seed(1:sedsz))
!        call random_number(rand_num)

!  mod(i,6)+ii*(mod(j,5))
             do j=1,m
                do i=1,j
                   hterm=mod(i,mi)+ii*(mod(j,mj))
                   H(i,j)=hterm
!                   H(j,i)=conjg(hterm)
                   H(j,i)=1.0+hterm
!                   write(6,*) ""
!                   write(6,*) "---- i=",i,"j=",j,"H(i,j) =",H(i,j)
!                   write(6,*) ""
                enddo
             enddo
!! diagonals
!           do i=1,m
!              H(i,i)= md*mod(i,mi)
!           enddo

!! diagonals
        do i=1,m
!              H(i,i)= md*mod(i,mi)+ii*(mod(i,mi))
!              H(i,i)= md*i + ((-1)**i)*ii*rand_diag(i)
            H(i,i)= md*i !+ ((-1)**i)*ii*rand_diag(i)
! Example from : "The Lanczos Alg with partial reorthogonalization"!              
!               H(i,i) = 1.0d0*i**2
!               write(6,*) H(i,i)
        enddo   

        ! Initial vectors for the Lanczos iteration
        iA = 1
        jA = 1
        alpha = 1.0_wp
        beta = 0.0_wp

        if (myrow .eq. 0 .and. mycol .eq. 0) then
          
          ! Set tolerance
          eps1 = sqrt(eps)

          ! Initialize other arrays
          zalpha=0.0_wp
          zbeta=0.0_wp
          zgamma=0.0_wp
          zborth = 0.0_wp
          fnorm = 0.0_wp
          pm1=0.0_wp
          qm1=0.0_wp
          p=0.0_wp
          q=0.0_wp
          s=0.0_wp
          tp=0.0_wp
          tq=0.0_wp
          pp=0.0_wp
          qp=0.0_wp
          rand_num=0.0_wp

!  Generating the test vector , v_1
          call random_seed(put=seed(1:sedsz))
          call random_number(rand_num)
!        s(1)=1.0d0+1.0*ii
!        s(mhalf)=1.0d0+1.0*ii
          s=1.0d0+1.0*ii
!          s=1.0_wp*rand_num+0.5_wp*ii*rand_num
          u0=dot_product(s,s)
          u0=sqrt(u0)
          p=s/u0
          q=p

!        if (myrow .eq. 0 .and. mycol .eq. 0) then
          Planc = 0.0_wp
          Qlanc = 0.0_wp
          Planc(:,1) = p
          Qlanc(:,1) = q

          ! Start the counting of number of orthogs
          orthtot = 0
        endif ! Only processor (0,0) handles arrays and does computations

        ! Start of Lanczos loop
        do k=1,mlanc
                if (myrow .eq. 0 .and. mycol .eq. 0) then
                  print *
                  print *,"******** m, k = ",m,k,"***********"
                  print *
                endif

!                if (k==1) then
!                     if (myrow .eq. 0 .and. mycol .eq. 0) then
!                       zbeta(k)=0.0_wp
!                       zgamma(k)=0.0_wp
!                     endif ! Only processor (0,0) handles arrays and does computations
!                end if 

                if (k .ne.1) then
                     if (myrow .eq. 0 .and. mycol .eq. 0) then
                       pm1 = Planc(:,k-1)
                       qm1 = Qlanc(:,k-1)
                     endif ! Only processor (0,0) handles arrays and does computations
                else
                     if (myrow .eq. 0 .and. mycol .eq. 0) then
                       pm1 = 0.0_wp
                       qm1 = 0.0_wp
                       zbeta(k)=0.0_wp
                       zgamma(k)=0.0_wp
                     endif ! Only processor (0,0) handles arrays and does computations

                end if 

        !  Checking biorthogonality with the frobenius norm, 
        ! and carry out Lanczos procedure

                if (myrow .eq. 0 .and. mycol .eq. 0) then
                    allocate(Pj(m,k),Qj(m,k))
                    Pj = 0.0_wp
                    Qj = 0.0_wp
                    Pj = Planc(:,1:k)
                    Qj = Qlanc(:,1:k)   
                    call zcheck_orthog_frob(Pj,Qj,m,k,frob)
                    fnorm(k) = frob
                    deallocate(Pj,Qj)

                    tq=matmul(H,q) !doing A*v_j
                    tp=matmul(transpose(conjg(H)),p) ! doing A^H*w_j
                    term=dot_product(p,tq) ! Following the SLEPc paper for Hermitian case
        !   Including this Paige-like regularization, from hermitian lanczos, for better convergence
                    qp=tq-zgamma(k)*qm1 ! v^_(j+1)
                    pp=tp-conjg(zbeta(k))*pm1 ! v^_(j+1)

                    term=dot_product(p,qp) ! Following the SLEPc paper for Hermitian case
                    zalpha(k) = term
                    qp=qp-zalpha(k)*q ! u_(j+1)  Paige-like regularization
                    pp=pp-conjg(zalpha(k))*p ! u_(j+1), Paige-like reg
                endif ! Only processor(0,0) does computations

                ! Doing rebiorthogonalization
                if (biortho) then

                        ! Get local level of orthog
                        if(myrow .eq. 0 .and. mycol .eq. 0)then
                            allocate(Pj(m,k),Qj(m,k),pw_jk(k),qw_jk(k))
                            Pj = 0.0
                            Qj = 0.0
                            pw_jk = 0.0
                            qw_jk = 0.0
                            Pj = Planc(:,1:k)
                            Qj = Qlanc(:,1:k)

                !!  Calculate local orthog
                            vp = 0.0
                            vq = 0.0
                            term1 = dot_product(qp,pp) ! following netlib Alg 7.13
                            lnorm1 = sqrt(abs(term1))  ! following yambo paper
                            lnorm2 = conjg(term1)/lnorm1
                            vp=(1.0/conjg(lnorm2))*pp ! w_(j+1)=w^_(j+1)/beta_(j+1)
                            vq=(1.0/lnorm1)*qp! v_(j+1)=v^_(j+1)/delta_(j+1)
                            pw_jk = matmul(conjg(transpose(Qj)),vp)
                            qw_jk = matmul(conjg(transpose(Pj)),vq)

                ! deallocate Pj, Qj
                            deallocate(Pj,Qj)

                !!  Check for partial orthog
                            allocate(porthopt(k),qorthopt(k))
                            porthopt = 0
                            qorthopt = 0
                            poterm = 0
                            qoterm = 0
                ! P's
                            do i=1,k
                              term = pw_jk(i)
                              term_re = abs(real(term))
                              term_im = abs(aimag(term))
                              if ((term_re .ge. eps1).or.(term_im .ge. eps1)) then
                                porthopt(i) = 1
                              endif
                            enddo
                            poterm = sum(porthopt)

                ! Q's
                            do i=1,k
                              term = qw_jk(i)
                              term_re = abs(real(term))
                              term_im = abs(aimag(term))
                              if ((term_re .ge. eps1).or. (term_im .ge. eps1)) then
                                qorthopt(i) = 1
                              endif
                            enddo
                            qoterm = sum(qorthopt)
                            write(6,*) "---# vecs that lost biortho",poterm,qoterm,"---"
                            write(6,*) " "
                            deallocate(pw_jk,qw_jk)

                !! Get the index of these vectors
                             if (poterm .ne. qoterm) then
                               oterm = max(poterm,qoterm)
                                
                !  Make the indices the same if the num of orthog is diff
                               if (poterm .ge. qoterm) then
                                 qorthopt = qorthopt + (porthopt-qorthopt)
                               else
                                 porthopt = porthopt + (qorthopt-porthopt)
                               endif
                               poterm = oterm
                               qoterm = oterm
                             endif

                             allocate(porthindx(poterm),qorthindx(qoterm))
                             porthindx = 0
                             qorthindx = 0
                            
                ! P's      
                            iterm = 1
                            do i=1,k
                              if(porthopt(i) .eq. 1) then
                                porthindx(iterm) = i
                                iterm = iterm + 1
                              endif
                            enddo
                ! Q's      
                            iterm = 1
                            do i=1,k
                              if(qorthopt(i) .eq. 1) then
                                qorthindx(iterm) = i
                                iterm = iterm + 1
                              endif
                            enddo
                            deallocate(porthopt,qorthopt)

                            !
                        endif ! Only processor(0,0) does computations

                       ! Broadcast the value of oterm
                       if(myrow .eq. 0 .and. mycol .eq. 0)then
                          orth_ord = 0
                          oterm = max(poterm,qoterm)
                          orth_ord(1,1) = oterm
                          print *,"oterm =",oterm,orth_ord(1,1)
                          print *,"poterm=",poterm,"qoterm=",qoterm
                          call igebs2d(icontxt,"All"," ",1,1,orth_ord,1)
                       else
                          call igebr2d(icontxt,"All"," ",1,1,orth_ord,1,0,0)
                       endif 
                       
                       oterm = orth_ord(1,1)
                       print *," # of vecs that lost orthog =",oterm
                       print *

                        if (oterm .ge. 1) then
                               ! Determine block size_k/n
                                if(oterm .le. block_size_thresh)then
                                   block_size_n = oterm
                                else
                                   block_size_n = block_size_opt
                                endif

                                sys_order_m = m
                                sys_order_n = orth_ord(1,1)


                                mloc_A = numroc(sys_order_m,block_size_m,myrow,0,nprow)
                                nloc_A = numroc(sys_order_n,block_size_n,mycol,0,npcol)
                                mloc_b = numroc(sys_order_m,block_size_m,myrow,0,nprow)
                                nloc_b = numroc(1,1,mycol,0,npcol)

                                ! Allocate global arrays
                                if(myrow .eq. 0 .and. mycol .eq. 0)then
                                  allocate(Pj(sys_order_m,sys_order_n),&
                                           Qj(sys_order_m,sys_order_n),&
                                           pp_g(sys_order_m,1),&
                                           qp_g(sys_order_m,1),&
                                           pp_orth(sys_order_m,1),&
                                           qp_orth(sys_order_m,1))

                                  Pj = 0.0_wp
                                  Qj = 0.0_wp
                                  do i=1,sys_order_n
                                    iterm = porthindx(i)
                                    jterm = qorthindx(i)
                                    Pj(:,i) = Planc(:,iterm)
                                    Qj(:,i) = Qlanc(:,jterm)
                                  enddo
                                  pp_g = 0.0_wp
                                  qp_g = 0.0_wp
                                  pp_orth = 0.0_wp
                                  qp_orth = 0.0_wp
                                  pp_g(:,1) = pp
                                  qp_g(:,1) = qp

                                  ! deallocate local orthog arrays
                                  deallocate(porthindx,qorthindx)
                                endif ! Only processor(0,0) does computations

                                ! Allocate space for local matrices
                                allocate(Pj_loc(mloc_A,nloc_A),&
                                         Qj_loc(mloc_A,nloc_A),&
                                         T_Qj_loc(mloc_A,mloc_A))


                                ! 1-D arrays
                                allocate(pp_loc(mloc_b,nloc_b),&
                                         qp_loc(mloc_b,nloc_b),&
                                         pp_orth_loc(mloc_b,nloc_b),&
                                         qp_orth_loc(mloc_b,nloc_b))

                                ! Initialize descriptors
                                lld_g = max(numroc(sys_order_m,sys_order_m,myrow,0,nprow),1)
                                lld = max(mloc_A,1)

                        !        if(myrow .eq. 0 .and. mycol .eq. 0)then
        !                          print *
        !                          print *
        !                          print *, "myrow = ",myrow," mycol =",mycol
        !                          print *, "mloc_A = ",mloc_A," nloc_A =",nloc_A
        !                          print *, "mloc_b = ",mloc_b," nloc_b =",nloc_b
        !                          print *, "block_size_m =",block_size_m,"block_size_n =",block_size_n
        !                          print *, " lld_g, lld"
        !                          print *, lld_g,lld
        !                          print *

                                ! Descriptors for global arrays

                                ! 2-D
                                !call descinit(desc_A_g,sys_order_m,sys_order_n,sys_order_m,sys_order_n,0,0,&
                                !              icontxt,lld_g,info)
                                call descinit(desc_Pj,sys_order_m,sys_order_n,sys_order_m,sys_order_n,0,0,&
                                              icontxt,lld_g,info)
                                call descinit(desc_Qj,sys_order_m,sys_order_n,sys_order_m,sys_order_n,0,0,&
                                              icontxt,lld_g,info)

                                ! 1-D
                                call descinit(desc_pp_g,sys_order_m,1,sys_order_m,1,0,0,&
                                              icontxt,lld_g,info) 
                                call descinit(desc_qp_g,sys_order_m,1,sys_order_m,1,0,0,&
                                              icontxt,lld_g,info) 

                                ! Local descriptors

                                ! 2-D
                                !call descinit(desc_A,sys_order_m,sys_order_n,block_size_m,block_size_n,0,0,&
                                !              icontxt,lld,info)
                                call descinit(desc_Pj_loc,sys_order_m,sys_order_n,&
                                           block_size_m,block_size_n,0,0,icontxt,lld,info)
                                call descinit(desc_Qj_loc,sys_order_m,sys_order_n,&
                                           block_size_m,block_size_n,0,0,icontxt,lld,info)
                                ! Temp matrix(m,m)
                                !call descinit(desc_T_Pj_loc,sys_order_m,sys_order_m,block_size_m,block_size_m,0,0,&
                                !              icontxt,lld,info)
                                call descinit(desc_T_Qj_loc,sys_order_m,sys_order_m,&
                                          block_size_m,block_size_m,0,0,icontxt,lld,info)

                                ! 1-D
                                call descinit(desc_pp_loc,sys_order_m,1,block_size_m,&
                                          1,0,0,icontxt,lld,info) 
                                call descinit(desc_qp_loc,sys_order_m,1,block_size_m,&
                                          1,0,0,icontxt,lld,info) 

                                ! orthogonalized vectors
                                call descinit(desc_pp_orth_loc,sys_order_m,1,block_size_m,&
                                          1,0,0,icontxt,lld,info) 
                                call descinit(desc_qp_orth_loc,sys_order_m,1,block_size_m,&
                                          1,0,0,icontxt,lld,info) 

                               ! Distribute the matrix Pj and Qj and have each processor print theirs locally.

                               ! 2-D arrays
                               ! call pzgeadd('n',sys_order_m,sys_order_n,alpha,A_g,iA,jA,desc_A_g,beta,A,iA,jA,desc_A)
                                call pzgeadd('n',sys_order_m,sys_order_n,alpha,Pj,iA,jA,&
                                         desc_Pj,beta,Pj_loc,iA,jA,desc_Pj_loc)
                                call pzgeadd('n',sys_order_m,sys_order_n,alpha,Qj,iA,jA,&
                                         desc_Qj,beta,Qj_loc,iA,jA,desc_Qj_loc)

                               ! 1-D arrays      
                                call pzgeadd('n',sys_order_m,1,alpha,pp_g,iA,jA,desc_pp_g,beta,&
                                             pp_loc,iA,jA,desc_pp_loc)
                                call pzgeadd('n',sys_order_m,1,alpha,qp_g,iA,jA,desc_qp_g,beta,&
                                             qp_loc,iA,jA,desc_qp_loc)

                     ! For pp_orth Pj*Qj^H
                     !subroutine mmdag_v(A,B,C,x,y,alpha,beta,sys_order_m,sys_order_n, &
                     !           desc_A,desc_B,desc_C,desc_x,desc_y)

                                call mmdag_v(Pj_loc,Qj_loc,T_Qj_loc,pp_loc,pp_orth_loc,alpha,beta,&
                                             sys_order_m,sys_order_n,desc_Pj_loc,desc_Qj_loc,&
                                             desc_T_Qj_loc,desc_pp_loc,desc_pp_orth_loc)
        !
                             ! For qp_orth Qj*Pj^H
                                T_Qj_loc = 0.0_wp
                                call mmdag_v(Qj_loc,Pj_loc,T_Qj_loc,qp_loc,qp_orth_loc,alpha,beta,&
                                             sys_order_m,sys_order_n,desc_Qj_loc,desc_Pj_loc,&
                                             desc_T_Qj_loc,desc_qp_loc,desc_qp_orth_loc)

                            ! Collect all the local parts of pp_orth and qq_orth
                            ! and then have them on processor (0,0)
                                call descinit(desc_pp_orth,sys_order_m,1,sys_order_m,1,0,0,icontxt,lld_g,info)
                                call descinit(desc_qp_orth,sys_order_m,1,sys_order_m,1,0,0,icontxt,lld_g,info)
                              
                                call pzgemr2d(sys_order_m,1,pp_orth_loc,1,1,desc_pp_orth_loc,pp_orth,&
                                           1,1,desc_pp_orth,icontxt)
                                call pzgemr2d(sys_order_m,1,qp_orth_loc,1,1,desc_qp_orth_loc,qp_orth,&
                                           1,1,desc_qp_orth,icontxt)


                                ! Deallocate space
                                if(myrow .eq. 0 .and. mycol .eq. 0)then
                                  
                                  ! Finish the orthogonalization process
                                  print *
                                  print *," Finished the re-orthogonalization process."
                                  print *
                                  pp = pp - pp_orth(:,1)
                                  qp = qp - qp_orth(:,1)
                                  orthtot = orthtot + 2*oterm
                                  deallocate(Pj,Qj,pp_g,qp_g,pp_orth,qp_orth)
                                endif
                                deallocate(Pj_loc,Qj_loc,T_Qj_loc,pp_loc,qp_loc,&
                                           pp_orth_loc,qp_orth_loc)
                        else  
                                if(myrow .eq. 0 .and. mycol .eq. 0)then
                                   deallocate(porthindx,qorthindx)
                                endif
                        endif ! End of oterm .ge. 1 condition
                endif ! End of re-biorthogonalization condition

                ! Update coeffs and vectors
                if(myrow .eq. 0 .and. mycol .eq. 0) then
                    term1 = dot_product(qp,pp) ! following netlib Alg 7.13

        !  A rough measure of the biorthogonality
                    zborth(k) = term1

                    lnorm1 = sqrt(abs(term1))  ! following yambo paper
        !          lnorm1 = sqrt(abs(term1*conjg(term1)))  ! following yambo paper
                    lnorm2 = conjg(term1)/lnorm1
                    if (k .ne. mlanc) then
                      zbeta(k+1)= lnorm1 ! following yambo paper
                      zgamma(k+1)=lnorm2 ! following Y. Saad's paper
                    endif

        !  Updating v_(j-1) for the next run (j+1) 
!                    qm1=q
!                    pm1=p

        !  Updating v_(j+1), to be used at start of (j+1) iteration, can use the tol here also
                    if (k .ne. mlanc) then
                      q=(1.0_wp/lnorm1)*qp! v_(j+1)=v^_(j+1)/delta_(j+1)
                      p=(1.0_wp/conjg(lnorm2))*pp ! w_(j+1)=w^_(j+1)/beta_(j+1)

        !  Add new vectors 
                      Planc(:,k+1) = p
                      Qlanc(:,k+1) = q
                    endif
                endif ! ----- Only processor (0,0)  Does the Lanczos iterations


        enddo ! End of Lanczos loop

!!  Saving alpha, beta to array
        if(myrow .eq. 0 .and. mycol .eq. 0) then
           open(10,file='par_lanc_nh_mat_pro_6_procs.dat')
           open(11,file='par_fnorm_nh_mat_pro_6_procs.dat')
           rewind(10)

           rewind(11)
           write(6,*) "----- Printing alpha, beta and gamma -------"
           print *
           print *," Total # of orthogs = ",orthtot
           print *
           do i=1,mlanc
!             write(6,*) ""
!             write(6,*) "lanc iteration =",i
!             write(6,*) ""
             write(10,'(1i3,6e16.8)') i,zalpha(i),zbeta(i),zgamma(i)
             write(11,'(1i3,2e16.8)') i,fnorm(i),abs(aimag(zborth(i)))
!             write(6,*) ""
!             write(6,*)'alpha=',zalpha(i)
!             write(6,*)""
!             write(6,*)'beta(i)=',zbeta(i)
!             write(6,*)""
!             write(6,*)'delta(i)=',zgamma(i)
!             write(6,*) ""
           enddo
           close(11)
           close(10)
        endif ! ----- Only processor (0,0)  Does the writing of files
        ! exit BLACS
        call blacs_exit(0)
        contains
        subroutine get_parameters(sys_order_m,sys_order_n,block_size_m,&
                                          block_size_n,nproc)
                integer, intent(out) :: sys_order_m,sys_order_n,block_size_m,&
                                                   block_size_n,nproc
                integer :: mypnum
                open(10,file='pmmv.in',status='old')
                rewind(10)
                read(10,*) sys_order_m, sys_order_n
                read(10,*) block_size_m, block_size_n
                close(10)
                call blacs_pinfo(mypnum,nproc)
        end subroutine get_parameters

        subroutine init_grid(sys_order_m,sys_order_n,block_size_m,block_size_n,nprow,npcol,&
                             myrow,mycol,icontxt,desc_A,desc_B,desc_C,desc_x,desc_y,&
                             mloc_A,nloc_A,mloc_B,nloc_B,mloc_C,nloc_C,mloc_x,nloc_x)

                 integer, intent(in) :: sys_order_m,sys_order_n,block_size_m,block_size_n
                 integer, intent(out) :: nprow,npcol,myrow,mycol,icontxt
                 integer, intent(out) :: mloc_A,nloc_A,mloc_B,nloc_B,mloc_C,nloc_C, &
                                         mloc_x,nloc_x
                 integer, dimension(9), intent(out) :: desc_A,desc_B,desc_C,desc_x,desc_y
                 integer :: numroc,info

                 ! Set processor grid topology. Square is best
!                 nprow=1; npcol=nproc
!                 nprow = 2; npcol = nproc/nprow
                 nprow=int(sqrt(real(nproc))); npcol=int(nproc/nprow)


                 call blacs_get(-1,0,icontxt) !blacs_get(icontxt,what,val), when what=0/2, icontxt is ignored
                 call blacs_gridinit(icontxt,'R',nprow,npcol)
                 call blacs_gridinfo(icontxt,nprow,npcol,myrow,mycol)

                 if(myrow .eq. 0 .and. mycol .eq. 0) then
                   print *
!                   print *, " system dim: (",sys_order_m,sys_order_n,")
                   print *, "--------- Setting up the system ---------"
                   print * ," icontxt =",icontxt
                   print *, " system dim: (",sys_order_m,sys_order_n,"  )"
                   print *, " nprow = ",nprow," npcol = ",npcol
                   print *, " block size: (",block_size_m,block_size_n,"  )"
                   print *
                 endif

                 ! Get dimensions of arrays locally owned
                 ! numroc(n,nb,iproc,isrcproc,nprocs)

                 mloc_A = numroc(sys_order_m,block_size_m,myrow,0,nprow)
                 nloc_A = numroc(sys_order_n,block_size_n,mycol,0,npcol)
                 mloc_B = mloc_A
                 mloc_C = mloc_A
                 nloc_B = nloc_A
                 nloc_C = nloc_A
                 mloc_x = numroc(sys_order_m,block_size_m,myrow,0,nprow)
                 nloc_x = numroc(1,1,mycol,0,npcol)
                 print *
                 print *, "myrow = ",myrow," mycol =",mycol
                 print *, "mloc_A = ",mloc_A," nloc_A =",nloc_A
                 print *, "mloc_B = ",mloc_B," nloc_B =",nloc_B
                 print *, "mloc_x = ",mloc_x," nloc_x =",nloc_x
!                 print *, "block_size =",block_size
                 print * 

                 ! Set up descriptor vectors for matrix a and vector b
                 ! call descinit(desc,m,n,mb,nb,irsrc,icsrc,ictxt,lld,info)
                 ! lld is determined by max(1,numroc(n,mb,iproc,isrcproc,nprocs)), for the rows

                 ! Matrices
                 call descinit(desc_A,sys_order_m,sys_order_n,block_size_m,block_size_n,0,0, &
                               icontxt,mloc_A,info)
                 call descinit(desc_B,sys_order_m,sys_order_n,block_size_m,block_size_n,0,0, &
                               icontxt,mloc_B,info)
                 call descinit(desc_C,sys_order_m,sys_order_n,block_size_m,block_size_n,0,0, &
                               icontxt,mloc_C,info)
                 ! Vectors
                 call descinit(desc_x,sys_order_m,1,block_size_m,1,0,0, &
                               icontxt,mloc_x,info)
                 call descinit(desc_y,sys_order_m,1,block_size_m,1,0,0, &
                               icontxt,mloc_x,info)

                 if (info .ne. 0)print*, 'descinit info: ',info
        end subroutine init_grid

        subroutine init_system(A_g,B_g,x_g,A,B,C,x,y,sys_order_m,sys_order_n, &
                               block_size_m,block_size_n, &
                               mloc_A,nloc_A,mloc_B,nloc_B, &
                               mloc_C,nloc_C,mloc_x,nloc_x, &
                               myrow,mycol,nprow,npcol)
                  implicit none

                  ! Real arrays
!                  real(kind=wp), dimension(:,:), intent(out) :: A
!                  real(kind=wp), dimension(:,:), intent(out) :: B
!                  real(kind=wp), dimension(:,:), intent(out) :: C
!                  real(kind=wp), dimension(:,:), intent(out) :: x
!                  real(kind=wp), dimension(:,:), intent(out) :: y

                  ! Complex arrays
                  complex(kind=wp), dimension(:,:), intent(out) :: A
                  complex(kind=wp), dimension(:,:), intent(out) :: B
                  complex(kind=wp), dimension(:,:), intent(out) :: C
                  complex(kind=wp), dimension(:,:), intent(out) :: x
                  complex(kind=wp), dimension(:,:), intent(out) :: y

                  integer, intent(in) :: sys_order_m,sys_order_n,mloc_A,nloc_A,mloc_B,nloc_B,&
                          mloc_C,nloc_C,mloc_x,nloc_x,block_size_m,block_size_n
                  integer, intent(in) :: myrow,mycol,nprow,npcol
                  integer :: indxl2g,i,i_loc,j,j_loc

                  ! Real arrays
!                  real(kind=wp), intent(in) :: A_g(sys_order_m,sys_order_n),&
!                                               B_g(sys_order_m,sys_order_n),&
!                                               x_g(sys_order_m,1)

                  ! Comlex arrays
                  complex(kind=wp), intent(in) :: A_g(sys_order_m,sys_order_n),&
                                               B_g(sys_order_m,sys_order_n),&
                                               x_g(sys_order_m,1)
                  ! Set up portion of system matrix owned by each processor
                  ! indxl2g(indxloc,nb,iproc,isrcproc,nprocs)

                  ! Initialize the parts of global matrices A owned locally
                  do i_loc=1,mloc_A
                     do j_loc=1,nloc_A
                        i = indxl2g(i_loc,block_size_m,myrow,0,nprow)
                        j = indxl2g(j_loc,block_size_n,mycol,0,npcol)
!                        A(i_loc,j_loc) = 1.0_wp 
!                        B(i_loc,j_loc) = 5.0_wp 
!                        C(i_loc,j_loc) = 0.0_wp 
                        A(i_loc,j_loc) = A_g(i,j)
                        B(i_loc,j_loc) = B_g(i,j)
                        C(i_loc,j_loc) = 0.0_wp
                      enddo
                   enddo

                   
!                   print *
!                   print *, "myrow = ",myrow," mycol =",mycol
!                   print *, "mloc_A = ",mloc_A," nloc_A =",nloc_A
!                   print *, "mloc_B = ",mloc_B," nloc_B =",nloc_B
!                   print *, "mloc_C = ",mloc_C," nloc_C =",nloc_C
!                   print * 

                   ! Initialize parts of global vector owned locally
                   do j_loc=1,nloc_x
                     j=indxl2g(j_loc,1,mycol,0,npcol)
                     do i_loc=1,mloc_x
                       i = indxl2g(i_loc,block_size_m,myrow,0,nprow)
!                       x(i_loc,j_loc) = 1.0_wp
                       x(i_loc,j_loc) = x_g(i,j)
                       y(i_loc,j_loc) = 0.0_wp
!                       x(i_loc) = 1.0_wp
!                       y(i_loc) = 0.0_wp
                      enddo
                   enddo
        end subroutine init_system

        subroutine mm_v(A,B,C,x,y,alpha,beta,sys_order_m,sys_order_n, &
                        desc_A,desc_B,desc_C,desc_x,desc_y)

                     ! Does matrix * vector using pdgemv or pzgemv

!                     ! Real arrays
!                     real(kind=wp), dimension(:,:), intent(in) :: A
!                     real(kind=wp), dimension(:,:), intent(in) :: B
!                     real(kind=wp), dimension(:,:), intent(inout) :: C
!                     real(kind=wp), dimension(:,:), intent(in) :: x
!                     real(kind=wp), dimension(:,:), intent(out) :: y
!                     real(kind=wp) :: alpha,beta

                     ! Complex arrays
                     complex(kind=wp), dimension(:,:), intent(in) :: A
                     complex(kind=wp), dimension(:,:), intent(in) :: B
                     complex(kind=wp), dimension(:,:), intent(inout) :: C
                     complex(kind=wp), dimension(:,:), intent(in) :: x
                     complex(kind=wp), dimension(:,:), intent(out) :: y
                     complex(kind=wp) :: alpha,beta

                     integer, intent(in) :: sys_order_m,sys_order_n
                     integer, dimension(9), intent(in) :: desc_A,desc_B,desc_C,desc_x,desc_y
                     integer :: iA,jA,iB,jB,iC,jC,ix,jx,nrhs
                     integer :: info
                     integer, dimension(sys_order_m+block_size_m) :: ipvt

                     ! Using the whole matrices and vectors
                     iA=1; jA=1; iB=1; jB=1; iC=1; jC=1; ix=1; jx=1
                     
!                     print *
!                     print *, " alpha = ",alpha," beta = ",beta
!                     print *, "sys_order_m =",sys_order_m,"sys_order_n =",sys_order_n
!                     print *

                     ! Holds up execution of all processes within the indicated scope
                     ! seems to be present ony with the INTEL BLACS
                     call blacs_barrier(icontxt, 'ALL')

                     ! pdgemm routine
                     ! Routine to do general matrix times matrix operation
                     ! C = alpha*A*B + beta*C,  A(m,k), B(k,n) --> C(m,n)
                     ! transa - form of matrix A used, 'N' for normal, 'T' for transpose
                     ! 'C' for conjugate transpose
                     ! transb - form of matrix B used
                     ! m - the num of rows of global matrix C, if transa='N' its num of rows in A
                     ! n - the num of cols of global matrix C, if transb='N' its num of cols in B
                     ! k- if transa='N',the num of cols of global matrix A, otherwise num of cols in A.
                     !     If transb='N' its num of rows in B, otherwise, num of cols in B
                     ! alpha - global scalar alpha in equation above.
                     ! a - local part of global matrix A, also identifies 1st elem of local array
                     ! iA - the row index of global matrix A, identifying the first row to use in calc
                     ! jA - the col index of global matrix A, indentifying the first col to use
                     ! desc_a - descriptor for matrix a
                     ! b - local part of global matrix B, also identifies 1st elem of local array
                     ! iC - the row index of global matrix B, identifying the first row to use in calc
                     ! jC - the col index of global matrix B, indentifying the first col to use
                     ! desc_b - descriptor for matrix b
                     ! beta - global scalar beta in equation above.
                     ! c - local part of global matrix B, also identifies 1st elem of local array
                     ! iC - the row index of global matrix C, identifying the first row to use in calc
                     ! jC - the col index of global matrix C, indentifying the first col to use
                     ! desc_c - descriptor for matrix c
                     !
                     ! ON RETURN
                     ! a - the updated local part of global matrix A, thats factorized
                     ! ipvt - the local part of the global vector ipvt, which contains the pivot indices
                     ! info - if 0, factorization is succesful, otherwise A is singular, global

                     ! call pdgemm(transa,transb,m,n,k,alpha,a,ia,ja,desc_a,b,ib,jb,desc_b,beta,
                     !              c,ic,jc,desc_c)
                     call pzgemm('n','n',sys_order_m,sys_order_n,sys_order_n,alpha,A,iA,jA,desc_A,B,iB,jB,desc_B,beta,&
                                 C,iC,jC,desc_C)
                     if (info .ne. 0)print*, 'pdgemm/pzgemm info: ',info
!                     print*, 'PE:',myrow,mycol,' C=',C
                    
                     ! Now do the matrix*v multiplication
                     ! call pdgemv(transa,transb,m,n,k,alpha,a,ia,ja,desc_a,b,ib,jb,desc_b,beta,
                     !              c,ic,jc,desc_c)
                     call pzgemv('n',sys_order_m,sys_order_n,alpha,C,iC,jC,desc_C,x,ix,jx,desc_x,1,&
                                 beta,y,ix,jx,desc_y,1)
                     if (info .ne. 0)print*, 'pzgemv info: ',info

!                     print*, 'PE:',myrow,mycol,' y=',y
        end subroutine mm_v
        
        subroutine mmdag_v(A,B,C,x,y,alpha,beta,sys_order_m,sys_order_n, &
                        desc_A,desc_B,desc_C,desc_x,desc_y)

                     ! Does matrix * matrix*vector using pzgemv and pzgemm
                     ! Uses the hermitian conjugate for matrix B instead.
                     ! A(m,k), B(n,k) --> C(m,n) 

!                     ! Real arrays
!                     real(kind=wp), dimension(:,:), intent(in) :: A
!                     real(kind=wp), dimension(:,:), intent(in) :: B
!                     real(kind=wp), dimension(:,:), intent(inout) :: C
!                     real(kind=wp), dimension(:,:), intent(in) :: x
!                     real(kind=wp), dimension(:,:), intent(out) :: y
!                     real(kind=wp) :: alpha,beta

                     ! Complex arrays
                     complex(kind=wp), dimension(:,:), intent(in) :: A
                     complex(kind=wp), dimension(:,:), intent(in) :: B
                     complex(kind=wp), dimension(:,:), intent(inout) :: C
                     complex(kind=wp), dimension(:,:), intent(in) :: x
                     complex(kind=wp), dimension(:,:), intent(out) :: y
                     complex(kind=wp) :: alpha,beta

                     integer, intent(in) :: sys_order_m,sys_order_n
                     integer, dimension(9), intent(in) :: desc_A,desc_B,desc_C,desc_x,desc_y
                     integer :: iA,jA,iB,jB,iC,jC,ix,jx,nrhs
                     integer :: info
                     integer, dimension(sys_order_m+block_size_m) :: ipvt

                     ! Using the whole matrices and vectors
                     iA=1; jA=1; iB=1; jB=1; iC=1; jC=1; ix=1; jx=1
                     
!                     print *
!                     print *, " alpha = ",alpha," beta = ",beta
!                     print *, "sys_order_m =",sys_order_m,"sys_order_n =",sys_order_n
!                     print *
!                     print *,"desc_C ="
!                     print *,desc_C
!                     print *
!                     print *,"desc_y ="
!                     print *,desc_y
!                     print *
                     ! Holds up execution of all processes within the indicated scope
                     ! seems to be present ony with the INTEL BLACS
!                     call blacs_barrier(icontxt, 'ALL')

                     ! pdgemm routine
                     ! Routine to do general matrix times matrix operation
                     ! C = alpha*A*B + beta*C,  A(m,k), B(k,n) --> C(m,n)
                     ! transa - form of matrix A used, 'N' for normal, 'T' for transpose
                     ! 'C' for conjugate transpose
                     ! transb - form of matrix B used
                     ! m - the num of rows of global matrix C, if transa='N' its num of rows in A
                     ! n - the num of cols of global matrix C, if transb='N' its num of cols in B
                     ! k- if transa='N',the num of cols of global matrix A, otherwise num of cols in A.
                     !     If transb='N' its num of rows in B, otherwise, num of cols in B
                     ! alpha - global scalar alpha in equation above.
                     ! a - local part of global matrix A, also identifies 1st elem of local array
                     ! iA - the row index of global matrix A, identifying the first row to use in calc
                     ! jA - the col index of global matrix A, indentifying the first col to use
                     ! desc_a - descriptor for matrix a
                     ! b - local part of global matrix B, also identifies 1st elem of local array
                     ! iC - the row index of global matrix B, identifying the first row to use in calc
                     ! jC - the col index of global matrix B, indentifying the first col to use
                     ! desc_b - descriptor for matrix b
                     ! beta - global scalar beta in equation above.
                     ! c - local part of global matrix B, also identifies 1st elem of local array
                     ! iC - the row index of global matrix C, identifying the first row to use in calc
                     ! jC - the col index of global matrix C, indentifying the first col to use
                     ! desc_c - descriptor for matrix c
                     !
                     ! ON RETURN
                     ! a - the updated local part of global matrix A, thats factorized
                     ! ipvt - the local part of the global vector ipvt, which contains the pivot indices
                     ! info - if 0, factorization is succesful, otherwise A is singular, global

                     ! call pdgemm(transa,transb,m,n,k,alpha,a,ia,ja,desc_a,b,ib,jb,desc_b,beta,
                     !              c,ic,jc,desc_c)

                     call pzgemm('n','c',sys_order_m,sys_order_m,sys_order_n,alpha,A,iA,jA,desc_A,B,&
                                 iB,jB,desc_B,beta,C,iC,jC,desc_C)
!                     if (info .ne. 0)print*, 'pdgemm/pzgemm info: ',info
!                     print*, 'PE:',myrow,mycol,' C=',C
                    
                     ! Now do the matrix*v multiplication
                     ! call pdgemv(transa,transb,m,n,k,alpha,a,ia,ja,desc_a,b,ib,jb,desc_b,beta,
                     !              c,ic,jc,desc_c)
!                     call pzgemv('n',sys_order_m,sys_order_n,alpha,C,iC,jC,desc_C,x,ix,jx,desc_x,1,&
!                                 beta,y,ix,jx,desc_y,1)
                     ! For doing transb = 'C'
                     call pzgemv('n',sys_order_m,sys_order_m,alpha,C,iC,jC,desc_C,x,ix,jx,desc_x,1,&
                                 beta,y,ix,jx,desc_y,1)
!                     if (info .ne. 0)print*, 'pzgemv info: ',info

!                     print*, 'PE:',myrow,mycol,' y=',y
        end subroutine mmdag_v

        subroutine zcheck_orthog_frob(A,B,m,k,frob)
!  Calculates the Frobenius norm of a complex matrix manually.
        implicit none
        integer :: i,j,k,m,f
        real(kind=wp) :: frob,fterm
        complex(kind=wp) :: A(m,k),B(m,k) 
        complex(kind=wp),allocatable :: iden(:,:) 

        allocate(iden(k,k))
        iden = 0.0_wp
!          write(6,*) "shape Pj=",shape(Pj)
!          write(6,*) " "  
!  Set the identity matrix
        do f = 1,k
          iden(f,f) = 1.0_wp
        enddo    
        iden = iden-matmul(conjg(transpose(A)),B)

!  Start calc of frobenius norm    
        fterm = 0.0_wp
        do f=1,k
          do j=1,k
             fterm = fterm + (abs(iden(f,j)))**2
!                 write(6,*) "i=",f,"j=",j,"fterm_ij=",fterm
!                 write(6,*) " "
          enddo
        enddo 
        frob=sqrt(fterm) 
        write(6,*) " " 
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "frob=",frob
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "
        deallocate(iden)
        return
        end subroutine zcheck_orthog_frob
        
        subroutine lanczos_non_herm(A,m,mlanc,zalpha,zbeta,zgamma,& 
                                   zborth,fnorm,paige,biortho)
        implicit none
!--------------------------------------------------------------------------------
! Follows the non-hermitian eigen vaalue problem algorithm from 
! https://netlib.org/utk/people/JackDongarra/etemplates/node245.html#bd:lanalg
! Actual algorithm only uses two recurrence terms
!-------------------------------------------------------------------------------
        integer :: i,j,k,l,m,n,f,mlanc,mhalf
        real(kind=wp) :: b0,temp1,tol
        real(kind=wp) :: u0,csgn,rand_num(m),fnorm(mlanc),fterm,frob
        real(kind=wp),parameter :: pi=3.141592654d0
        complex(kind=wp),parameter :: ii=cmplx(0.0_wp,1.0_wp)
        integer,parameter :: sedsz=33
        integer :: seed(sedsz) 
        complex(kind=wp) ::zalpha(mlanc),zbeta(mlanc),zgamma(mlanc), &
                     &  lnorm1,lnorm2,term,term1,zborth(mlanc),tleft,&
                     &   tright
        complex(kind=wp):: pm1(m),qm1(m),p(m),q(m),s(m),tq(m),tp(m), &
                      & pp(m),qp(m) ! u is now vp/ v^_(j+1) ,wp is
                                    !  w^_(j+1)
        
        complex(kind=wp) :: A(m,m)
        logical :: paige, biortho  ! to decide whether or not to do a Paige like step

!  Matrices for calculating biorthogonality
        complex(kind=wp) :: Planc(m,mlanc), Qlanc(m,mlanc)
        complex(kind=wp), allocatable :: iden(:,:),Pj(:,:),Qj(:,:), &
                                        iden1(:,:),iden2(:,:)
! ScaLAPACK variable declarations
        integer :: numroc,indxl2g
        integer :: desc_Pj(9),desc_Qj(9),desc_Pj_loc(9),desc_Qj_loc(9),desc_T_Pj_loc(9),&
                   desc_T_Qj_loc(9),desc_pp(9),desc_qp(9),&
                   desc_pp_loc(9),desc_qp_loc(9),desc_pp_orth_loc(9),desc_qp_orth_loc(9),&
                   pp_orth(m),qp_orth(m)
        integer :: block_size_thresh,block_size_opt,block_size_m,block_size_k,&
                   mloc_Pj,nloc_Pj,mloc_pp,nloc_pp,i_loc,j_loc,icontxt,info,lld,&
                   lld_g,lld_loc,rsrc,csrc
        complex(kind=wp),allocatable :: Pj_loc(:,:),Qj_loc(:,:),T_Pj_loc(:,:),T_Qj_loc(:,:),&
                                        pp_loc(:,:),qp_loc(:,:),pp_orth_loc(:,:),qp_orth_loc(:,:)
        complex(kind=wp) :: alpha,beta
!   Allocating the arrays needed: naming convention follows Y Saad
!    paper, "The Lanczos Biorthogonalization Algorithm and other Oblique
!    Projection Methods For Solving Large Unsymmetric Systems. " 
!   vqm1 = q_(j-1),v=v_j, tq=A*q_j, tp=A^H*p_j, qp=p^_(j+1)=A*q_j - gamma_j*q_(j-1)
!   pp = p^_(j+1)=A^H*p_j-conjg(alpha_j)*p_j-conjg(beta_j)*p_(j-1).
!   Adapting Saad's method to the compex case and choosing beta*delta
!   accordinf to the yambo paper "Implementation and testing of
!   Lanczos-based algorithms for Random-Phase Approximation
!   eigenproblems."
 
!        allocate(A(m,m))
!        allocate(alpha(mlanc),beta(mlanc))
!        allocate(vm1(m),v(m),s(m),t(m),u(m) ! the temp vectors
        
!  Printing out matrix elements
!        do i=1,m
!          do j=1,m
!             write(6,*) ""
!             write(6,*) "--- i = ",i,"j = ",j,"A(i,j) = ",A(i,j)
!             write(6,*) ""
!          enddo
!        enddo
!        write(6,*) "||A||_F =",norm2(real(A,kind=8))
!   Initializing arrays to 0

        ! Read in the values related to block_size
        ! block_size_thresh block_size_opt block_size_m

        open(10,file='block.in',status='old')
        rewind(10)
        read(10,*) block_size_thresh,block_size_opt,block_size_m
        close(10)

        zalpha=0.0_wp
        zbeta=0.0_wp
        zgamma=0.0_wp
        zborth = 0.0_wp
        fnorm = 0.0_wp
        pm1=0.0_wp
        qm1=0.0_wp
        p=0.0_wp
        q=0.0_wp
        s=0.0_wp
        tp=0.0_wp
        tq=0.0_wp
        pp=0.0_wp
        qp=0.0_wp
        rand_num=0.0_wp
        mhalf=int(m/2)
!        one=1.0

!  Generating the test vector , v_1
        call random_seed(put=seed(1:sedsz))
        call random_number(rand_num)
!        s(1)=1.0d0+1.0*ii
!        s(mhalf)=1.0d0+1.0*ii
!        s=1.0d0+1.0*ii
        s=1.0_wp*rand_num+0.5_wp*ii*rand_num
!        write(6,*) ""
!        write(6,*) "--- The vector u is ",""
!        do i=1,m
!          write(6,*) "u(i) = ",u(i)
!        enddo
!        u0=DZNRM2(m,u,1)
!        u0=ZDOTC(m,u,1,u,1)
        u0=dot_product(s,s)
        u0=sqrt(u0)
!        write(6,*) "--- The Euclidean norm =",u0,u0**2
!        write(6,*) "--- The Euclidean norm =",u0
        p=s/u0
        q=p
!        do i=1,m
!          write(6,*)""
!          write(6,*) "v(i) = ",v(i)
!          write(6,*) ""
!        enddo
!        s(1)=u(1)/u0
!        write(6,*) "--- s(1) = ",s(1),"u(1) = ",u(1)

!  Collecting q_1 and p_1 to add to P and Q
!        rsrc = 0
!        csrc = 0
        if (myrow .eq. 0 .and. mycol .eq. 0) then
          Planc = 0.0_wp
          Qlanc = 0.0_wp
          Planc(:,1) = p
          Qlanc(:,1) = q
        endif ! Only processor (0,0) handles arrays and does computations

!----------  Starting the Lanczos iteration ------------------------
        do k=1,mlanc
          if (k==1) then
             if (myrow .eq. 0 .and. mycol .eq. 0) then
               zbeta(k)=0.0_wp
               zgamma(k)=0.0_wp
             endif ! Only processor (0,0) handles arrays and does computations
          end if 
          write(6,*) "******** m = ",m," Lanczos iteration =",k,"***********"

!  Checking biorthogonality with the frobenius norm
!          allocate(iden(k,k),Pj(m,k),Qj(m,k))

          if (myrow .eq. 0 .and. mycol .eq. 0) then
            allocate(Pj(m,k),Qj(m,k))
!          iden = 0.0
            Pj = 0.0_wp
            Qj = 0.0_wp
            Pj = Planc(:,1:k)
            Qj = Qlanc(:,1:k)   
!          write(6,*) "Pj=",Pj
!          write(6,*) " "
!          write(6,*) "Qj=",Qj
!          write(6,*) " "  
            call zcheck_orthog_frob(Pj,Qj,m,k,frob)
!!  Set the identity matrix
!          do f = 1,k
!            iden(f,f) = 1.0d0
!          enddo    
!          iden = iden-matmul(conjg(transpose(Pj)),Qj)
!!          write(6,*) "iden=",iden
!!          write(6,*) " "
!!  Start calc of frobenius norm    
!          fterm = 0.0
!          do f=1,k
!            do j=1,k
!                 fterm = fterm + (abs(iden(f,j)))**2
!!                 write(6,*) "i=",f,"j=",j,"fterm_ij=",fterm
!!                 write(6,*) " "
!            enddo
!          enddo 
!          fterm=sqrt(fterm)     
!            write(6,*) " "
!            write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^" 
!            write(6,*) "frob =",frob !,"alt="
!!          write(6,*) sqrt(norm2(real(matmul(conjg(iden),iden))))
!            write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^"
!            write(6,*) " "
            fnorm(k) = frob
            deallocate(Pj,Qj)
          endif ! Only processor(0,0) does computations

!          deallocate(iden,Pj,Qj)
!          write(6,*)"alpha=",alpha(k),"beta =",beta(k),"delta=",delta(k)
!          write(6,*) ""
!          term=dot_product(p,q) ! Following the SLEPc paper for Hermitian casep
!          write(6,*) " k=1",term
          if (myrow .eq. 0 .and. mycol .eq. 0) then
            tq=matmul(A,q) !doing A*v_j
            tp=matmul(transpose(conjg(A)),p) ! doing A^H*w_j
            term=dot_product(p,tq) ! Following the SLEPc paper for Hermitian case
!            write(6,*) ""
!            write(6,*) "====== Before Lanczos, overlap ========"
!            write(6,*) " <p|A|q> =", term
!            write(6,*) ""
!   Including this Paige-like regularization, from hermitian lanczos, for better convergence
!          if (paige) then
            do i=1,m
!            write(6,*) ""
!            write(6,*) " --- i = ",i,"------"
!            write(6,*) "before updating, t(i)=",t(i)
              qp(i)=tq(i)-zgamma(k)*qm1(i) ! v^_(j+1)
              pp(i)=tp(i)-conjg(zbeta(k))*pm1(i) ! v^_(j+1)
!            write(6,*) ""
!            write(6,*) " After updating, t(i) ="
!            write(6,*) ""
!           write(6,*)"u_(j+1)(i)=",u(i)
!            write(6,*) ""
!            write(6,*) "s(i)=",s(i),"r(i)=",r(i)
!            write(6,*) ""
            enddo

            term=dot_product(p,qp) ! Following the SLEPc paper for Hermitian case
            zalpha(k) = term
!            write(6,*) ""
!            write(6,*) "alpha(k) =",zalpha(k)
!          write(6,*) "After calculating alpha, coeffs are "
!             write(6,*) ""
!             write(6,*)'alpha=',zalpha(i)
!             write(6,*)""
!             write(6,*)'beta(i)=',zbeta(i)
!             write(6,*)""
!             write(6,*)'delta(i)=',zdelta(i)
!             write(6,*) ""
!          do i=1,m
!            write(6,*) ""
!            write(6,*) " --- i = ",i,"------"
!            write(6,*) "before updating, t(i)=",t(i)
            if (paige) then
              do i=1,m
                qp(i)=qp(i)-zalpha(k)*q(i) ! u_(j+1)  Paige-like regularization
                pp(i)=pp(i)-conjg(zalpha(k))*p(i) ! u_(j+1), Paige-like reg
              enddo
            else
              do i=1,m
                qp(i)=tq(i)-zalpha(k)*q(i) ! u_(j+1) , netlib Alg 7.42
                pp(i)=tp(i)-conjg(zalpha(k))*p(i) ! u_(j+1), netlib Alg 7.42
              enddo
            endif
          endif ! Only processor(0,0) does computations

!            write(6,*) ""
!            write(6,*) " After updating, t(i) ="
!            write(6,*) ""
!            write(6,*)"u_(j+1)(i)=",u(i)
!            write(6,*) ""
!            write(6,*) "--- the result of daxpy is "
!            write(6,*) ""
!          enddo

!!  Re-biorthogonalizing the Lanczos vectors, full reorthog
          if (biortho) then  
            if (myrow .eq. 0 .and. mycol .eq. 0) then
              print *
              print *," Startin re-orthogonalization process."
              print *
            endif ! Only processor(0,0) does computations
                  
             ! Doing matrix*matrix*vector
             alpha = 1.0_wp
             beta = 0.0_wp
!!  Alternative way to do Gram-Schmidt with (I-Q_k*(P_k)^H)
!            allocate(iden1(m,m),iden2(m,m),Pj(m,k),Qj(m,k))
!            iden1 = 0.0
!            iden2 = 0.0
!            Pj = 0.0
!            Qj = 0.0
!            Pj = Planc(:,1:k)
!            Qj = Qlanc(:,1:k)   
!!          write(6,*) "Pj=",Pj
!!          write(6,*) " "
!!          write(6,*) "Qj=",Qj
!!          write(6,*) " "  
!!  Set the identity matrix
!            do f = 1,m
!              iden1(f,f) = 1.0d0
!              iden2(f,f) = 1.0d0
!            enddo  
!
!! For Hermitian case, works well  
!            iden1 = iden1-matmul(Qj,conjg(transpose(Pj)))
!            iden2 = iden2-matmul(Pj,conjg(transpose(Qj)))
!
!! re biorthogonalizing
!            qp = matmul(iden1,qp)
!            pp = matmul(iden2,pp)
!            deallocate(iden1,iden2,Pj,Qj)
           
            ! Alternative formulation, using ScaLAPACK multiplication
            ! Get the current Lanczos vectors, global arrays
            if (myrow .eq. 0 .and. mycol .eq. 0) then
              allocate(Pj(m,k),Qj(m,k))
              Pj = 0.0_wp
              Qj = 0.0_wp
              Pj = Planc(:,1:k)
              Qj = Qlanc(:,1:k)
            endif ! ----- Only processor (0,0)  Does the Lanczos iterations


            ! Determine block_size_k, since this changes with iteration.
            if(k .le. block_size_thresh) then
               block_size_k = k
               print *
               print *, " block_size_k =",block_size_k
            else
               block_size_k = block_size_opt
               print *
               print *, " block_size_k =",block_size_k
            endif

            ! Get the local sizes owned by each processor 
            ! Pj and Qj will have the same distribution of submatrices,
            ! and similarly, pp and qp will have same block distribution.
            
            mloc_Pj = numroc(m,block_size_m,myrow,0,nprow)
            nloc_Pj = numroc(k,block_size_k,mycol,0,npcol)
            mloc_pp = numroc(m,block_size_m,myrow,0,nprow)
            nloc_pp = numroc(1,1,mycol,0,npcol)
            lld_g = max(numroc(m,m,myrow,0,nprow),1)
            lld = max(mloc_Pj,1)
            lld_loc = max(lld_g,lld)
            print *
            print *,"myrow=",myrow,"mycol=",mycol
            print *
            print *,"mloc_Pj=",mloc_Pj,"nloc_Pj=",nloc_Pj
            print *
            print *,"mloc_pp=",mloc_pp,"nloc_pp=",nloc_pp
            print *
            print *,"lld_g=",lld_g,"lld=",lld,"lld_loc=",lld_loc
            print *

            ! Set up descriptors for the global matrices and vectors,existing only on
            ! processor(0,0)
            call descinit(desc_Pj,m,k,m,k,0,0,&
                          icontxt,lld_g,info)
!            print *
!            print *,"parameters: desc_Pj, m, k, m, k, rsrc, csrc, icontxt, lld, info"
!            print *
!            print *,"desc_Pj",m,k,m,k,desc_Pj(7),desc_Pj(8),icontxt,lld_g,info
!            print * 
!            print *,"desc_Pj"
!            print *
            print *,desc_Pj
            print *

!            call descinit(desc_Qj,m,k,m,k,0,0,&
!                          icontxt,m,info) 

!            call descinit(desc_pp,m,1,m,1,0,0,&
!                          icontxt,m,info) 
!            call descinit(desc_qp,m,1,m,1,0,0,&
!                          icontxt,m,info) 
!
!            ! Set up descriptors for the local matrices and vectors, showing how
!            ! they'll be distributed from the global arrays. 
!            ! The temporary product matrix Tj = Pj*Qj^H
!            call descinit(desc_Pj_loc,m,k,block_size_m,block_size_k,rsrc,csrc,&
!                          icontxt,lld,info) 
!            print *
!            print *,"parameters: desc_Pj, m, k, m, k, rsrc, csrc, icontxt, lld_loc, info"
!            print *
!            print *,"desc_Pj",m,k,m,k,rsrc,csrc,icontxt,lld_loc,info
!            print * 
!            print *,"desc_Pj"
!            print *,desc_Pj
!            print *
!            call descinit(desc_Qj_loc,m,k,block_size_m,block_size_k,0,0,&
!                          icontxt,mloc_Pj,info) 
!            call descinit(desc_T_Pj_loc,m,m,block_size_m,block_size_m,0,0,&
!                          icontxt,mloc_Pj,info) 
!            call descinit(desc_T_Qj_loc,m,m,block_size_m,block_size_m,0,0,&
!                          icontxt,mloc_Pj,info) 
!
!            ! Vectors
!            call descinit(desc_pp_loc,m,1,block_size_m,1,0,0,&
!                          icontxt,mloc_Pj,info) 
!            call descinit(desc_qp_loc,m,1,block_size_m,1,0,0,&
!                          icontxt,mloc_Pj,info) 
!
!            ! And for the orthogonalized vecs, local part
!            call descinit(desc_pp_orth_loc,m,1,block_size_m,1,0,0,&
!                          icontxt,mloc_Pj,info) 
!            call descinit(desc_qp_orth_loc,m,1,block_size_m,1,0,0,&
!                          icontxt,mloc_Pj,info) 

            ! Allocate the local arrays, remember vectors are also (m,1) array
            ! Matrices
            allocate(Pj_loc(mloc_Pj,nloc_Pj),Qj_loc(mloc_Pj,nloc_Pj),&
                     T_Pj_loc(mloc_Pj,mloc_Pj),T_Qj_loc(mloc_Pj,mloc_Pj))

            ! vectors
            allocate(pp_loc(mloc_pp,nloc_pp),qp_loc(mloc_pp,nloc_pp),&
                     pp_orth_loc(mloc_pp,nloc_pp),qp_orth_loc(mloc_pp,nloc_pp))

            ! Initialize the values of the local arrays and vectors from the global arrays

            ! Matrices
            Pj_loc = 0.0_wp
            Qj_loc = 0.0_wp
            T_Pj_loc = 0.0_wp
            T_Qj_loc = 0.0_wp
!!            do i_loc=1,mloc_Pj
!!              do j_loc=1,nloc_Pj
!!                i = indxl2g(i_loc,block_size_m,myrow,0,nprow)
!!                j = indxl2g(j_loc,block_size_n,mycol,0,npcol)
!!                Pj_loc(i_loc,j_loc) = Pj(i,j)
!!                Qj_loc(i_loc,j_loc) = Qj(i,j)
!!              enddo
!!            enddo

            ! Vectors
            pp_loc = 0.0_wp
            qp_loc = 0.0_wp
            pp_orth = 0.0_wp
            qp_orth = 0.0_wp

            ! Now have processor (0,0) distribute Pj and Qj to the Pj_loc
            ! pzgeadd - performs sum operation for two distributed general matrices
            ! sub(C) := beta*sub(C) + alpha*op(sub(A))
            ! call pzgeadd(trans,m,n,alpha,a,ia,ja,desc_aa,beta,c,ic,jc,desc_c)
!!            if (myrow .eq. 0 .and. mycol .eq. 0) then
!            call pzgeadd('n',m,k,alpha,Pj,1,1,desc_Pj,beta,Pj_loc,1,1,desc_Pj_loc)
!            call pzgeadd('n',m,k,alpha,Qj,1,1,desc_Qj,beta,Qj_loc,1,1,desc_Qj_loc)
!
!            call pzgeadd('n',m,1,alpha,pp,1,1,desc_pp,beta,pp_loc,1,1,desc_pp_loc)
!            call pzgeadd('n',m,1,alpha,qp,1,1,desc_qp,beta,qp_loc,1,1,desc_qp_loc)
!!            endif ! ----- Only processor (0,0)  Does the Lanczos iterations

!!            do j_loc=1,nloc_pp
!!              j = indxl2g(j_loc,1,mycol,0,npcol)
!!              do i_loc=1,mloc_pp
!!                i = indxl2g(i_loc,block_size_m,myrow,0,nprow)
!!                pp_loc(i_loc,j_loc) = pp(i)
!!                qp_loc(i_loc,j_loc) = qp(i)
!!               enddo
!!             enddo


             ! For pp_orth Pj*Qj^H
             !subroutine mmdag_v(A,B,C,x,y,alpha,beta,sys_order_m,sys_order_n, &
             !           sys_order_k,desc_A,desc_B,desc_C,desc_x,desc_y)

!             call mmdag_v(Pj_loc,Qj_loc,T_Pj_loc,pp_loc,pp_orth_loc,alpha,beta,m,m,k,&
!                          desc_Pj_loc,desc_Qj_loc,desc_T_Pj_loc,desc_pp_loc,desc_pp_orth_loc)
!
!             ! For qp_orth Qj*Pj^H
!             call mmdag_v(Qj_loc,Pj_loc,T_Qj_loc,qp_loc,qp_orth_loc,alpha,beta,m,m,k,&
!                          desc_Pj_loc,desc_Qj_loc,desc_T_Qj_loc,desc_qp_loc,desc_qp_orth_loc)
!
!             ! Collect all the local parts of pp_orth and qq_orth
!             ! and then have them on processor (0,0)
!             call descinit(desc_pp_orth,m,1,m,1,0,0,icontxt,m,info)
!             call descinit(desc_qp_orth,m,1,m,1,0,0,icontxt,m,info)
!              
!             call pzgemr2d(m,1,pp_orth_loc,1,1,desc_pp_orth_loc,pp_orth,&
!                           1,1,desc_pp_orth,icontxt)
!             call pzgemr2d(m,1,qp_orth_loc,1,1,desc_qp_orth_loc,qp_orth,&
!                           1,1,desc_qp_orth,icontxt)

!!             ! Have processor (0,0) broadcast the global vectors pp_orth,qq_orth to
!!             ! all procesors in the grid.
!!             if(myrow .eq. 0 .and. mycol .eq. 0) then
!!               call zgebs2d(icontxt,"All"," ",m,1,pp_orth,m)
!!               call zgebs2d(icontxt,"All"," ",m,1,qp_orth,m)
!!             else
!!               call zgebr2d(icontxt,"All"," ",m,1,pp_orth,m,0,0)
!!               call zgebr2d(icontxt,"All"," ",m,1,qp_orth,m,0,0)
!!             endif
            
            ! Deallocate arrays 
            deallocate(Pj_loc,Qj_loc,T_Pj_loc,T_Qj_loc)
            deallocate(pp_loc,qp_loc,pp_orth_loc,qp_orth_loc)

            if(myrow .eq. 0 .and. mycol .eq. 0) then
              deallocate(Pj,Qj)
            ! vectors
             ! Finish the orthogonalization process
              print *
              print *," Finished the re-orthogonalization process."
              print *
              pp = pp !- pp_orth
              qp = qp !- qp_orth
            endif ! ----- Only processor (0,0)  Does the Lanczos iterations

          endif
!!  New norm
!          term = dot_product(pp,pp) ! following netlib Alg 7.13
!          pp = (1/sqrt(term))*pp ! renormalizing s, so it doesnt get too large
          if(myrow .eq. 0 .and. mycol .eq. 0) then
            term1 = dot_product(qp,pp) ! following netlib Alg 7.13

!  A measure of the biorthogonality
            zborth(k) = term1

!            term = dot_product(qp,qp) ! following netlib Alg 7.13
!            write(6,*) ""
!            write(6,*) "(r^_(j+1),r^_(j+1))=",term
!            write(6,*) "norm=",sqrt(term) 
!            write(6,*) ""
!            term = dot_product(pp,pp) ! following netlib Alg 7.13
!            write(6,*) ""
!            write(6,*) "(s^_(j+1),s^_(j+1))=",term
!            write(6,*) "norm=",sqrt(term) 
!            write(6,*) ""
!            write(6,*) ""
!            write(6,*) "(r^_(j+1),s^_(j+1))=",term1 
!            write(6,*) ""
            lnorm1 = sqrt(abs(term1))  ! following yambo paper
!          lnorm1 = sqrt(abs(term1*conjg(term1)))  ! following yambo paper
            lnorm2 = conjg(term1)/lnorm1
!            write(6,*) "|(r^_(j+1),s^_(j+1))|=",lnorm1 !,"sign =",csgn
!!          write(6,*) "(v^_(j+1),w^_(j+1))=",lnorm1
!            write(6,*) "lnorm2 =",lnorm2
!            write(6,*) ""
            if (k .ne. mlanc) then
              zbeta(k+1)= lnorm1 ! following yambo paper
              zgamma(k+1)=lnorm2 ! following Y. Saad's paper
!            write(6,*) ""
!            write(6,*) "for iteration ",k,"beta =",lnorm
!            write(6,*) ""
            endif

!  Updating v_(j-1) for the next run (j+1) 
            do i=1,m
              qm1(i)=q(i)
              pm1(i)=p(i)
            enddo

!  Updating v_(j+1), to be used at start of (j+1) iteration, can use the tol here also
            if (k .ne. mlanc) then
              do i=1,m
!            s(i)=(1/beta(k+1))*s(i)
                q(i)=(1.0_wp/lnorm1)*qp(i)! v_(j+1)=v^_(j+1)/delta_(j+1)
                p(i)=(1.0_wp/conjg(lnorm2))*pp(i) ! w_(j+1)=w^_(j+1)/beta_(j+1)
              enddo

!            write(6,*) "======= Norm of updated lanczos vectors======"
!            write(6,*) " q =",q
!            write(6,*) " "
!            write(6,*)"||q||=",sqrt(dot_product(q,q))
!            write(6,*) " "
!            write(6,*) " p =",p
!            write(6,*) " " 
!            write(6,*) "||p||=",sqrt(dot_product(p,p))
!            write(6,*) " "
!  Add new vectors 
              Planc(:,k+1) = p
              Qlanc(:,k+1) = q
            endif
          endif ! ----- Only processor (0,0)  Does the Lanczos iterations

!  Checking the local biorthogonality
!          term = dot_product(q,p)
!          write(6,*) "========== biorthogonality ========="
!          write(6,*) " < p_i, q_i > =",term
!          write(6,*) " "
        enddo ! END OF LANCZOS LOOP
        return
        end subroutine lanczos_non_herm
        
        subroutine lanc_non_herm_pro_1(A,m,mlanc,zalpha,zbeta,zgamma,& 
                                   zborth,fnorm,paige,biortho)
        implicit none
!--------------------------------------------------------------------------------
! Follows the non-hermitian eigen vaalue problem algorithm from 
! https://netlib.org/utk/people/JackDongarra/etemplates/node245.html#bd:lanalg
! Actual algorithm only uses two recurrence terms.
!
! Does the partial reorthogonalization on the unnormalized pp and qp
!-------------------------------------------------------------------------------
        integer :: i,j,k,l,m,n,f,mlanc,mhalf,iterm,jterm
        real(kind=8) :: b0,temp1,tol,eps1,term_re,term_im
        real(kind=8) :: u0,csgn,rand_num(m),fnorm(mlanc),fterm,frob
        real,parameter :: pi=3.141592654d0,eps=1.11E-16
        complex,parameter :: ii=cmplx(0.d0,1.d0)
        integer,parameter :: sedsz=33
        integer :: seed(sedsz),orthtot,poterm,qoterm,chk,oterm 
        complex(kind=8) ::zalpha(mlanc),zbeta(mlanc),zgamma(mlanc), &
                     &  lnorm1,lnorm2,term,term1,zborth(mlanc)
        complex(kind=8):: pm1(m),qm1(m),p(m),q(m),s(m),tq(m),tp(m), &
                      & pp(m),qp(m),qo(m),po(m),vp(m),vq(m),wp(m),wq(m) ! u is now vp/ v^_(j+1) ,wp is
                                    !  w^_(j+1)
        
        complex(kind=8) :: A(m,m)
        logical :: paige, biortho  ! to decide whether or not to do a Paige like step

!  Matrices for calculating biorthogonality
        integer, allocatable :: porthopt(:),qorthopt(:),porthindx(:),&
                                qorthindx(:)
        complex(kind=8) :: Planc(m,mlanc), Qlanc(m,mlanc)
        complex(kind=8), allocatable :: iden(:,:),Pj(:,:),Qj(:,:), &
                        iden1(:,:),iden2(:,:),&
                        pw_jk(:),qw_jk(:)

                                      

!   Allocating the arrays needed: naming convention follows Y Saad
!    paper, "The Lanczos Biorthogonalization Algorithm and other Oblique
!    Projection Methods For Solving Large Unsymmetric Systems. " 
!   vqm1 = q_(j-1),v=v_j, tq=A*q_j, tp=A^H*p_j, qp=p^_(j+1)=A*q_j - gamma_j*q_(j-1)
!   pp = p^_(j+1)=A^H*p_j-conjg(alpha_j)*p_j-conjg(beta_j)*p_(j-1).
!   Adapting Saad's method to the compex case and choosing beta*delta
!   accordinf to the yambo paper "Implementation and testing of
!   Lanczos-based algorithms for Random-Phase Approximation
!   eigenproblems."
 
!        allocate(A(m,m))
!        allocate(alpha(mlanc),beta(mlanc))
!        allocate(vm1(m),v(m),s(m),t(m),u(m) ! the temp vectors
        
!  Printing out matrix elements
!        do i=1,m
!          do j=1,m
!             write(6,*) ""
!             write(6,*) "--- i = ",i,"j = ",j,"A(i,j) = ",A(i,j)
!             write(6,*) ""
!          enddo
!        enddo
!        write(6,*) "||A||_F =",norm2(real(A,kind=8))
!   Initializing arrays to 0
        zalpha=0.0
        zbeta=0.0
        zgamma=0.0
        zborth = 0.0
        fnorm = 0.0
        pm1=0.0
        qm1=0.0
        p=0.0
        q=0.0
        s=0.0
        tp=0.0
        tq=0.0
        pp=0.0
        qp=0.0
        po = 0.0
        qo = 0.0
        rand_num=0.0
        mhalf=int(m/2)

! Set tolerance
        eps1 = sqrt(eps)
!        eps1 = sqrt(eps) !1000*eps

!  Generating the test vector , v_1
        call random_seed(put=seed(1:sedsz))
        call random_number(rand_num)
!        s(1)=1.0d0+1.0*ii
!        s(mhalf)=1.0d0+1.0*ii
!        s=1.0d0+1.0*ii
        s=1.0d0*rand_num+0.5*ii*rand_num
!        write(6,*) ""
!        write(6,*) "--- The vector u is ",""
!        do i=1,m
!          write(6,*) "u(i) = ",u(i)
!        enddo
!        u0=DZNRM2(m,u,1)
!        u0=ZDOTC(m,u,1,u,1)
        u0=dot_product(s,s)
        u0=sqrt(u0)
!        write(6,*) "--- The Euclidean norm =",u0,u0**2
!        write(6,*) "--- The Euclidean norm =",u0
        p=s/u0
        q=p
!        do i=1,m
!          write(6,*)""
!          write(6,*) "v(i) = ",v(i)
!          write(6,*) ""
!        enddo
!        s(1)=u(1)/u0
!        write(6,*) "--- s(1) = ",s(1),"u(1) = ",u(1)

!  Collecting q_1 and p_1 to add to P and Q
        Planc = 0.0
        Qlanc = 0.0
        Planc(:,1) = p
        Qlanc(:,1) = q
 
        write(6,*) " "
        write(6,*) " ****** eps1 =",eps1,"*******"
        write(6,*) " "
! total # of orthogonalizations
       orthtot = 0
!----------  Starting the Lanczos iteration ------------------------
        do k=1,mlanc
          if (k==1) then
             zbeta(k)=0.0
             zgamma(k)=0.0
          end if

          if (k .ne. 1)then
            pm1 = Planc(:,k-1)
            qm1 = Qlanc(:,k-1)
          else
            pm1 = 0.0
            qm1 = 0.0
          endif

          write(6,*) "******** Lanczos iteration =",k,"***********"

!  Checking biorthogonality with the frobenius norm
!          allocate(iden(k,k),Pj(m,k),Qj(m,k))
          allocate(Pj(m,k),Qj(m,k))
!          iden = 0.0
          Pj = 0.0
          Qj = 0.0
          Pj = Planc(:,1:k)
          Qj = Qlanc(:,1:k)   
!          write(6,*) "Pj=",Pj
!          write(6,*) " "
!          write(6,*) "Qj=",Qj
!          write(6,*) " "  
          call zcheck_orthog_frob(Pj,Qj,m,k,frob)
!!  Set the identity matrix
!          do f = 1,k
!            iden(f,f) = 1.0d0
!          enddo    
!          iden = iden-matmul(conjg(transpose(Pj)),Qj)
!!          write(6,*) "iden=",iden
!!          write(6,*) " "
!!  Start calc of frobenius norm    
!          fterm = 0.0
!          do f=1,k
!            do j=1,k
!                 fterm = fterm + (abs(iden(f,j)))**2
!!                 write(6,*) "i=",f,"j=",j,"fterm_ij=",fterm
!!                 write(6,*) " "
!            enddo
!          enddo 
!          fterm=sqrt(fterm)     
          write(6,*) " "
          write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^" 
          write(6,*) "frob =",frob !,"alt="
!          write(6,*) sqrt(norm2(real(matmul(conjg(iden),iden))))
          write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^"
          write(6,*) " "
          fnorm(k) = frob
          deallocate(Pj,Qj)
!          deallocate(iden,Pj,Qj)
!          write(6,*)"alpha=",alpha(k),"beta =",beta(k),"delta=",delta(k)
!          write(6,*) ""
!          term=dot_product(p,q) ! Following the SLEPc paper for Hermitian casep
!          write(6,*) " k=1",term
          tq=matmul(A,q) !doing A*v_j
          tp=matmul(transpose(conjg(A)),p) ! doing A^H*w_j
          term=dot_product(p,tq) ! Following the SLEPc paper for Hermitian case
          write(6,*) ""
          write(6,*) "====== Before Lanczos, overlap ========"
          write(6,*) " <p|A|q> =", term
          write(6,*) ""
!   Including this Paige-like regularization, from hermitian lanczos, for better convergence
!          if (paige) then
          do i=1,m
!            write(6,*) ""
!            write(6,*) " --- i = ",i,"------"
!            write(6,*) "before updating, t(i)=",t(i)
             qp(i)=tq(i)-zgamma(k)*qm1(i) ! v^_(j+1)
             pp(i)=tp(i)-conjg(zbeta(k))*pm1(i) ! v^_(j+1)
!            write(6,*) ""
!            write(6,*) " After updating, t(i) ="
!            write(6,*) ""
!           write(6,*)"u_(j+1)(i)=",u(i)
!            write(6,*) ""
!            write(6,*) "s(i)=",s(i),"r(i)=",r(i)
!            write(6,*) ""
          enddo

          term=dot_product(p,qp) ! Following the SLEPc paper for Hermitian case
          zalpha(k) = term
          write(6,*) ""
          write(6,*) "alpha(k) =",zalpha(k)
!          write(6,*) "After calculating alpha, coeffs are "
!             write(6,*) ""
!             write(6,*)'alpha=',zalpha(i)
!             write(6,*)""
!             write(6,*)'beta(i)=',zbeta(i)
!             write(6,*)""
!             write(6,*)'delta(i)=',zdelta(i)
!             write(6,*) ""
!          do i=1,m
!            write(6,*) ""
!            write(6,*) " --- i = ",i,"------"
!            write(6,*) "before updating, t(i)=",t(i)
          if (paige) then
            do i=1,m
              qp(i)=qp(i)-zalpha(k)*q(i) ! u_(j+1)  Paige-like regularization
              pp(i)=pp(i)-conjg(zalpha(k))*p(i) ! u_(j+1), Paige-like reg
            enddo
          else
            do i=1,m
              qp(i)=tq(i)-zalpha(k)*q(i) ! u_(j+1) , netlib Alg 7.42
              pp(i)=tp(i)-conjg(zalpha(k))*p(i) ! u_(j+1), netlib Alg 7.42
            enddo
          endif
!            write(6,*) ""
!            write(6,*) " After updating, t(i) ="
!            write(6,*) ""
!            write(6,*)"u_(j+1)(i)=",u(i)
!            write(6,*) ""
!            write(6,*) "--- the result of daxpy is "
!            write(6,*) ""
!          enddo


! Checking for and doing partial reorthogonality
          if (k .ne. mlanc) then
!            allocate(pw_jk(k),qw_jk(k))
            allocate(Pj(m,k),Qj(m,k),pw_jk(k),qw_jk(k))
!            allocate(iden1(m,m),iden2(m,m))
!            iden1 = 0.0
!            iden2 = 0.0
            Pj = 0.0
            Qj = 0.0
            pw_jk = 0.0
            qw_jk = 0.0
            Pj = Planc(:,1:k)
            Qj = Qlanc(:,1:k)
            write(6,*) "shape Pj/Qj =",shape(Pj)
            write(6,*) " "

!!  Calculate local orthog
            vp = 0.0
            vq = 0.0
            term1 = dot_product(qp,pp) ! following netlib Alg 7.13
            lnorm1 = sqrt(abs(term1))  ! following yambo paper
            lnorm2 = conjg(term1)/lnorm1
            vp=(1.0/conjg(lnorm2))*pp ! w_(j+1)=w^_(j+1)/beta_(j+1)
            vq=(1.0/lnorm1)*qp! v_(j+1)=v^_(j+1)/delta_(j+1)
            pw_jk = matmul(conjg(transpose(Qj)),vp)
            qw_jk = matmul(conjg(transpose(Pj)),vq)
!!            estnorm(k) = max(maxval(abs(pw_jk)),maxval(abs(pw_jk)))
            write(6,*) " "
            write(6,*) "orthog, <p_j+1,q_k>,max",maxval(abs(pw_jk))
            write(6,*) " "
            write(6,*) pw_jk
            write(6,*) " "
            write(6,*) "orthog, <q_j+1,p_k>,max",maxval(abs(qw_jk))
            write(6,*) " "
            write(6,*) qw_jk
            write(6,*) " "

! deallocate Pj, Qj
            deallocate(Pj,Qj)

!!  Check for partial orthog
            allocate(porthopt(k),qorthopt(k))
            porthopt = 0
            qorthopt = 0
            poterm = 0
            qoterm = 0
! P's
            do i=1,k
              term = pw_jk(i)
              term_re = abs(real(term))
              term_im = abs(aimag(term))
!              write(6,*) "i=",i,"term =",term
!              write(6,*) " "
              if ((term_re .ge. eps1).or.(term_im .ge. eps1)) then
                write(6,*) " "
                write(6,*) " vector",i,"lost orth, pw_jk=",term,eps1
                write(6,*) " "
                porthopt(i) = 1
              endif
            enddo
            poterm = sum(porthopt)

! Q's
            do i=1,k
              term = qw_jk(i)
              term_re = abs(real(term))
              term_im = abs(aimag(term))
              if ((term_re .ge. eps1).or. (term_im .ge. eps1)) then
                write(6,*) " "
                write(6,*) " vector",i,"lost orth, qw_jk=",term,eps1
                write(6,*) " "
                qorthopt(i) = 1
              endif
            enddo
            qoterm = sum(qorthopt)
            write(6,*) "---# vecs that lost biortho",poterm,qoterm,"---"
            write(6,*) " "
!            deallocate(porthopt,qorthopt)
!            deallocate(Pj,Qj,pw_jk,qw_jk) 
!
!! Get the index of these vectors
             if (poterm .ne. qoterm) then
!               write(6,*) " "
!               write(6,*) "poterm =",poterm,"qoterm=",qoterm
!               write(6,*) " "
!               write(6,*) "orig porthopt ="
!               write(6,*) " "
!               write(6,*) porthopt
!               write(6,*) " "
!               write(6,*) "orig qorthopt ="
!               write(6,*) " "
!               write(6,*) qorthopt
!               write(6,*) " "
               oterm = max(poterm,qoterm)
                
!  Make the indices the same if the num of orthog is diff
               if (poterm .ge. qoterm) then
                 qorthopt = qorthopt + (porthopt-qorthopt)
!                 write(6,*) " new qorthopt "
!                 write(6,*) " "
!                 write(6,*) qorthopt
!                 write(6,*) " "
               else
                 porthopt = porthopt + (qorthopt-porthopt)
!                 write(6,*) " new porthopt "
!                 write(6,*) " "
!                 write(6,*) porthopt
!                 write(6,*) " "
               endif
               poterm = oterm
               qoterm = oterm
             endif

             allocate(porthindx(poterm),qorthindx(qoterm))
             porthindx = 0
             qorthindx = 0
            
! P's      
            iterm = 1
            do i=1,k
              if(porthopt(i) .eq. 1) then
                porthindx(iterm) = i
                iterm = iterm + 1
              endif
            enddo
! Q's      
            iterm = 1
            do i=1,k
              if(qorthopt(i) .eq. 1) then
                qorthindx(iterm) = i
                iterm = iterm + 1
              endif
            enddo
            deallocate(porthopt,qorthopt)

! Carry out rebiorthog wrt those vectors that lost it
!            poterm = 2
!            qoterm=2

! Construct the Pj and Qj
            if ((poterm .ge. 1) .or. (qoterm .ge. 1)) then
               write(6,*) " doing partial 2-sided GS as biorthog failed"
               write(6,*) " "
               allocate(Pj(m,poterm),Qj(m,qoterm),iden1(m,m),iden2(m,m))
               iden1 = 0.0
               iden2 = 0.0
               Pj = 0.0
               Qj = 0.0
               do i=1,poterm
                 iterm = porthindx(i)
                 jterm = qorthindx(i)
                 Pj(:,i) = Planc(:,iterm)
                 Qj(:,i) = Qlanc(:,jterm)
               enddo
 
!               frob = 0.0

!  Set the identity matrix
               do f = 1,m
                  iden1(f,f) = 1.0d0
                  iden2(f,f) = 1.0d0
               enddo  

! For Hermitian case, works well  
               iden1 = iden1-matmul(Qj,conjg(transpose(Pj)))
               iden2 = iden2-matmul(Pj,conjg(transpose(Qj)))

! re biorthogonalizing
               qp = matmul(iden1,qp)
               pp = matmul(iden2,pp)
               deallocate(iden1,iden2,Pj,Qj)
            endif   !
!            deallocate(porthindx,qorthindx)

! reassign pp and qp
!            pp = wp
!            qp = wq
            deallocate(porthindx,qorthindx,pw_jk,qw_jk)
!            deallocate(iden1,iden2,Pj,Qj,pw_jk,qw_jk)
          endif  ! re biorthog condition

! 
          term1 = dot_product(qp,pp) ! following netlib Alg 7.13

!  A measure of the biorthogonality
          zborth(k) = term1

          term = dot_product(qp,qp) ! following netlib Alg 7.13
          write(6,*) ""
          write(6,*) "(r^_(j+1),r^_(j+1))=",term
          write(6,*) "norm=",sqrt(term) 
          write(6,*) ""
          term = dot_product(pp,pp) ! following netlib Alg 7.13
          write(6,*) ""
          write(6,*) "(s^_(j+1),s^_(j+1))=",term
          write(6,*) "norm=",sqrt(term) 
          write(6,*) ""
          write(6,*) ""
          write(6,*) "(r^_(j+1),s^_(j+1))=",term1 
          write(6,*) ""
          lnorm1 = sqrt(abs(term1))  ! following yambo paper
          lnorm2 = conjg(term1)/lnorm1
          write(6,*) "|(r^_(j+1),s^_(j+1))|=",lnorm1 !,"sign =",csgn
!          write(6,*) "(v^_(j+1),w^_(j+1))=",lnorm1
          write(6,*) "lnorm2 =",lnorm2
          write(6,*) ""
          if (k .ne. mlanc) then
            zbeta(k+1)= lnorm1 ! following yambo paper
            zgamma(k+1)=lnorm2 ! following Y. Saad's paper
!            write(6,*) ""
!            write(6,*) "for iteration ",k,"beta =",lnorm
!            write(6,*) ""
          endif

!!  Updating v_(j-1) for the next run (j+1) 
!          do i=1,m
!            qm1(i)=q(i)
!            pm1(i)=p(i)
!          enddo

!  Updating v_(j+1), to be used at start of (j+1) iteration, can use the tol here also
          if (k .ne. mlanc) then
            do i=1,m
!            s(i)=(1/beta(k+1))*s(i)
              q(i)=(1.0/lnorm1)*qp(i)! v_(j+1)=v^_(j+1)/delta_(j+1)
              p(i)=(1.0/conjg(lnorm2))*pp(i) ! w_(j+1)=w^_(j+1)/beta_(j+1)
            enddo

!  Add new vectors 
            Planc(:,k+1) = p
            Qlanc(:,k+1) = q

! Check for orthog after partial rebiorthog
            call zcheck_orthog_frob(Planc(:,1:k+1),Qlanc(:,1:k+1),& 
                                       m,k+1,frob)
            write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(6,*) "Lost semi orthogonality at lanc iter",k
            write(6,*) " "
            write(6,*) "Norm after GS = ",frob
            write(6,*) " "
            write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(6,*) " "
          endif

! update total orthogs
          orthtot = orthtot + poterm + qoterm
        enddo ! END OF LANCZOS LOOP
!  Print out final total orthogonalizations
        write(6,*) " "
        write(6,*) "****** TOTAL ORTHOGS =",orthtot,"******"
        write(6,*) " "
        write(6,*) " "
        return
        end subroutine lanc_non_herm_pro_1
end program planczos_minprbo






