        program lanczos_test
        implicit none

!
! Code to test the implementation of the Hermitian and non-Hermitian
! Lanczos algorithm for tri-diagonalization of a matrix. It carries out, full,
! and partial rebiorthogonalization strategies.
! Subsequently, the calculation of the diagonal matrix elements of the resolvent (w-H)^-1, in terms of
! continued fractions is done and compared to the exact solution.
! Follows the papers: "Lanczos methods in SLEPc"
!                     

        integer, parameter :: dp = kind(1.0d0)

        ! Matrix parameters
        integer, parameter :: m=1000,mhalf=int(m/2), mlanc=300
        integer :: i,j,k,l,n,info,itemp,prbo,F_j,orthtot,orthfreq
        real(kind=8) :: b0,temp1,tol,md, mc,mf,frob,tau,sgn_star
        real(kind=8) :: u0,DZNRM2,ZDOTC,lnorm,DLANGE,rand_num(m),rand_num2(m),eps1 
        real,parameter :: pi=3.141592654d0,eps=1.11E-16
        integer,parameter :: sedsz=9
        integer           :: seed(sedsz),seed2(sedsz+3),herm,mi,mj,chknorm(mlanc)
        complex,parameter :: ii=complex(0.d0,1.d0)
        logical :: reorth, paige, full
        character (10) :: mult_type,prbo_type

        ! Non allocatable        
        complex(kind=8) ::zalpha(mlanc),zbeta(mlanc+1),zdelta(mlanc+1), &
                &           zgamma(mlanc+1),zborth(mlanc),Planc(m,mlanc)
        real(kind=8) ::alpha(mlanc),beta(mlanc+1),delta(mlanc+1),&
                       borth(mlanc),d_j(mlanc+1),v0(m)
        complex(kind=8) :: z(m,m), &
                &     vm1(m), v(m),s(m),t(m),u(m),w(m),v_h(m),tmp_vec(m)

        complex(kind=8) :: hterm, c2,H(m,m),Hw(m,m),F(m,m),R(mhalf,mhalf), &
                            C(mhalf,mhalf),Lch(m,m),Iden(m,m),temp_mat(m,m),&
                            temp_mat1(m,m),zterm,zterm1,zterm2,zterm3,zterm4,zterm5,&
                            al,be,arg,rp,rm,rlim,zsum1,zsum2,Iden_lanc(mlanc,mlanc),&
                            H_diag(m,m),H_off_diag(m,m)

        ! Real arrays, allocatable 
        real(kind=8), allocatable ::Htest(:,:),fnorm(:),estnorm(:),& ! fnorm - frobenius norm
                                    HGS(:,:),orth_lev(:)             ! orth_lev - the level of orthogonality/biorthogonality
        real(kind=dp), allocatable :: tmp_array1(:,:),tmp_vec1(:),tmp_vec2(:)
        complex(kind=dp), allocatable :: ztmp_array1(:,:),ztmp_vec1(:),ztmp_vec2(:),ztmp_mat(:,:),ztmp_mat2(:,:),ztmp_mat3(:,:),&
                                         g_jn0_vec_p(:),g_jn0_vec_m(:),ztmp_vec3(:),ztmp_vec4(:,:),ztmp_vec5(:),ztmp_vec6(:)
        real(kind=8) :: mat_ij,Hmat(m,m)
        integer :: dim1, dim2,nnz
        character (10) ::  dum
 
        ! Arrays and variables for DGETRI, DGETRF
        integer :: lda, lwork,omega_max
        complex(kind=8), allocatable :: work(:),work1(:),resolv_exact(:,:),resolv_lanc(:,:),&
                                        inner_prod(:,:),tmp_array(:,:)
        integer, allocatable :: ipiv(:), ipiv1(:)
        real(kind=8) :: omega, eta,omega_scal,rterm,rterm1,rterm2,rterm3,spect_range,del_spect,vbnorm

        ! Choosing Matrix type
        herm=3  !! hermitian=1, non-symm,real= 2, nonherm =3

        ! Reorthog flags and parameters
        ! for prbo, 1 -> prbo 1, 2 -> prbo 2, minprbo 
        reorth = .true.
        full = .false.
        paige = .true.
        prbo_type = 'prbo 2'
        eps1 = sqrt(eps) ! tolerance for reorthogonalization
        print *
        print *, " reorth = ",reorth," full = ",full," paige = ",paige," PRBO type = ",prbo_type
        print *

        ! Ways of carrying out matrix-matrix-vector product for reorthogonalization
        ! options : pref - using explicit matrix mult, oneshot - using implicit matrix mult and doing it all in one go
        ! naive - doing it the simplest way with dot products.
        ! default case with mult_type = ' ' is 'pref'
        mult_type = ' '
        
        print *
        print *, " mult type = ",mult_type
        print *

        ! Initializing and assembling arrays 
        allocate(fnorm(mlanc),estnorm(mlanc))

        ! Fill in the identity matrix in Lanczos subspace and in Hamiltonian space
        do i=1,m
           Iden(i,i) = 1.0d0
        enddo

        do i=1,mlanc
           Iden_lanc(i,i) = 1.0d0
        enddo


        Hmat = 0.0d0
        H = 0.0d0
        Iden = 0.0d0
        Iden_lanc = 0.0d0
        temp_mat = 0.0d0
        temp_mat1 = 0.0d0
        Planc = 0.0d0
        v = 0.0d0
        v0 = 0.0d0

        alpha = 0.0d0
        beta = 0.0d0
        delta = 0.0d0
        zalpha = 0.0d0
        zbeta = 0.0d0
        zdelta = 0.0d0
        zgamma = 0.0d0
        fnorm = 0.0d0
        sgn_star = 1.0d0
        

        ! Generating the random initial vector , v_0 
        seed(1)=100312; seed(2)=4318992; seed(3)=21313984; seed(4)=856721; seed(5)=123332
        seed(6)=9012212; seed(7)=5123568; seed(8)=1023145; seed(9)=2133457

        rand_num = 0.0d0
        rand_num2 = 0.0d0
        call random_seed(put=seed(1:sedsz))
        call random_seed(put=seed2(1:sedsz+3))
        call random_number(rand_num)
        call random_number(rand_num2)

        do i=1,m
          v(i) = 1.0d0*rand_num(i)+0.5d0*ii*rand_num(i)
          v0(i) = rand_num(i)
!          v(i) = 1.0d0
!          v0(i) = 1.0d0
        enddo

        ! Dielectric function parameters

        ! Energy parameters
        eta = 0.1d0 ! For non Hermitian matrix, needs to be 0.0 as Lanczos coeffs are complex with generally nonzero imag parts
        omega_scal = 0.05d0
        omega_max =10 

        ! Spectral range (difference between maximum and minimum eigenvlaues) parameters.
        ! Realistic solid state Hamiltonians have narrow spectral ranges
        spect_range = 22.0d0
        del_spect = spect_range/(1.0d0*m)
        print *
        print *, " Spectral range = ",spect_range," step size = ",del_spect
        print *

        ! File to saving the Lanczos coefficients, loss of orthogonality, dielectric function to file

         ! hermitian
!        open(10,file='lanc_coeffs_herm_no_m_1000_mlanc_100_spec_range_22.dat')
!        open(14,file='fnorm_herm_no_m_1000_mlanc_100_spec_range_22.dat')
!        open(11,file='resolvent_exact_herm_m_1000_spec_range_22.dat')
!        open(111,file='resolvent_lanc_herm_no_m_1000_mlanc_100_spec_range_22.dat')

!        open(10,file='lanc_coeffs_herm_full_m_1000_mlanc_100_spec_range_22.dat')
!        open(14,file='fnorm_herm_full_m_1000_mlanc_100_spec_range_22.dat')
!        open(11,file='resolvent_exact_herm_m_1000_spec_range_22.dat')
!        open(111,file='resolvent_lanc_herm_full_m_1000_mlanc_100_.dat')

         ! non-symm, real
!        open(10,file='lanc_coeffs_non_symm_no_m_1000_mlanc_200_spec_range_22.dat')
!        open(14,file='fnorm_non_symm_no_m_1000_mlanc_200_spec_range_22.dat')
!        open(11,file='resolvent_exact_non_symm_m_1000_spec_range_22.dat')
!        open(111,file='resolvent_lanc_non_symm_no_m_1000_mlanc_200_spec_range_22.dat')

        ! non hermitian
!        open(10,file='lanc_coeffs_non_herm_no_m_1000_mlanc_200_spec_range_22.dat')
!        open(14,file='fnorm_non_herm_no_m_1000_mlanc_200_spec_range_22.dat')
!        open(11,file='resolvent_exact_non_herm_m_1000_spec_range_22.dat')
!        open(111,file='resolvent_lanc_non_herm_no_m_1000_mlanc_200_spec_range_22.dat')

!        open(10,file='lanc_coeffs_non_herm_full_m_1000_mlanc_200_spec_range_22.dat')
!        open(14,file='fnorm_non_herm_full_m_1000_mlanc_200_spec_range_22.dat')
!        open(11,file='resolvent_exact_non_herm_m_1000_spec_range_22.dat')
!        open(111,file='resolvent_lanc_non_herm_full_m_1000_mlanc_200_spec_range_22.dat')
!        open(222,file='orthinfo_non_herm_full_m_1000_mlanc_200_spec_range_22.dat')

!        open(10,file='lanc_coeffs_non_herm_minprbo_reps_m_1000_mlanc_200_spec_range_22.dat')
!        open(14,file='fnorm_non_herm_minprbo_reps_m_1000_mlanc_200_spec_range_22.dat')
!        open(11,file='resolvent_exact_non_herm_m_1000_spec_range_22.dat')
!        open(111,file='resolvent_lanc_non_herm_minprbo_reps_m_1000_mlanc_200_spec_range_22.dat')
!        open(222,file='orthinfo_non_herm_minprbo_reps_m_1000_mlanc_200_spec_range_22.dat')

!        open(10,file='lanc_coeffs_non_herm_prbo_1_reps_m_1000_mlanc_200_spec_range_22.dat')
!        open(14,file='fnorm_non_herm_prbo_1_reps_m_1000_mlanc_200_spec_range_22.dat')
!        open(11,file='resolvent_exact_non_herm_m_1000_spec_range_22.dat')
!        open(111,file='resolvent_lanc_non_herm_prbo_1_reps_m_1000_mlanc_200_spec_range_22.dat')
!        open(222,file='orthinfo_non_herm_prbo_1_reps_m_1000_mlanc_200_spec_range_22.dat')

        open(10,file='lanc_coeffs_non_herm_prbo_2_reps_m_1000_mlanc_300_spec_range_22.dat')
        open(14,file='fnorm_non_herm_prbo_2_reps_m_1000_mlanc_300_spec_range_22.dat')
        open(11,file='resolvent_exact_non_herm_m_1000_spec_range_22.dat')
        open(111,file='resolvent_lanc_non_herm_prbo_2_reps_m_1000_mlanc_300_spec_range_22.dat')
        open(222,file='orthinfo_non_herm_prbo_2_reps_m_1000_mlanc_300_spec_range_22.dat')
        rewind(10)
        rewind(14)
        rewind(11)
        rewind(111)
        rewind(222)

        ! file orthinfo saves info about reorthogonalization
        ! first number is Lanczos iteration, second is the number of rebiorthogs done for each basis set
        ! third term is the total, fourth is the cumulative freq of rebiorthogs, and last number
        ! the cumulative total number of rebiorthogs
        write(222,*) "   # iter,    oterm,    2*oterm,    orthfreq,     orthtot "

        if (herm==1) then
           print *
           write(6,*) "===== Doing Hermitian, case ====="
           print *

             
           ! diagonals
           do i=1,m

         ! Based on example from : "The Lanczos Alg with partial reorthogonalization"  H(i,i) = 1.0d0*i**2
               H(i,i) = del_spect*i**1
               Iden(i,i) = 1.0d0
          enddo

           ! Tridiagonalize the Hamiltonian
           call lanczos_herm(H,v,m,mlanc,alpha,beta,fnorm,Planc,full)


          ! Saving alpha, beta to array, for unnormalized Lanczos, beta -->beta squared. So easier to save sqrt(beta)
           write(6,*) "----- Printing alpha and beta -------"
           print *
           do i=1,mlanc
              write(6,*) ""
              write(10,'(1f16.10,2x,1f16.10)') alpha(i),beta(i)
              write(14,*) i,fnorm(i)
              write(6,*) 'iter = ',i,'alpha=',alpha(i),'beta(i)=',beta(i)
              write(6,*) ""
           enddo
           write(10,'(1f16.10,2x,1f16.10)') 0.00,beta(mlanc+1)
           close(10)
           close(14)
           close(15)
        elseif (herm==2) then
           print *
           write(6,*) "===== Doing non-symmetric, real case ====="
           print *


           ! diagonal
           do i=1,m
!              Hmat(i,i)= md*mod(i,mi)

          ! Based on example from : "The Lanczos Alg with partial reorthogonalization"  H(i,i) = 1.0d0*i**2
               Hmat(i,i) = del_spect*i**1
               Iden(i,i) = 1.0d0
          enddo

           ! Tridiagonalize the Hamiltonian
           call lanczos_non_symm(Hmat,v0,m,mlanc,alpha,beta,delta)

        ! Saving alpha, beta to array
        write(6,*) "----- Printing alpha and beta -------"
           do i=1,mlanc
             write(6,*) ""
!          write(10,'(3f16.10)') alpha(i),beta(i),delta(i)
             write(10,'(3e16.8)') alpha(i),beta(i),delta(i)
             write(6,*)'alpha=',alpha(i),'beta(i)=',beta(i), &
                &  'delta(i)=',delta(i)
             write(6,*) ""
           write(10,'(3e16.8)') 0.00,beta(mlanc+1),delta(mlanc+1)
           enddo
           close(10)
        else 
           print *
           write(6,*) "===== Doing non-Hermitian, case ====="
           print *

           ! diagonal
           do i=1,m
               H(i,i)=del_spect*i**1
!               H(i,i) = hterm
               Iden(i,i) = 1.0d0
          enddo


           ! Tridiagonalize the Hamiltonian
           if (reorth) then
              if (full) then
                 write(6,*) " "
                 write(6,*) " Doing full rebiorthogonalization "
                 call lanczos_non_herm(H,v,m,mlanc,zalpha,zbeta,zdelta,& 
                                   zborth,fnorm,paige,reorth,mult_type)
              else

                 select case (prbo_type)
                     case ('prbo 1')
                         write(6,*) " "
                         write(6,*) " Doing partial rebiorthogonalization of type ",prbo_type
                         prbo = 1
                         call lanczos_non_herm_prbo(H,v,m,mlanc,zalpha,zbeta,zdelta,& 
                                   zborth,fnorm,paige,reorth,eps1,prbo)

                     case ('prbo 2')
                         write(6,*) " "
                         write(6,*) " Doing partial rebiorthogonalization of type ",prbo_type
                         prbo = 2
                         call lanczos_non_herm_prbo(H,v,m,mlanc,zalpha,zbeta,zdelta,& 
                                   zborth,fnorm,paige,reorth,eps1,prbo)
                     case ('minprbo')
                         write(6,*) " "
                         write(6,*) " Doing partial rebiorthogonalization of type ",prbo_type
                         call lanczos_non_herm_minprbo(H,v,m,mlanc,zalpha,zbeta,&
                                  zdelta,zborth,fnorm,paige,reorth,eps1,mult_type)

                     case default
                         call lanczos_non_herm_minprbo(H,v,m,mlanc,zalpha,zbeta,&
                                  zdelta,zborth,fnorm,paige,reorth,eps1,mult_type)
                 end select 
              endif
           else
              write(6,*) " "
              write(6,*) "Not doing any re biorthogonalization"
              call lanczos_non_herm(H,v,m,mlanc,zalpha,zbeta,zdelta,& 
                                   zborth,fnorm,paige,reorth,mult_type)
           endif

           !  Saving alpha, beta to array
!           write(6,*) "----- Printing alpha, beta and gamma -------"
           do i=1,mlanc
!             write(6,*) ""
!             write(6,*) "lanc iteration =",i
!             write(6,*) ""
             write(10,'(6e16.8)') zalpha(i),zbeta(i),zdelta(i)
             write(14,'(1i3,2e16.8)') i,fnorm(i),abs(aimag(zborth(i)))
!             write(6,*) ""
!             write(6,*)'alpha=',zalpha(i)
!             write(6,*)""
!             write(6,*)'beta(i)=',zbeta(i)
!             write(6,*)""
!             write(6,*)'delta(i)=',zgamma(i)
!             write(6,*) ""
           enddo
           write(10,'(6e16.8)') 0.0d0,0.0d0,zbeta(mlanc+1),zdelta(mlanc+1)
           close(14)
           close(10)
        
      endif  ! herm type

!        deallocate(fnorm,estnorm)

!! Exact method of calculating resolvent matrix elements      
!! by inverting the matrix: (w*I + i*eta - H)

!! Resolvent: <v_0|(w*I + i*eta - H)^-1|v_0>
!! Dielectric function: e(w) = 1.0 - <v_0|(w*I + i*eta - H)^-1|v_0>

      ! allocate array to save matrix elem of resolvent
      lda = m
      lwork = 64*m
      allocate(work(lwork),work1(lwork),ipiv(m),ipiv1(m))

      allocate(resolv_exact(0:omega_max,3),resolv_lanc(0:omega_max,3))
      resolv_exact = 0.0d0
      resolv_lanc = 0.0d0

      ! Hamiltonian, taking care of real matrices
      if ((herm .eq. 0) .or. (herm .eq. 2)) then
             H = Hmat
             v = v0
      endif 

   !  Calculating norm of the initial vector , v_0
      rterm = dot_product(v,v)
      print *
      print *," <v_0|v_0> =",rterm, "|| v_0 || = ",sqrt(rterm)
      print *


      print *
      print *, "------ Exact method of calculating resolvent ----"
      print *
      do i=0,omega_max
        omega = omega_scal*i
        print *
        print *, " i = ",i,"omega = ",omega
        print *

        ! Positive and negative energies
        temp_mat = (omega+ii*eta)*Iden - H
        temp_mat1 = (-omega+ii*eta)*Iden - H

        ! Do LU factorization of A, A = L*U
        call zgetrf(m,m,temp_mat,lda,ipiv,info)
        call zgetrf(m,m,temp_mat1,lda,ipiv1,info)

        ! Compute inverse of A
        call zgetri(m,temp_mat,lda,ipiv,work,lwork,info)
        call zgetri(m,temp_mat1,lda,ipiv1,work1,lwork,info)

        ! Compute (w*I + i*eta - H)^-1|v_0>
        t = 0.0d0
        u = 0.0d0
        t = matmul(temp_mat,v)
        u = matmul(temp_mat1,v)


        ! Calculate resolvent matrix element
        !  <v_0|(w*I + i*eta - H)^-1|v_0>

        zterm = dot_product(v,t)
        zterm1 = dot_product(v,u)
        resolv_exact(i,1) = zterm
        resolv_exact(i,2) = zterm1
        resolv_exact(i,3) = zterm+zterm1
        print *
        print *, "omega = ",omega," eta = ",eta
        print *,"<v_0|(w*I + i*eta - H)^-1|v_0> = ",zterm
!        print *
!        print *,"<v_0|(-w*I + i*eta - H)^-1|v_0> = ",zterm1
!        print *
!        print *,"R(w) + R(-w) =",zterm+zterm1
        print*
      enddo ! omega loop 
      deallocate(work,work1,ipiv,ipiv1)

!! Approximate method of calculating diagonal resolvent matrix elements      
! by continued fraction from Lanczos coeffs: 
! <v_0|(w*I + i*eta - H)^-1|v_0> ~ <v_0|v_0>*1/(w+i*eta-a_1 - b^2_2/(w+i*eta-a_2 ...))

      print *
      print *, "------ Approximate method of calculating diagonal element of resolvent: Lanczos alg ----"
      print *

      ! Taking care of non Hermitian case and non symmetric case where we have (z)alpha, (z)beta and (z)delta
      if ((herm .eq. 0) .or. (herm .eq. 1) ) then
             delta = beta
      endif 

      do i=0,omega_max
        omega = omega_scal*i
        print *
        print *, " i = ",i,"omega = ",omega," eta =",eta
        print *

        ! Positive and negative energies
        zterm = omega + ii*eta
        zterm1 = -omega + ii*eta

        if (herm .ne. 3) then

                ! Get asymptote of continued frac
                arg = (omega - alpha(mlanc))**2 - 4.0d0*beta(mlanc+1)*delta(mlanc+1)
                arg = sqrt(arg)
                rp = 0.5d0*(omega - alpha(mlanc) + arg)
                rm = 0.5d0*(omega - alpha(mlanc) - arg)

                if (aimag(rp) .lt. 0.0d0) then
                   rlim = rp
                else
                   rlim = rm
                endif

                al = zterm - alpha(mlanc) - rlim
                be = zterm1 - alpha(mlanc) - rlim
                print *
                print *,"omega = ",omega, "asymptotes: al, be ",al,be
                print *

                ! Adding up the terms in continued frac

                do j=mlanc,1,-1
                  al = zterm - alpha(j) - beta(j+1)*delta(j+1)/al 
                  be = zterm1 - alpha(j) - beta(j+1)*delta(j+1)/be 
                enddo
         else
                ! Get asymptote of continued frac
                arg = (omega - zalpha(mlanc))**2 - 4.0d0*zbeta(mlanc+1)*zdelta(mlanc+1)
!                print *
!                print *," zterm = ",zterm,"alpha =",zalpha(mlanc),"beta =",zbeta(mlanc+1),"zdelta =",zdelta(mlanc+1)
!                print *
                arg = sqrt(arg)
                rp = 0.5d0*(omega - zalpha(mlanc) + arg)
                rm = 0.5d0*(omega - zalpha(mlanc) - arg)

                if (aimag(rp) .lt. 0.0d0) then
                   rlim = rp
                else
                   rlim = rm
                endif

                al = zterm - zalpha(mlanc) - rlim
                be = zterm1 - zalpha(mlanc) - rlim
                print *
                print *, "asymptotes: al, be ",al,be
                print *

                ! Adding up the terms in continued frac

                do j=mlanc,1,-1
                  al = zterm - zalpha(j) - zbeta(j+1)*zdelta(j+1)/al 
                  be = zterm1 - zalpha(j) - zbeta(j+1)*zdelta(j+1)/be 
                enddo
        endif

        ! add this resolvent to array
        zterm3 = dot_product(v,v)

        resolv_lanc(i,1) = zterm3/al !+ zterm1/be
        resolv_lanc(i,2) = zterm3/be !+ zterm1/be
        resolv_lanc(i,3) = zterm3/al + zterm3/be

        print *
        print *,"mlanc=",mlanc, "omega = ",omega
        print*
        print *,"<v_0|(w*I + i*eta - H)^-1|v_0> ~ ",zterm3/al
        print *
!        print *,"<v_0|(w*I - i*eta + H)^-1|v_0> ~ ",zterm3/be
!        print *
!        print *,"R(w) + R(-w) =",zterm3/al + zterm3/be
!        print *
      enddo  

      !  Saving omega, resolvent, dielectric fxn calculations to file 
      do i=0,omega_max
        zterm = resolv_exact(i,1)
        zterm1 = resolv_exact(i,2)
        zterm2 = resolv_exact(i,3)
        zterm3 = resolv_lanc(i,1)
        zterm4 = resolv_lanc(i,2)
        zterm5 = resolv_lanc(i,3)
        omega = omega_scal*i

        ! Save just -1*resolvent
!        write(11,'(7e16.8)') omega,-real(zterm),-aimag(zterm),-real(zterm1),-aimag(zterm1),-real(zterm2),-aimag(zterm2)
!        write(111,'(7e16.8)') omega,-real(zterm3),-aimag(zterm3),-real(zterm4),-aimag(zterm4),-real(zterm5),-aimag(zterm5)
        
        ! Save dielectric function
        write(11,'(7e16.8)') omega,real(1.0d0-zterm),aimag(1.0d0-zterm),real(1.0d0-zterm1),aimag(1.0d0-zterm1), &
                             real(1.0d0-zterm2),aimag(1.0d0-zterm2)
        write(111,'(7e16.8)') omega,real(1.0d0-zterm3),aimag(1.0d0-zterm3),real(1.0d0-zterm4),aimag(1.0d0-zterm4), &
                             real(1.0d0-zterm5),aimag(1.0d0-zterm5)
      enddo
      close(11)
      close(111)
!      close(222)

      deallocate(fnorm,estnorm)
        stop
        end program lanczos_test


!======================================================================
!                      Functions and subroutines
!
!======================================================================
      
        subroutine check_orthog_frob(A,m,k,frob)
!  Calculates the loss of orthog of a real matrix with Frobenius norm manually.
        implicit none
        integer :: i,j,k,m,f
        real(kind=8) :: frob,A(m,k),fterm
        real(kind=8),allocatable :: iden(:,:) 

        allocate(iden(k,k))
        iden = 0.0
!          write(6,*) "shape Pj=",shape(Pj)
!          write(6,*) " "  
!  Set the identity matrix
        do f = 1,k
          iden(f,f) = 1.0d0
        enddo    
        iden = iden-matmul(transpose(A),A)

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
        write(6,*) " " 
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "frob=",frob,"alt",norm2(real(iden))
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "
!          fnorm(k) = fterm 
        deallocate(iden)
        return
        end subroutine check_orthog_frob


        subroutine zcheck_orthog_frob(A,B,m,k,frob)
!  Calculates the loss of orthog of a real matrix with Frobenius norm manually.
        implicit none
        integer :: i,j,k,m,f
        real(kind=8) :: frob,fterm
        complex(kind=8) :: A(m,k),B(m,k) 
        complex(kind=8),allocatable :: iden(:,:) 

        allocate(iden(k,k))
        iden = 0.0
!          write(6,*) "shape Pj=",shape(Pj)
!          write(6,*) " "  
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

        subroutine lanc_coeff_to_matrix(alpha,beta,delta,mlanc,T)
!  Takes in the real Lanczos coefficients, a_i, b_i, c_i, and constructs explicitly the
! tridiagonal complex, square matrix T, which has these coefficients
! This assumes that in beta, delta, beta(1)=delta(1)=0.0

        implicit none
        integer :: i,j,k,m,n,mlanc
        real(kind=8) :: alpha(mlanc),beta(mlanc+1),delta(mlanc+1),rterm1,rterm2,rterm3
        complex(kind=8) :: T(mlanc,mlanc) 
        complex(kind=8),allocatable :: iden(:,:) 

        ! Initialize T
         
        do j=1,mlanc
          do i=1,j
             if (abs(i-j) .eq. 1) then
                rterm1=beta(j)
                rterm2=delta(j)
                T(i,j)=rterm1
                T(j,i)=rterm2
!                   write(6,*) ""
!                   write(6,*) "---- i=",i,"j=",j,"T(i,j) =",T(i,j)
             endif     
          enddo
        enddo

!! diagonals
        do i=1,mlanc
           T(i,i)= alpha(i)
        enddo

        return
        end subroutine lanc_coeff_to_matrix

        subroutine zlanc_coeff_to_matrix(zalpha,zbeta,zdelta,mlanc,T)
!  Takes in the Lanczos coefficients, a_i, b_i, c_i, stored in complex arrays and constructs explicitly the
! tridiagonal complex, square matrix T, which has these coefficients
! This assumes that in beta, delta, beta(1)=delta(1)=0.0

        implicit none
        integer :: i,j,k,m,n,mlanc
        complex(kind=8) :: zalpha(mlanc),zbeta(mlanc+1),zdelta(mlanc+1),zterm1,zterm2,zterm3
        complex(kind=8) :: T(mlanc,mlanc) 
        complex(kind=8),allocatable :: iden(:,:) 

        ! Initialize T
         
        do j=1,mlanc
          do i=1,j
             if (abs(i-j) .eq. 1) then
                zterm1=zbeta(j)
                zterm2=zdelta(j)
                T(i,j)=zterm1
                T(j,i)=zterm2
!                   write(6,*) ""
!                   write(6,*) "---- i=",i,"j=",j,"T(i,j) =",T(i,j)
             endif     
          enddo
        enddo

!! diagonals
        do i=1,mlanc
           T(i,i)= zalpha(i)
        enddo

        return
        end subroutine zlanc_coeff_to_matrix



        subroutine gram_schmidt_class(A,m,n,B,fnorm)
! Classical Gram-Schmidt from paper: "The loss of orthogonality in the
! Gram-Schmidt orthogonalization process"
! *** CANT MAKE A AND B THE SAME!!!!
! Takes in the matrix A, and reorthogonalizes it to matrix B.
        implicit none
        integer :: i,j,k,l,m,f,n
        real(kind=8) :: b0,temp1,tol,fnorm(n),fterm,lterm,term
        real(kind=8) :: u0,lnorm,DLANGE
        real(kind=8) :: rand_num(m)
        real,parameter :: pi=3.141592654d0
        complex,parameter :: ii=complex(0.d0,1.d0)
        
        real(kind=8):: vm1(m),a_j(m),s(m),q_j(m),v_j(m),w_j(m),q_k(m)
        
        real(kind=8) :: A(m,n)

!  Matrices for calculating biorthogonality
        real(kind=8) :: B(m,n) ! since doing orthogonalization
!       from 1 to j-1
        real(kind=8), allocatable :: iden(:,:),Pj(:,:)

        vm1=0.0
        a_j=0.0
        q_k=0.0
        q_j=0.0
        v_j=0.0
        w_j = 0.0

! Initial vector        
        a_j=A(:,1)
        u0=dot_product(a_j,a_j)
        u0=sqrt(u0)
!        write(6,*) "--- The Euclidean norm =",u0
        q_j=a_j/u0
        B(:,1) = q_j


!----------  Starting the Gram-Schmidt iteration ------------------------
        do k=1,n
          write(6,*) " "
          write(6,*) "******** Gram-Schmidt iteration =",k,"***********"

! Traditional formulation : "Loss of orthogonality in the Gram-Schmidt
! Orthogonalization process"

!  1st orthogonlization
          a_j = A(:,k)
          v_j = a_j

          do j=1,k-1
             q_k = B(:,j) ! Needs to be the orthogonal matrix
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

!!  New norm
            lnorm = dot_product(v_j,v_j)
            lnorm = sqrt(lnorm)
!            write(6,*) "norm v_j=",lnorm
!            write(6,*) " "  
        
            do i=1,m
              q_j(i)=(1/lnorm)*v_j(i)
            enddo

!!  Add new vectors 
            B(:,k) = q_j

! Alternate formulation of Gram-Schmidt , "Rounding error analysis of
! the classical Gram-Schmidt orthogonalization process"

!            allocate(iden(m,m),Pj(m,k-1))
!            Pj = Planc(:,1:k-1)
!!            write(6,*) "shape Pj=",shape(Pj)
!!            write(6,*) " "  
!!  Set the identity matrix
!            do f = 1,m
!              iden(f,f) = 1.0d0
!            enddo    
!            iden = iden-matmul(Pj,transpose(Pj))
!
!! u_j
!            v_j = matmul(iden,a_j)             
!
!!  New norm
!            lnorm = dot_product(v_j,v_j)
!            lnorm = sqrt(lnorm)
!            write(6,*) "norm v_j=",lnorm
!            write(6,*) " "  
!        
!            do i=1,m
!              q_j(i)=(1/lnorm)*q_j(i)
!            enddo
!
!!  Add new vectors 
!            Planc(:,k) = q_j
!            deallocate(iden,Pj)

!  Checking biorthogonality with the frobenius norm and spectral norm
            allocate(iden(k,k),Pj(m,k))
            iden = 0.0
            Pj = 0.0
            Pj = B(:,1:k)
!            write(6,*) "shape Pj=",shape(Pj)
!            write(6,*) " "  
!  Set the identity matrix
            do f = 1,k
              iden(f,f) = 1.0d0
            enddo    
!            write(6,*) " I"
!            write(6,*) iden
!            write(6,*) " "
            iden = iden-matmul(transpose(Pj),Pj)
!            write(6,*) " I - Q^T*Q has shape ",shape(iden)
!            write(6,*) " "
!            write(6,*) iden
!            write(6,*) " "
!            write(6,*) " Q^T*Q"
!            write(6,*) " "
!            write(6,*) matmul(transpose(Pj),Pj)
!            write(6,*) " "
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
            write(6,*) "frob=",fterm !norm2(real(iden,kind=8)) !,"L2 nrm=",term
            write(6,*) "^^^^^^^^^^^^^^^^^^"
            write(6,*) " "
            fnorm(k) = fterm 

! deallocate necessary arrays
            deallocate(iden,Pj)

        enddo ! END OF GRAM-SCHMIDT LOOP
        return
        end subroutine gram_schmidt_class
        

        subroutine gram_schmidt_class_2(A,m,n,B,frob)
! Classical Gram-Schmidt from paper: "The loss of orthogonality in the
! Gram-Schmidt orthogonalization process"
! *** CANT MAKE A AND B THE SAME!!!!
! Just returns the orthog matrix B and the frobenius norm after the GS
! has completed
        implicit none
        integer :: i,j,k,l,m,f,n
        real(kind=8) :: b0,temp1,tol,fnorm(n),fterm,lterm,term,frob
        real(kind=8) :: u0,lnorm,DLANGE
        real(kind=8) :: rand_num(m)
        real,parameter :: pi=3.141592654d0
        complex,parameter :: ii=complex(0.d0,1.d0)
        
        real(kind=8):: vm1(m),a_j(m),s(m),q_j(m),v_j(m),w_j(m),q_k(m)
        
        real(kind=8) :: A(m,n)

!  Matrices for calculating biorthogonality
        real(kind=8) :: B(m,n) ! since doing orthogonalization
!       from 1 to j-1
        real(kind=8), allocatable :: iden(:,:),Pj(:,:)


!   Initializing arrays to 0
        vm1=0.0
        a_j=0.0
        q_k=0.0
        q_j=0.0
        v_j=0.0
        w_j = 0.0
        B = 0.0

! Initial vector        
        a_j=A(:,1)
        u0=dot_product(a_j,a_j)
        u0=sqrt(u0)
!        write(6,*) "--- The Euclidean norm =",u0
        q_j=a_j/u0
        B(:,1) = q_j


!----------  Starting the Gram-Schmidt iteration ------------------------
        do k=1,n
!          write(6,*) " "
!          write(6,*) "******** Gram-Schmidt iteration =",k,"***********"

! Traditional formulation : "Loss of orthogonality in the Gram-Schmidt
! Orthogonalization process"

!  1st orthogonlization
          a_j = A(:,k)
          v_j = a_j

          do j=1,k-1
             q_k = B(:,j) ! Needs to be the orthogonal matrix
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

!!  New norm
            lnorm = dot_product(v_j,v_j)
            lnorm = sqrt(lnorm)
!            write(6,*) "norm v_j=",lnorm
!            write(6,*) " "  
        
            do i=1,m
              q_j(i)=(1/lnorm)*v_j(i)
            enddo

!!  Add new vectors 
            B(:,k) = q_j


!  Checking biorthogonality with the frobenius norm and spectral norm
            allocate(iden(k,k),Pj(m,k))
            iden = 0.0
            Pj = 0.0
            Pj = B(:,1:k)
!            write(6,*) "shape Pj=",shape(Pj)
!            write(6,*) " "  
!  Set the identity matrix
            do f = 1,k
              iden(f,f) = 1.0d0
            enddo    
!            write(6,*) " I"
!            write(6,*) iden
!            write(6,*) " "
            iden = iden-matmul(transpose(Pj),Pj)
!            write(6,*) " I - Q^T*Q has shape ",shape(iden)
!            write(6,*) " "
!            write(6,*) iden
!            write(6,*) " "
!            write(6,*) " Q^T*Q"
!            write(6,*) " "
!            write(6,*) matmul(transpose(Pj),Pj)
!            write(6,*) " "
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
            frob = fterm
            

! deallocate necessary arrays
            deallocate(iden,Pj)

        enddo ! END OF GRAM-SCHMIDT LOOP
        return
        end subroutine gram_schmidt_class_2


        subroutine gram_schmidt_class_sing(A,m,n,p,q,frob)
! Classical Gram-Schmidt from paper: "The loss of orthogonality in the
! Gram-Schmidt orthogonalization process"
! *** CANT MAKE A AND B THE SAME!!!!
! Just returns the vector p thats orthogonalized against all the vectors
! in A, to get q,  and the frobenius norm after the vector q has been
! added to the basis.
        implicit none
        integer :: i,j,k,l,m,f,n
        real(kind=8) :: b0,temp1,tol,fnorm(n),fterm,lterm,term,frob
        real(kind=8) :: u0,lnorm,DLANGE
        real(kind=8) :: rand_num(m)
        real,parameter :: pi=3.141592654d0
        complex,parameter :: ii=complex(0.d0,1.d0)
        
        real(kind=8):: vm1(m),a_j(m),s(m),q_j(m),v_j(m),w_j(m),q_k(m)
        
        real(kind=8) :: A(m,n)

!  Matrices for calculating biorthogonality
        real(kind=8) :: p(m),q(m) ! since doing orthogonalization
!       from 1 to j-1
        real(kind=8), allocatable :: iden(:,:),Pj(:,:)


        vm1=0.0
        a_j=0.0
        q_k=0.0
        q_j=0.0
        v_j=0.0
        w_j = 0.0
        q = 0.0

!  Alternative way to do Gram-Schmidt with (I-Q_k*(P_k)^H)
        allocate(iden(m,m))
        iden = 0.0
            
!  Set the identity matrix
        do f = 1,m
           iden(f,f) = 1.0d0
        enddo  

! For Hermitian case, works well  
        iden = iden-matmul(A,transpose(A))

! re biorthogonalizing
        q = matmul(iden,p)
        deallocate(iden)

!  Checking biorthogonality with the frobenius norm and spectral norm
        allocate(iden(n+1,n+1),Pj(m,n+1))
        iden = 0.0
        Pj = 0.0
        Pj(:,1:n) = A
        Pj(:,n+1) = q
!            write(6,*) "shape Pj=",shape(Pj)
!            write(6,*) " "  

!  Set the identity matrix
        do f = 1,n+1
           iden(f,f) = 1.0d0
        enddo    
!            write(6,*) " I"
!            write(6,*) iden
!            write(6,*) " "
        iden = iden-matmul(transpose(Pj),Pj)
!            write(6,*) " I - Q^T*Q has shape ",shape(iden)
!            write(6,*) " "
!            write(6,*) iden
!            write(6,*) " "
!            write(6,*) " Q^T*Q"
!            write(6,*) " "
!            write(6,*) matmul(transpose(Pj),Pj)
!            write(6,*) " "

!  Start calc of frobenius norm    
        fterm = 0.0
        do f=1,n+1
           do j=1,n+1
              fterm = fterm + (abs(iden(f,j)))**2
!                 write(6,*) "i=",f,"j=",j,"fterm_ij=",fterm
!                 write(6,*) " "
           enddo
        enddo 
        fterm=sqrt(fterm) 
        frob = fterm
            
        write(6,*) " " 
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "n = ",k,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

! deallocate necessary arrays
        deallocate(iden,Pj)

        return
        end subroutine gram_schmidt_class_sing


        subroutine gram_schmidt_class_sing_pro(A,m,n,o,p,q,indx,frob)
! Classical Gram-Schmidt from paper: "The loss of orthogonality in the
! Gram-Schmidt orthogonalization process"
! *** CANT MAKE A AND B THE SAME!!!!
! Just returns the vector p thats orthogonalized against all the vectors
! in A that p has lost orthog to, denoted,with index indx(n) to get q,  and the frobenius norm after the vector q has been
! added to the lanczosbasis.
        implicit none
        integer :: i,j,k,l,m,f,n,o,indx(o),iterm
        real(kind=8) :: b0,temp1,tol,fnorm(n),fterm,lterm,term,frob
        real(kind=8) :: u0,lnorm,DLANGE
        real(kind=8) :: rand_num(m)
        real,parameter :: pi=3.141592654d0
        complex,parameter :: ii=complex(0.d0,1.d0)
        
!        real(kind=8) :: alpha(mlanc),beta(mlanc)
        real(kind=8):: vm1(m),a_j(m),s(m),q_j(m),v_j(m),w_j(m),q_k(m)
        
        real(kind=8) :: A(m,n)

!  Matrices for calculating biorthogonality
        real(kind=8) :: p(m),q(m) ! since doing orthogonalization
!       from 1 to j-1
        real(kind=8), allocatable :: iden(:,:),Pj(:,:)


        vm1=0.0
        a_j=0.0
        q_k=0.0
        q_j=0.0
        v_j=0.0
        w_j = 0.0
        q = 0.0

!  Alternative way to do Gram-Schmidt with (I-Q_k*(P_k)^H)
        allocate(iden(m,m))
        iden = 0.0
            
!  Set the identity matrix
        do f = 1,m
           iden(f,f) = 1.0d0
        enddo  

! Form the matrix that holds the selected vectors
        allocate(Pj(m,o))
        Pj = 0.0
        do i=1,o
          iterm = indx(i)
          Pj(:,i) = A(:,iterm)
        enddo

! For Hermitian case, works well  
        iden = iden-matmul(Pj,transpose(Pj))

! re biorthogonalizing
        q = matmul(iden,p)
        deallocate(iden,Pj)

!  Checking biorthogonality with the frobenius norm and spectral norm
        allocate(iden(n+1,n+1),Pj(m,n+1))
        iden = 0.0
        Pj = 0.0
        Pj(:,1:n) = A(:,1:n)
        Pj(:,n+1) = q
!            write(6,*) "shape Pj=",shape(Pj)
!            write(6,*) " "  

!  Set the identity matrix
        do f = 1,n+1
           iden(f,f) = 1.0d0
        enddo    
!            write(6,*) " I"
!            write(6,*) iden
!            write(6,*) " "
        iden = iden-matmul(transpose(Pj),Pj)
!            write(6,*) " I - Q^T*Q has shape ",shape(iden)
!            write(6,*) " "
!            write(6,*) iden
!            write(6,*) " "
!            write(6,*) " Q^T*Q"
!            write(6,*) " "
!            write(6,*) matmul(transpose(Pj),Pj)
!            write(6,*) " "

!  Start calc of frobenius norm    
        fterm = 0.0
        do f=1,n+1
           do j=1,n+1
              fterm = fterm + (abs(iden(f,j)))**2
!                 write(6,*) "i=",f,"j=",j,"fterm_ij=",fterm
!                 write(6,*) " "
           enddo
        enddo 
        fterm=sqrt(fterm) 
        frob = fterm
            
        write(6,*) " " 
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "n = ",k,"frob=",frob !norm2(real(iden,kind=8)) !,"L2 nrm=",term
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

! deallocate necessary arrays
        deallocate(iden,Pj)

        return
        end subroutine gram_schmidt_class_sing_pro


        subroutine zgs_class_two_side(A,B,m,n,C,D)
! Two sided full Gram-Schmidt for complex matrices
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
    
!------- Starting the two sided Gram-Schmidt  -------
        do g=1,n
           write(6,*) " "
           write(6,*) "********* Gram-Schmidt iteration=",g,"*********"
           write(6,*) " "
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
              Pj = 0.0
              Qj = 0.0

! For Hermitian case, works well  
              iden1 = iden1-matmul(Qj,conjg(transpose(Pj)))
              iden2 = iden2-matmul(Pj,conjg(transpose(Qj)))

! re biorthogonalizing
              q_j = matmul(iden1,b_j)
              p_j = matmul(iden2,a_j)
            else
              q_j = b_j
              p_j = a_j
            endif

            deallocate(iden1,iden2,Pj,Qj)
          
!!  New norm
            zterm = dot_product(q_j,p_j) ! following netlib Alg 7.13
            lnorm1 = sqrt(abs(zterm))  ! following yambo paper
            lnorm2 = conjg(zterm)/lnorm1
            write(6,*) "lnorm2 =",lnorm2
            write(6,*) ""

! Normalizing and updating Lanczos vectors
            q_j = (1.0/lnorm1)*q_j
            p_j = (1.0/conjg(lnorm2))*p_j


            C(:,g) = p_j
            D(:,g) = q_j
        enddo ! loop over Lanczos columns for two-sided CGS

! Checking the level of orthog with Frobenius norm after two-sided GS
        call zcheck_orthog_frob(C,D,m,g,frob)
        write(6,*) " " 
        write(6,*) "^^^^^^^^^^^^^^^^^^"    
        write(6,*) "n = ",g,"frob=",frob 
        write(6,*) "^^^^^^^^^^^^^^^^^^"
        write(6,*) " "

        return
        end subroutine zgs_class_two_side


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


        
        write(6,*) "||A||_F =",norm2(A)
        write(6,*) " "
        write(6,*) " Doing Gram-Schmidt with reorthogonalization" 

!   Initializing arrays to 0
        vm1=0.0
        a_j=0.0
        q_k=0.0
        q_j=0.0
        v_j=0.0
        w_j=0.0
        t_j=0.0
        Planc = 0.0

       
! Normalizing initial vector 
        a_j=A(:,1)
        u0=dot_product(a_j,a_j)
        u0=sqrt(u0)
        write(6,*) "--- The Euclidean norm =",u0
        q_j=a_j/u0
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
!            write(6,*) "norm v_j=",lnorm
!            write(6,*) " "  

!  2nd orthogonalization
          t_j = v_j
          w_j = v_j
          do j=1,k-1
             q_k = Planc(:,j)
             term = dot_product(v_j,q_k)
!          write(6,*) ""
!          write(6,*) "j =",j,"alpha(k,j) =",term
!          write(6,*) ""
             do i=1,m
!            write(6,*) ""
!            write(6,*) " --- i = ",i,"------"
!            write(6,*) "before updating, t(i)=",t(i)
               w_j(i)=w_j(i)-term*q_k(i) 
!            write(6,*) ""
             enddo
           enddo
!!  New norm
            lnorm = dot_product(w_j,w_j)
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
 
        subroutine lanczos_herm_re(A,v0,m,mlanc,alpha,beta,fnorm)
        implicit none

! Carries out the Lanczos algorithm for a real matrix
! Starting with input vector v0

        integer :: i,j,k,l,m,n,mlanc,f
        real(kind=8) :: b0,temp1,tol,fnorm(mlanc),fterm,lterm
        real(kind=8) :: u0,lnorm,DLANGE
        real(kind=8) :: rand_num(m)
        real,parameter :: pi=3.141592654d0
        complex,parameter :: ii=complex(0.d0,1.d0)
        
        real(kind=8) :: alpha(mlanc),beta(mlanc),c(mlanc),work(mlanc)
        real(kind=8):: v0(m),vm1(m),v(m),s(m),t(m),u(m),q_j(m),q_k(m),q_jm1(m)
        
        real(kind=8) :: A(m,m)

!  Matrices for calculating biorthogonality
        real(kind=8) :: Planc(m,mlanc)
        real(kind=8), allocatable :: iden(:,:),Pj(:,:)
        integer,parameter :: sedsz=33
        integer :: seed(sedsz) 

        
        write(6,*) "||A||_F =",norm2(A)

!   Initializing arrays to 0
        alpha=0.0
        beta=0.0
        c=0.0
        q_jm1=0.0
        q_j=0.0
        t=0.0
        u=0.0

! Normalizing initial vector v0
        u0=dot_product(v0,v0)
        u0=sqrt(u0)
        write(6,*) "--- The Euclidean norm =",u0
        q_j=v0/u0
        beta(1)=0

!  Collecting q_1 and p_1 to add to P and Q
        Planc(:,1) = q_j


!----------  Starting the Lanczos iteration ------------------------
        do k=1,mlanc
          write(6,*) " "
          write(6,*) "******** Lanczos iteration =",k,"***********"

!  Checking orthogonality with the frobenius norm
          allocate(iden(k,k),Pj(m,k))
          iden = 0.0
          Pj = 0.0
          Pj = Planc(:,1:k)
          write(6,*) "shape Pj=",shape(Pj)
          write(6,*) " "  
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
          write(6,*) "frob=",fterm,"alt",norm2(real(iden))
          write(6,*) "^^^^^^^^^^^^^^^^^^"
          write(6,*) " "
          fnorm(k) = fterm 
          deallocate(iden,Pj)

          t = matmul(A,q_j) !doing A*v_j

          ! u_{j+1}
          u = t - beta(k)*q_jm1 
          
          alpha(k)=dot_product(q_j,u)
          write(6,*) ""
          write(6,*) "alpha(k) =",alpha(k)
          write(6,*) ""

          ! add second term to u_{j+1}
          u = u - alpha(k)*q_j 
          
!  New norm
          lnorm = dot_product(u,u)
          lnorm = sqrt(lnorm)
          if (k .ne. mlanc) then
          beta(k+1)=lnorm
!            write(6,*) ""
!            write(6,*) "for iteration ",k,"beta =",lnorm
!            write(6,*) ""
          endif

!  Updating v_(j-1) for the next run (j+1) 
          q_jm1 = q_j

!  Updating v_(j+1), to be used at start of (j+1) iteration, can use the tol here also
          if (k .ne. mlanc) then
             q_j = (1/lnorm)*u
         
!  Add new vectors 
            Planc(:,k+1) = q_j
          else
             
          endif

        enddo ! END OF LANCZOS LOOP

        return
        end subroutine lanczos_herm_re 

        subroutine lanczos_herm_re_full(A,v0,m,mlanc,alpha,beta,fnorm)
        implicit none
! Lanczos algorithm for real matrix, carrying out        
!  Full reorthogonalization
        integer :: i,j,k,l,m,n,mlanc,f
        real(kind=8) :: b0,temp1,tol,fnorm(mlanc),fterm,lterm
        real(kind=8) :: u0,lnorm,DLANGE
        real(kind=8) :: rand_num(m)
        real,parameter :: pi=3.141592654d0
        complex,parameter :: ii=complex(0.d0,1.d0)
        
        real(kind=8) :: alpha(mlanc),beta(mlanc),c(mlanc),work(mlanc)
        real(kind=8):: v0(m),vm1(m),v(m),s(m),t(m),u(m),q_j(m),q_k(m),q_jm1(m)
        
        real(kind=8) :: A(m,m)

!  Matrices for calculating reorthogonality
        real(kind=8) :: Planc(m,mlanc)
        real(kind=8), allocatable :: iden(:,:),Pj(:,:)
        integer,parameter :: sedsz=33
        integer :: seed(sedsz) 

        
        write(6,*) "||A||_F =",norm2(A)

!   Initializing arrays to 0
        alpha=0.0
        beta=0.0
        c=0.0
        q_jm1=0.0
        q_j=0.0
        t=0.0
        u=0.0

!  Normalizing the initial vector , v_0
        u0=dot_product(v0,v0)
        u0=sqrt(u0)
        write(6,*) "--- The Euclidean norm =",u0
        q_j=v0/u0
        beta(1)=0

!  Collecting q_1 and p_1 to add to P and Q
        Planc(:,1) = q_j


!----------  Starting the Lanczos iteration ------------------------
        do k=1,mlanc
          write(6,*) " "
          write(6,*) "******** Lanczos iteration =",k,"***********"

!  Checking loss of orthogonality with the frobenius norm
          allocate(iden(k,k),Pj(m,k))
          iden = 0.0
          Pj = 0.0
          Pj = Planc(:,1:k)
          write(6,*) "shape Pj=",shape(Pj)
          write(6,*) " "  
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
          write(6,*) "frob=",fterm,"alt",norm2(real(iden))
          write(6,*) "^^^^^^^^^^^^^^^^^^"
          write(6,*) " "
          fnorm(k) = fterm 
          deallocate(iden,Pj)

          ! doing A*v_j
          t=matmul(A,q_j) 
          
          ! u_(j+1)
          u = t - beta(k)*q_jm1

          alpha(k)=dot_product(q_j,u)
          write(6,*) ""
          write(6,*) "alpha(k) =",alpha(k)
          write(6,*) ""

          ! Adding second term to u_(j+1)
          u = u - alpha(k)*q_j 

          ! New norm
          lnorm = dot_product(u,u)
          lnorm = sqrt(lnorm)
          write(6,*) "Norm after Lanczos recurrence =",lnorm
          write(6,*) " "

!!  Starting Reorthogonalization, full reorthog
! Traditional way
!          do j=1,k
!            q_k = Planc(:,j)
!            lterm = dot_product(u,q_k)
!            do i=1,m
!              u(i) = u(i) -lterm*q_k(i)
!            enddo 
!          enddo

! Alternative method, using I-Q_k-1*Q^T_k-1
          allocate(iden(m,m),Pj(m,k-1))
          iden = 0.0
          Pj = 0.0
          Pj = Planc(:,1:k-1)
          write(6,*) "shape Pj=",shape(Pj)
          write(6,*) " "  

!  Set the identity matrix
          do f = 1,m
            iden(f,f) = 1.0d0
          enddo    
          iden = iden-matmul(Pj,transpose(Pj))

! do the Gram-Schmidt orthog
          u = matmul(iden,u)
          deallocate(iden,Pj)

!  New norm after full reorthog
          lnorm = dot_product(u,u)
          lnorm = sqrt(lnorm)
          write(6,*) "Norm after Lanczos recurrence =",lnorm
          write(6,*) " "

!  Update beta
          if (k .ne. mlanc) then
          beta(k+1)=lnorm
!            write(6,*) ""
!            write(6,*) "for iteration ",k,"beta =",lnorm
!            write(6,*) ""
          endif

!  Updating v_(j-1) for the next run (j+1) 
          q_jm1 = q_j

!  Updating v_(j+1), to be used at start of (j+1) iteration, can use the tol here also
          if (k .ne. mlanc) then
            q_j = (1/lnorm)*u

          !  Add new vectors 
            Planc(:,k+1) = q_j
          else
             
          endif

        enddo ! END OF LANCZOS LOOP

        return
        end subroutine lanczos_herm_re_full

 
        subroutine lanc_herm_re_minpro(A,v0,m,mlanc,alpha,beta,fnorm,&
                                           estnorm,chknorm)
        implicit none
! Lanczos algorithm with real matrix, doing
!  minimal Partial reorthogonalization(minpro) following closely the paper
! "The Lanczos ALgorithm with Partial Reorthogonalization", H D Simon
        integer :: i,j,k,l,m,n,mlanc,f,chk
        real(kind=8) :: b0,temp1,tol,fnorm(mlanc),fterm,lterm, &
                        estnorm(mlanc),term,frob,eps1,reps
        real(kind=8) :: u0,lnorm,DLANGE
        real(kind=8) :: rand_num(m)
        real,parameter :: pi=3.141592654d0,eps=1.1E-16
        complex,parameter :: ii=complex(0.d0,1.d0)
        
        real(kind=8) :: alpha(mlanc),beta(mlanc),c(mlanc),work(mlanc)
        real(kind=8):: v0(m),vm1(m),v(m),s(m),t(m),u(m),q_j(m),q_k(m),&
                       q_jm1(m),q_jp1(m),w(m)
        
        real(kind=8) :: A(m,m)

!  Matrices for calculating reorthogonality
        real(kind=8) :: Planc(m,mlanc)
        real(kind=8), allocatable :: iden(:,:),Pj(:,:),Pj2(:,:),w_jk(:)
        integer,parameter :: sedsz=33
        integer :: seed(sedsz),chknorm(mlanc),oterm,iterm,orthtot 
        integer, allocatable :: orthindx(:),orthopt(:)
        
        
        reps = sqrt(eps)
        eps1 = reps !1000*eps
        write(6,*) "||A||_F =",norm2(A)
        write(6,*) "eps=",eps,"eps*N=",eps1
        write(6,*) " "

!   Initializing arrays to 0
        alpha=0.0
        beta=0.0
        c=0.0
        q_jm1=0.0
        q_j=0.0
        q_jp1=0.0
        t=0.0
        u=0.0
        v = 0.0
        w = 0.0
        Planc = 0.0

!  Normalizing the initial vector , v_0
        u0=dot_product(v0,v0)
        u0=sqrt(u0)
        write(6,*) "--- The Euclidean norm =",u0
        q_j=v0/u0
        beta(1)=0

!  Collecting q_1 and p_1 to add to P and Q
        Planc(:,1) = q_j


!----------  Starting the Lanczos iteration ------------------------

! total # of orthogonalizations
        orthtot = 0
        do k=1,mlanc
          write(6,*) " "
          write(6,*) "******** Lanczos iteration =",k,"***********"
          
          if (k .ne. 1)then
            q_jm1 = Planc(:,k-1)
          else
            q_jm1 = 0.0
          endif
          q_j = Planc(:,k)

!  Checking loss of orthogonality with the frobenius norm
          allocate(Pj(m,k))
!          iden = 0.0
          Pj = 0.0
          Pj = Planc(:,1:k)
          write(6,*) "shape Pj=",shape(Pj)
          write(6,*) " "

!  Calculate Frobenius norm  
          call check_orthog_frob(Pj,m,k,frob)
          fnorm(k) = frob
          deallocate(Pj)

          t=matmul(A,q_j) !doing A*q_j
          
          ! u_(j+1)
          u = t - beta(k)*q_jm1 
          alpha(k)=dot_product(q_j,u)
          write(6,*) ""
          write(6,*) "alpha(k) =",alpha(k)
          write(6,*) ""
          
          ! u_(j+1)
          u = u - alpha(k)*q_j 

          !  New norm
          lnorm = dot_product(u,u)
          lnorm = sqrt(lnorm)
          write(6,*) "Norm after Lanczos recurrence =",lnorm
          write(6,*) " "

!  Update beta
          if (k .ne. mlanc) then
          beta(k+1)=lnorm
!            write(6,*) ""
!            write(6,*) "for iteration ",k,"beta =",lnorm
!            write(6,*) ""
          endif

!!  Updating v_(j-1) for the next run (j+1) 
!          do i=1,m
!            q_jm1(i)=q_j(i)
!          enddo

!  Updating v_(j+1), to be used at start of (j+1) iteration, can use the tol here also
          if (k .ne. mlanc) then
            do i=1,m
              q_jp1(i)=(1/lnorm)*u(i)
!              q_j(i)=u(i)
            enddo

!  Maintain local orthog and add new  new vectors 
!            term = dot_product(q_jp1,q_j)
!            write(6,*) " "
!            write(6,*) " Local orthog. =",term
!            write(6,*) " "
!            q_jp1 = q_jp1 -term*q_j
            Planc(:,k+1) = q_jp1

!  Checking norm after adding a new Lanczos vector
            write(6,*) " ~~~~~ added new lanczos vec ~~~~~" 
            call check_orthog_frob(Planc(:,1:k+1),m,k+1,frob)
            write(6,*) " "
          else
             
          endif


! Alternative method, using I-Q_k-1*Q^T_k-1
          if (k .ne. mlanc) then
            allocate(iden(m,m),Pj(m,k),w_jk(k))
            iden = 0.0
            Pj = 0.0
            w_jk = 0.0
            Pj = Planc(:,1:k)
            write(6,*) "shape Pj=",shape(Pj)
            write(6,*) " "  

! Alternative check of orthogonality: <q_j+1,q_k> , k=1,...,j
            v =Planc(:,k+1)
            w_jk = matmul(transpose(Pj),v)
            estnorm(k) = maxval(abs(w_jk)) 
            write(6,*) " "
            write(6,*) "orthog, <q_j+1,q_k>, max",maxval(abs(w_jk))
            write(6,*) " "
            write(6,*) w_jk
            write(6,*) " "

! Check for partial orthog
            allocate(orthopt(k))
            orthopt=0
            do i=1,k
              term = abs(w_jk(i))
              if (term .ge. eps1) then
!                write(6,*) " "
!                write(6,*) " vector",i,"lost orth, w_jk=",term,eps1
!                write(6,*) " "
                orthopt(i) = 1
              endif
            enddo
            oterm = sum(orthopt)
!            write(6,*) "---# vecs that lost orthog",sum(orthopt),"---"
!            write(6,*) " "

!  Get the index of these vecs
            allocate(orthindx(oterm))
            orthindx = 0
            iterm = 1
            do i=1,k
              if (orthopt(i) .eq. 1) then
                orthindx(iterm) = i
                iterm = iterm + 1
              endif
            enddo
!            write(6,*) " The vecs to orthog against are:"
!            write(6,*) " "
!            write(6,*) orthindx
            deallocate(orthopt)

!!  Carry out full reorthog with those vecs that have lost semiorthog
            if (oterm .ge. 1) then
              write(6,*) "doing partial Gram-Schmidt as orthog failed"
              write(6,*) " "

! First do classic GS on the Pj,
              v = 0.0
              w = 0.0
              v = Planc(:,k+1)
              frob = 0.0
              call gram_schmidt_class_sing_pro(Planc(:,1:k),m,&
                                      k,oterm,v,w,orthindx,frob)
!              write(6,*) "After CGS, frobenius norm =",frob

!  Update Planc(:,k+1)
              Planc(:,k+1) = w
              call check_orthog_frob(Planc(:,1:k+1),m,k+1,frob)
              write(6,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
              write(6,*) "Lost partial orthogonality at lanc iter",k
              write(6,*) " "
              write(6,*) "Norm after GS = ",frob
              write(6,*) " "
              write(6,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
              write(6,*) " "
            endif
            deallocate(orthindx)
            chknorm(k) = oterm
            
! update total of orthogonalizations
            orthtot = orthtot + oterm

!! Check for semi-orthogonality, sqrt(eps)
!            chk = 0
!            do i=1,k
!              term = abs(w_jk(i))
!              if (term .ge. sqrt(eps)) then
!                write(6,*) " "
!                write(6,*) " vector",i,"lost orth, w_jk=",term,sqrt(eps)
!                write(6,*) " "
!                chk = chk+1
!              endif
!            enddo
!            write(6,*) "---# of vecs that lost semi orthog ",chk,"---"
!            write(6,*) " "
!!            chknorm(k) = chk
!!!  Carry out full reorthog if some vecs have lost semiorthog
!            if (chk .ge. 1) then
!              write(6,*) "doing full Gram-Schmidt as semi-orthog failed"
!              write(6,*) " "
!              v = 0.0
!              w = 0.0
!              v = Planc(:,k+1)
!              frob = 0.0
!              call gram_schmidt_class_sing(Planc(:,1:k),m,k,v,w,frob)
!              write(6,*) "After CGS, frobenius norm =",frob
!              frob = 0.0
!
!!  Update Planc(:,k+1)
!              Planc(:,k+1) = w
!              call check_orthog_frob(Planc(:,1:k+1),m,k+1,frob)
!              write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!              write(6,*) "Lost semi orthogonality at lanc iter",k
!              write(6,*) " "
!              write(6,*) "Norm after GS = ",frob
!              write(6,*) " "
!              write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!              write(6,*) " "
!
!            endif

            deallocate(iden,Pj,w_jk)
          endif  ! condition k .ne. mlanc
        enddo ! END OF LANCZOS LOOP

!  Print out final total orthogonalizations
        write(6,*) " "
        write(6,*) "****** TOTAL ORTHOGS =",orthtot,"******"
        write(6,*) " "
        write(6,*) " "
        return
        end subroutine lanc_herm_re_minpro



        subroutine lanczos_herm(A,v0,m,mlanc,alpha,beta,fnorm,Planc,full)
! Does a Lanczos for complex Hermitian matrix
! For a good working example: mi=6, mj=5, m=1000, mlanc=10,
! choose s=1.0+1.0*ii as starting lanczos vectors.
! The estimated spectrum will be close to the exact for the extreme eigen values
! Does full reorthog

        implicit none
        integer :: i,j,k,l,m,n,mlanc,f
        real(kind=8) :: b0,temp1,tol,fnorm(mlanc),fterm,lterm
        real(kind=8) :: u0,lnorm,DLANGE
        real,parameter :: pi=3.141592654d0
        complex,parameter :: ii=complex(0.d0,1.d0)
        
        real(kind=8) :: alpha(mlanc),beta(mlanc+1),c(mlanc+1),work(mlanc)
        complex(kind=8):: v0(m),vm1(m),v(m),s(m),t(m),u(m)
        
        complex(kind=8) :: A(m,m)

!  Matrices for calculating biorthogonality
        complex(kind=8) :: Planc(m,mlanc)
        complex(kind=8), allocatable :: iden(:,:),Pj(:,:),temp_mat(:,:)
        logical :: full

!   Allocating the arrays needed: naming convention follows SLEPc
!   Technical paper, "Lanczos methods in SLEPc", 
!   vm1 = v_(j-1),v=v_j, t=A*v_j,u=u_(j+1)=A*v_j - beta_j*v_(j-1)
 

!   Initializing arrays to 0
        Planc = 0.0
        alpha=0.0
        beta=0.0
        c=0.0
        vm1=0.0
        s=0.0
        t=0.0
        u=0.0

! Normalizing the initial vector , v_1
        u0=dot_product(v0,v0)
        u0=sqrt(u0)
!        write(6,*) "--- The Euclidean norm =",u0,u0**2
!        print *
!        print *, "******* Doing the Hermitian Lanczos algorithm *********"
!        print *
!        write(6,*) "--- The Euclidean norm =",u0
        v=v0/u0

        beta(1)=0

!  Collecting q_1 and p_1 to add to P and Q
        Planc(:,1) = v


!----------  Starting the Lanczos iteration ------------------------
        do k=1,mlanc
          write(6,*) " "
          write(6,*) "******** Lanczos iteration =",k,"***********"

!  Checking biorthogonality with the frobenius norm
          allocate(iden(k,k),Pj(m,k))
          Pj = Planc(:,1:k)
          iden = 0.0d0
!          write(6,*) "shape Pj=",shape(Pj)
!          write(6,*) " "  
!  Set the identity matrix
          do f = 1,k
            iden(f,f) = 1.0d0
          enddo    
          iden = iden-matmul(conjg(transpose(Pj)),Pj)

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
          write(6,*) "Loss of orthog. =",fterm
          write(6,*) "^^^^^^^^^^^^^^^^^^"
          write(6,*) " "
          fnorm(k) = fterm 
          deallocate(iden,Pj)

          t=matmul(A,v) !doing A*v_j
           
         ! u_(j+1)
          u = t - beta(k)*vm1
          alpha(k)=dot_product(v,u)
!!          write(6,*) ""
!          write(6,*) "alpha(k) =",alpha(k)
!          write(6,*) ""

          ! Adding second term to u_(j+1)
          u = u - alpha(k)*v 
!          print *
!          print *,"<u,u> = ",dot_product(u,u)
!          print *

          if (full) then
!            print *
!            print *,"----- Doing full reorthog -----"
!            print *      
            allocate(iden(m,m),Pj(m,k))
            iden = 0.0d0
            Pj = 0.0d0
            Pj = Planc(:,1:k)

            ! Set the identity matrix
            do i=1,m
              iden(i,i) = 1.0d0
            enddo
            iden = iden - matmul(Pj,conjg(transpose(Pj)))

            ! Reorthogonalization
            u = matmul(iden,u)
            deallocate(iden,Pj)
          end if

!  New norm
          lnorm = dot_product(u,u)
          lnorm = sqrt(lnorm)
!          if (k .ne. mlanc) then
          beta(k+1)=lnorm
!!            write(6,*) ""
!          write(6,*) "for iteration ",k,"beta =",lnorm
!!            write(6,*) ""
!!          endif

!  Updating v_(j-1) for the next run (j+1) 
          do i=1,m
            vm1(i)=v(i)
          enddo

!  Updating v_(j+1), to be used at start of (j+1) iteration, can use the tol here also
!          print *
!          print *,"<u,u> = ",dot_product(u,u)
!          print *
          if (k .ne. mlanc) then
            v = (1/lnorm)*u
!  Add new vectors 
            Planc(:,k+1) = v
          else
             
          endif
!          print *
!          print *,"<v_j+1,v_j+1> = ",dot_product(v,v)
!          print *
        enddo ! END OF LANCZOS LOOP

        return
        end subroutine lanczos_herm 

        
        subroutine lanczos_non_symm(A,v0,m,mlanc,alpha,beta,delta)
        implicit none
! Carries out Lanczos algorithm for a real nonsymmetric matrix.

        integer :: i,j,k,l,m,n,mlanc,f
        real(kind=8) :: b0,temp1,tol,fnorm(mlanc),fterm
        real(kind=8) :: u0,lnorm
        real,parameter :: pi=3.141592654d0
        complex,parameter :: ii=complex(0.d0,1.d0)
        
        real(kind=8) ::alpha(mlanc),beta(mlanc+1),delta(mlanc+1), &
                     &  lnorm1,lnorm2,nsgn,one,term,A(m,m)
        
        real(kind=8) :: v0(m),vm1(m),wm1(m),v(m),w(m),s(m),tv(m),tw(m), &
                                   vp(m),wp(m),term1

!   Allocating the arrays needed: naming convention follows Y Saad
!    paper, "The Lanczos Biorthogonalization Algorithm and other Oblique
!    Projection Methods For Solving Large Unsymmetric Systems. " 
!   vm1 = v_(j-1),v=v_j, tv=A*v_j, tw=A^H*w_j, vp=v^_(j+1)=A*v_j - beta_j*v_(j-1)
!   wp = w^_(j+1)=A^H*w_j-alpha_j*w_j-delta_j*w_(j-1)
 
       
!   Initializing arrays to 0
        alpha=0.0
        beta=0.0
        delta=0.0
        vm1=0.0
        wm1=0.0
        v=0.0
        w=0.0
        tv=0.0
        tw=0.0
        vp=0.0
        wp=0.0
        one=1.0

!  Normalizing the initial vector , v_0
        u0=dot_product(v0,v0)
        u0=sqrt(u0)
!        write(6,*) "--- The Euclidean norm =",u0
        v=v0/u0
        w=v


!----------  Starting the Lanczos iteration ------------------------
        do k=1,mlanc
          if (k==1) then
             beta(k)=0.0
             delta(k)=0.0
          end if
!          write(6,*) "******** Lanczos iteration =",k,"***********"
!          write(6,*)"alpha=",alpha(k),"beta =",beta(k),"delta=",delta(k)
!          write(6,*) ""
          tv=matmul(A,v) !doing A*v_j
          tw=matmul(transpose(A),w) ! doing A^T*w_j

          ! v^_(j+1), w^_(j+1)
          vp = tv - beta(k)*vm1 
          wp = tw - delta(k)*wm1 

          !          term=dot_product(w,tv) ! Following the Y. Saad paper
          term=dot_product(w,tv) ! Following the SLEPc paper for Hermitian case
!          term = sqrt(term1*conjg(term1))
          alpha(k) = term
          write(6,*) ""
!          write(6,*) "alpha(k) =",alpha(k)
!          write(6,*) "After calculating alpha, coeffs are "
          write(6,*)"alpha=",alpha(k),"beta =",beta(k),"delta=",delta(k)
          write(6,*) ""

          ! u_(j+1), w_(j+1)
          vp = vp - alpha(k)*v 
          wp = wp - alpha(k)*w

          !  New norm
          term1 = dot_product(vp,wp)
          lnorm1 = dot_product(vp,wp) ! following Y. Saad's paper
          nsgn = sign(one,lnorm1)
!          nsgn = 1.0 
!          write(6,*) ""
!          write(6,*) "(v^_(j+1),w^_(j+1))=",term1,"lnorm1=",lnorm1, &
!                   & "sign =",nsgn
!          write(6,*) ""
          lnorm1 = sqrt(abs(lnorm1))
          lnorm2 = nsgn*lnorm1  ! following Y. Saad's paper
!          lnorm2 = dot_product(wp,wp)  ! following yambo paper
!          lnorm2 = sqrt(abs(lnorm2))
          delta(k+1) = lnorm1 ! following Y. Saad's paper
!            delta(k+1)= lnorm1/lnorm2 ! following yambo paper
          beta(k+1)=lnorm2 ! following Y. Saad's paper
!            write(6,*) ""
!            write(6,*) "for iteration ",k,"beta =",lnorm
!            write(6,*) ""

!  Updating v_(j-1) for the next run (j+1) 
          vm1 = v
          wm1 = w

!  Updating v_(j+1), to be used at start of (j+1) iteration, can use the tol here also
          if (k .ne. mlanc) then
            v = (1.0/lnorm1)*vp ! v_(j+1)=v^_(j+1)/delta_(j+1)
            w = (1.0/lnorm2)*wp ! w_(j+1)=w^_(j+1)/beta_(j+1)
          else
             
          endif

        enddo ! END OF LANCZOS LOOP
        return
        end subroutine lanczos_non_symm


        subroutine lanczos_non_herm(A,v0,m,mlanc,zalpha,zbeta,zgamma,& 
                                   zborth,fnorm,paige,reorth,mult_type)
        implicit none
  
! Carries out Lanczos algorithm on a general non Hermitian complex matrix. 
! Uses complex coefficients.
!
!--------------------------------------------------------------------------------
! Follows the non-hermitian eigen vaalue problem algorithm from 
! https://netlib.org/utk/people/JackDongarra/etemplates/node245.html#bd:lanalg
! Actual algorithm only uses two recurrence terms
!-------------------------------------------------------------------------------
        integer :: i,j,k,l,m,n,f,mlanc,mhalf,orthfreq,orthtot
        real(kind=8) :: b0,temp1,tol,eps1
        real(kind=8) :: u0,csgn,rand_num(m),fnorm(mlanc),fterm,frob
        real,parameter :: pi=3.141592654d0
        complex,parameter :: ii=complex(0.d0,1.d0)
        integer,parameter :: sedsz=33
        integer :: seed(sedsz) 
        complex(kind=8) ::zalpha(mlanc),zbeta(mlanc+1),zgamma(mlanc+1), &
                     &  lnorm1,lnorm2,term,term1,zborth(mlanc),tleft,&
                     &   tright,pterm,qterm
        complex(kind=8):: v0(m),pm1(m),qm1(m),p(m),q(m),s(m),tq(m),tp(m), &
                      & pp(m),qp(m),pp_j(m),qp_j(m),p_k(m),q_k(m) ! u is now vp/ v^_(j+1) ,wp is
                                    !  w^_(j+1)
        
        complex(kind=8) :: A(m,m)
        logical :: paige, biortho,reorth  ! to decide whether or not to do a Paige like step 
        character (10) :: mult_type ! methods of doing the matrix-matrix-vector operation

!  Matrices for calculating biorthogonality
        complex(kind=8) :: Planc(m,mlanc), Qlanc(m,mlanc)
        complex(kind=8), allocatable :: iden(:,:),Pj(:,:),Qj(:,:), &
                                        iden1(:,:),iden2(:,:)
  ! Lapack routines
        complex( kind = kind( 1.0d0 ) ):: ZDOTC

!   Allocating the arrays needed: naming convention follows Y Saad
!    paper, "The Lanczos Biorthogonalization Algorithm and other Oblique
!    Projection Methods For Solving Large Unsymmetric Systems. " 
!   vqm1 = q_(j-1),v=v_j, tq=A*q_j, tp=A^H*p_j, qp=p^_(j+1)=A*q_j - gamma_j*q_(j-1)
!   pp = p^_(j+1)=A^H*p_j-conjg(alpha_j)*p_j-conjg(beta_j)*p_(j-1).
!   Adapting Saad's method to the complex case and choosing beta*delta
!   accordinf to the yambo paper "Implementation and testing of
!   Lanczos-based algorithms for Random-Phase Approximation
!   eigenproblems."
 
!   Initializing arrays to 0
        zalpha=0.0d0
        zbeta=0.0d0
        zgamma=0.0d0
        zborth = 0.0d0
        fnorm = 0.0d0
        pm1=0.0d0
        qm1=0.0d0
        p=0.0d0
        q=0.0d0
        tp=0.0d0
        tq=0.0d0
        pp=0.0d0
        qp=0.0d0



!  Normalizing the initial vector , v_0
        u0=dot_product(v0,v0)
        u0=sqrt(u0)
!        write(6,*) "--- The Euclidean norm =",u0,u0**2
        p=v0/u0
        q=p

!  Collecting q_1 and p_1 to add to P and Q
        Planc = 0.0d0
        Qlanc = 0.0d0
        Planc(:,1) = p
        Qlanc(:,1) = q

        ! orthogs info
        orthtot = 0
        orthfreq = 0
!----------  Starting the Lanczos iteration ------------------------
        do k=1,mlanc
          if (k==1) then
             zbeta(k)=0.0d0
             zgamma(k)=0.0d0
          end if
          write(6,*) "******** Lanczos iteration =",k,"***********"

!  Checking biorthogonality with the frobenius norm
          allocate(Pj(m,k),Qj(m,k))
!          iden = 0.0
          Pj = 0.0d0
          Qj = 0.0d0
          Pj = Planc(:,1:k)
          Qj = Qlanc(:,1:k)   
!          write(6,*) "Pj=",Pj
!          write(6,*) " "
!          write(6,*) "Qj=",Qj
!          write(6,*) " "  
          call zcheck_orthog_frob(Pj,Qj,m,k,frob)
          write(6,*) " "
          write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^" 
          write(6,*) "Loss of biorthog. =",frob 
          write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^"
          write(6,*) " "
          fnorm(k) = frob
          deallocate(Pj,Qj)

          tq=matmul(A,q) !doing A*v_j
          tp=matmul(transpose(conjg(A)),p) ! doing A^H*w_j
          term=dot_product(p,tq) ! Following the SLEPc paper for Hermitian case
!          write(6,*) ""
!          write(6,*) "====== Before Lanczos, overlap ========"
!          write(6,*) " <p|A|q> =", term
!          write(6,*) ""

          ! Calculating alpha
          qp = tq - zgamma(k)*qm1 ! v^_(j+1)
          pp = tp - conjg(zbeta(k))*pm1 ! v^_(j+1)

          term=dot_product(p,qp) ! Following the SLEPc paper for Hermitian case
          zalpha(k) = term
          write(6,*) ""
          write(6,*) "alpha(k) =",zalpha(k)
!          write(6,*) "After calculating alpha, coeffs are "
!             write(6,*) ""
!             write(6,*)'alpha=',zalpha(i)
          write(6,*)""
          write(6,*)'beta(k) =',zbeta(k)
          write(6,*)""
          write(6,*)'delta(k) =',zgamma(k)
          write(6,*) ""

!   Including this Paige-like regularization, from hermitian lanczos, for better convergence
          if (paige) then
              qp = qp - zalpha(k)*q ! u_(j+1)  Paige-like regularization
              pp = pp - conjg(zalpha(k))*p ! u_(j+1), Paige-like reg
          else
              qp = tq - zalpha(k)*q ! u_(j+1) , netlib Alg 7.42
              pp = tp - conjg(zalpha(k))*p ! u_(j+1), netlib Alg 7.42
          endif

!!  Re-biorthogonalizing the Lanczos vectors, full reorthog
          if (reorth) then         
           
          select case (mult_type) 
            case ('pref')
             print *
             print *," using matrix-matrix-vector case ",mult_type        

! Method 1 : Preferred/most robust way to do y to do Gram-Schmidt with (I-Q_k*(P_k)^H)
! Uses matrices and does the matrix-matrix-vector product
                    allocate(iden1(m,m),iden2(m,m),Pj(m,k),Qj(m,k))
                    iden1 = 0.0d0
                    iden2 = 0.0d0
                    Pj = 0.0d0
                    Qj = 0.0d0
                    Pj = Planc(:,1:k)
                    Qj = Qlanc(:,1:k)   
        !          write(6,*) "Pj=",Pj
        !          write(6,*) " "
        !          write(6,*) "Qj=",Qj
        !          write(6,*) " "  
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

            case ('oneshot')        
             print *
             print *," using matrix-matrix-vector case ",mult_type   

!! Method 2 : oneshot
!! Uses implicit multiplication of matrices and does not form intermediate products
!! Less costly, but for some more challenging non Hermitian matrices  will not work so well
                    allocate(Pj(m,k),Qj(m,k))
                    Pj = 0.0d0
                    Qj = 0.0d0
                    Pj = Planc(:,1:k)
                    Qj = Qlanc(:,1:k)   

                    pp_j=0.0d0
                    qp_j=0.0d0
                    if (k .ne. 1) then
                       pp_j = matmul(Pj,matmul(conjg(transpose(Qj)),pp) )
                       qp_j = matmul(Qj,matmul(conjg(transpose(Pj)),qp) )
!                       print *
!                       print *," <pp_j,qp_j> = ",dot_product(pp_j,qp_j)
                    ! Updating pp and qq
                       pp = pp - pp_j
                       qp = qp - qp_j
                    else
                       pp = pp
                       qp = qp
                    endif        
                    deallocate(Pj,Qj)
             case ('naive')      
             print *
             print *," using matrix-matrix-vector case ",mult_type    

!! Method 3 : the naive way, using just dot products
!! Uses implicit multiplication of matrices and does not form intermediate products
!! Less costly, but for some more challenging non Hermitian matrices  will not work so well
                    pp_j = pp
                    qp_j = qp

                    do j=1,k
                       p_k = Planc(:,j) ! initially A, B
                       q_k = Qlanc(:,j)
                       pterm = dot_product(pp_j,q_k)
                       qterm = dot_product(qp_j,p_k)
        !               pterm = ZDOTC(m,pp_j,1,q_k,1)
        !               qterm = ZDOTC(m,qp_j,1,p_k,1)
                       pp_j = pp_j - pterm*p_k
                       qp_j = qp_j - conjg(qterm)*q_k
                    enddo  
                    pp = pp_j
                    qp = qp_j
             case default
             print *
             print *," using matrix-matrix-vector case default = 'pref'",mult_type      

! Method 1 : Preferred/most robust way to do y to do Gram-Schmidt with (I-Q_k*(P_k)^H)
! Uses matrices and does the matrix-matrix-vector product
                    allocate(iden1(m,m),iden2(m,m),Pj(m,k),Qj(m,k))
                    iden1 = 0.0d0
                    iden2 = 0.0d0
                    Pj = 0.0d0
                    Qj = 0.0d0
                    Pj = Planc(:,1:k)
                    Qj = Qlanc(:,1:k)   
        !          write(6,*) "Pj=",Pj
        !          write(6,*) " "
        !          write(6,*) "Qj=",Qj
        !          write(6,*) " "  
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
            end select ! end case selection for matrix-matrix-vector product

            ! update orthogs info
            orthfreq = orthfreq + 1
            orthtot = orthtot +2*k
          endif ! rebiorthog condition

!!  New norm
          term1 = dot_product(qp,pp) ! following netlib Alg 7.13

!  A measure of the biorthogonality
          zborth(k) = term1

          term = dot_product(qp,qp) ! following netlib Alg 7.13
!          write(6,*) ""
!          write(6,*) "(r^_(j+1),r^_(j+1))=",term
!          write(6,*) "norm=",sqrt(term) 
!          write(6,*) ""
          term = dot_product(pp,pp) ! following netlib Alg 7.13
!          write(6,*) ""
!          write(6,*) "(s^_(j+1),s^_(j+1))=",term
!          write(6,*) "norm=",sqrt(term) 
!          write(6,*) ""
!          write(6,*) ""
!          write(6,*) "(r^_(j+1),s^_(j+1))=",term1 
!          write(6,*) ""
          lnorm1 = sqrt(abs(term1))  ! following yambo paper
          lnorm2 = conjg(term1)/lnorm1
!          write(6,*) "|(r^_(j+1),s^_(j+1))|=",lnorm1 !,"sign =",csgn
!          write(6,*) "lnorm2 =",lnorm2
!          write(6,*) ""
          zbeta(k+1)= lnorm1 ! following yambo paper
          zgamma(k+1)=lnorm2 ! following Y. Saad's paper
!            write(6,*) ""
!            write(6,*) "for iteration ",k,"beta =",lnorm
!            write(6,*) ""

!  Updating v_(j-1) for the next run (j+1) 
          qm1 = q
          pm1 = p

!  Updating v_(j+1), to be used at start of (j+1) iteration, can use the tol here also
          if (k .ne. mlanc) then
            q = (1.0d0/lnorm1)*qp! v_(j+1)=v^_(j+1)/delta_(j+1)
            p = (1.0d0/conjg(lnorm2))*pp ! w_(j+1)=w^_(j+1)/beta_(j+1)

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

!  Checking the local biorthogonality
!          term = dot_product(q,p)
!          write(6,*) "========== biorthogonality ========="
!          write(6,*) " < p_i, q_i > =",term
!          write(6,*) " "
          if (reorth) then
             write(222,*)k,k,2*k,orthfreq,orthtot
          else
             write(222,*)k,0,2*0,orthfreq,orthtot
          endif
                  
        enddo ! END OF LANCZOS LOOP
        write(6,*) " "
        write(6,*) "****** FULL RE ORTHOGS = ",(mlanc+1)*mlanc,"******"
        write(6,*) " "
        write(6,*) " "
        return
        end subroutine lanczos_non_herm



        subroutine lanczos_non_herm_minprbo(A,v0,m,mlanc,zalpha,zbeta,zgamma,& 
                                   zborth,fnorm,paige,reorth,eps1,mult_type)
        implicit none
!--------------------------------------------------------------------------------
! Carries out the non-hermitian Lanczos algorithm on a complex matrix. Does a minimal
! partial re biorthogonalization strategy (minPRBO).
! Follows the non-hermitian eigen value problem algorithm from 
! https://netlib.org/utk/people/JackDongarra/etemplates/node245.html#bd:lanalg
! Actual algorithm only uses two recurrence terms.
!
! Does the minimal partial reorthogonalization on the unnormalized pp and qp
!-------------------------------------------------------------------------------
        integer :: i,j,k,l,m,n,f,mlanc,mhalf,iterm,jterm
        real(kind=8) :: b0,temp1,tol,eps1,term_re,term_im
        real(kind=8) :: u0,csgn,rand_num(m),fnorm(mlanc),fterm,frob
        real,parameter :: pi=3.141592654d0,eps=1.11E-16
        complex,parameter :: ii=complex(0.d0,1.d0)
        integer,parameter :: sedsz=33
        integer :: seed(sedsz),orthtot,orthfreq,poterm,qoterm,chk,oterm 
        complex(kind=8) ::zalpha(mlanc),zbeta(mlanc+1),zgamma(mlanc+1), &
                     &  lnorm1,lnorm2,term,term1,zborth(mlanc),pterm,qterm
        complex(kind=8):: v0(m),pm1(m),qm1(m),p(m),q(m),s(m),tq(m),tp(m), &
                      & pp(m),qp(m),qo(m),po(m),vp(m),vq(m),wp(m),wq(m),pp_j(m),qp_j(m),p_k(m),q_k(m) ! u is now vp/ v^_(j+1) ,wp is
                                    !  w^_(j+1)
        
        complex(kind=8) :: A(m,m)
        logical :: paige, reorth  ! to decide whether or not to do a Paige like step
        character (10) :: mult_type ! methods of doing the matrix-matrix-vector operation

!  Matrices for calculating biorthogonality
        integer, allocatable :: porthopt(:),qorthopt(:),porthindx(:),&
                                qorthindx(:)
        complex(kind=8) :: Planc(m,mlanc), Qlanc(m,mlanc)
        complex(kind=8), allocatable :: iden(:,:),Pj(:,:),Qj(:,:), &
                        iden1(:,:),iden2(:,:),&
                        pw_jk(:),qw_jk(:)

  ! Lapack routines
        complex( kind = kind( 1.0d0 ) ):: ZDOTC
                                      

!   Allocating the arrays needed: naming convention follows Y Saad
!    paper, "The Lanczos Biorthogonalization Algorithm and other Oblique
!    Projection Methods For Solving Large Unsymmetric Systems. " 
!   vqm1 = q_(j-1),v=v_j, tq=A*q_j, tp=A^H*p_j, qp=p^_(j+1)=A*q_j - gamma_j*q_(j-1)
!   pp = p^_(j+1)=A^H*p_j-conjg(alpha_j)*p_j-conjg(beta_j)*p_(j-1).
 
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


!  Normalizing initial vector , v_0
        u0=dot_product(v0,v0)
        u0=sqrt(u0)
!        write(6,*) "--- The Euclidean norm =",u0
        p=v0/u0
        q=p

!  Collecting q_1 and p_1 to add to P and Q
        Planc = 0.0
        Qlanc = 0.0
        Planc(:,1) = p
        Qlanc(:,1) = q
 
        write(6,*) " "
        write(6,*) " ****** eps1 =",eps1,"*******"
        write(6,*) " "

! total # of orthogonalizations and freq
        orthtot = 0
        orthfreq = 0

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
          call zcheck_orthog_frob(Pj,Qj,m,k,frob)
          write(6,*) " "
          write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^" 
          write(6,*) "Loss of biorthog. =",frob !,"alt="
          write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^"
          write(6,*) " "
          fnorm(k) = frob
          deallocate(Pj,Qj)

          tq=matmul(A,q) !doing A*v_j
          tp=matmul(transpose(conjg(A)),p) ! doing A^H*w_j
          term=dot_product(p,tq) ! Following the SLEPc paper for Hermitian case
!          write(6,*) ""
!          write(6,*) "====== Before Lanczos, overlap ========"
!          write(6,*) " <p|A|q> =", term
!          write(6,*) ""

!   Including this Paige-like regularization, from hermitian lanczos, for better convergence
          qp = tq - zgamma(k)*qm1 ! v^_(j+1)
          pp = tp - conjg(zbeta(k))*pm1 ! v^_(j+1)

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
          if (paige) then
            qp = qp - zalpha(k)*q ! u_(j+1)  Paige-like regularization
            pp = pp - conjg(zalpha(k))*p ! u_(j+1), Paige-like reg
          else
            qp = tq - zalpha(k)*q ! u_(j+1) , netlib Alg 7.42
            pp = tp - conjg(zalpha(k))*p ! u_(j+1), netlib Alg 7.42
          endif


! Checking for and doing partial re biorthogonality
          if (k .ne. mlanc) then
            allocate(Pj(m,k),Qj(m,k),pw_jk(k),qw_jk(k))
            Pj = 0.0
            Qj = 0.0
            pw_jk = 0.0
            qw_jk = 0.0
            Pj = Planc(:,1:k)
            Qj = Qlanc(:,1:k)
!            write(6,*) "shape Pj/Qj =",shape(Pj)
!            write(6,*) " "

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
!            write(6,*) " "
!            write(6,*) "orthog, <p_j+1,q_k>,max",maxval(abs(pw_jk))
!            write(6,*) " "
!            write(6,*) pw_jk
!            write(6,*) " "
!            write(6,*) "orthog, <q_j+1,p_k>,max",maxval(abs(qw_jk))
!            write(6,*) " "
!            write(6,*) qw_jk
!            write(6,*) " "

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
!                write(6,*) " "
!                write(6,*) " vector",i,"lost orth, pw_jk=",term,eps1
!                write(6,*) " "
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
!                write(6,*) " "
!                write(6,*) " vector",i,"lost orth, qw_jk=",term,eps1
!                write(6,*) " "
                qorthopt(i) = 1
              endif
            enddo
            qoterm = sum(qorthopt)
            write(6,*) "---# vecs that lost biortho",poterm,qoterm,"---"
            write(6,*) " "
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

               ! update orthfreq
             if (oterm .gt. 0) then
                  orthfreq = orthfreq + 1
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


! Construct the Pj and Qj
            if ((poterm .ge. 1) .or. (qoterm .ge. 1)) then
               write(6,*) " doing partial 2-sided GS as biorthog failed"
               write(6,*) " "


               allocate(Pj(m,poterm),Qj(m,qoterm))
               Pj = 0.0
               Qj = 0.0
               do i=1,poterm
                 iterm = porthindx(i)
                 jterm = qorthindx(i)
                 Pj(:,i) = Planc(:,iterm)
                 Qj(:,i) = Qlanc(:,jterm)
               enddo
 
! Carry out rebiorthog wrt those vectors that lost it

                  select case (mult_type) 
                    case ('pref')
                     print *
                     print *," using matrix-matrix-vector case ",mult_type        

        ! Method 1 : Preferred/most robust way to do y to do Gram-Schmidt with (I-Q_k*(P_k)^H)
        ! Uses matrices and does the matrix-matrix-vector product
                            allocate(iden1(m,m),iden2(m,m))
                            iden1 = 0.0
                            iden2 = 0.0
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
                            deallocate(iden1,iden2)

                    case ('oneshot')        
                     print *
                     print *," using matrix-matrix-vector case ",mult_type   

        !! Method 2 : oneshot
        !! Uses implicit multiplication of matrices and does not form intermediate products
        !! Less costly, but for some more challenging non Hermitian matrices  will not work so well

                            pp_j=0.0
                            qp_j=0.0
                            if (k .ne. 1) then
                               pp_j = matmul(Pj,matmul(conjg(transpose(Qj)),pp) )
                               qp_j = matmul(Qj,matmul(conjg(transpose(Pj)),qp) )
        !                       print *
        !                       print *," <pp_j,qp_j> = ",dot_product(pp_j,qp_j)
                            ! Updating pp and qq
                               pp = pp - pp_j
                               qp = qp - qp_j
                            else
                               pp = pp
                               qp = qp
                            endif    

                     case ('naive')      
                     print *
                     print *," using matrix-matrix-vector case ",mult_type    

        !! Method 3 : the naive way, using just dot products
        !! Uses implicit multiplication of matrices and does not form intermediate products
        !! Less costly, but for some more challenging non Hermitian matrices  will not work so well
                            pp_j = pp
                            qp_j = qp

                            do j=1,poterm
                               p_k = Pj(:,j) ! initially A, B
                               q_k = Qj(:,j)
                               pterm = dot_product(pp_j,q_k)
                               qterm = dot_product(qp_j,p_k)
                !               pterm = ZDOTC(m,pp_j,1,q_k,1)
                !               qterm = ZDOTC(m,qp_j,1,p_k,1)
                               pp_j = pp_j - pterm*p_k
                               qp_j = qp_j - conjg(qterm)*q_k
                            enddo  
                            pp = pp_j
                            qp = qp_j
                     case default
                     print *
                     print *," using matrix-matrix-vector case default = 'pref'",mult_type      

        ! Method 1 : Preferred/most robust way to do y to do Gram-Schmidt with (I-Q_k*(P_k)^H)
        ! Uses matrices and does the matrix-matrix-vector product
                            allocate(iden1(m,m),iden2(m,m))
                            iden1 = 0.0
                            iden2 = 0.0
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
                            deallocate(iden1,iden2)
                    end select ! end case selection for matrix-matrix-vector product
                    deallocate(Pj,Qj)
            endif   !actual rebiorthog

            deallocate(porthindx,qorthindx,pw_jk,qw_jk)
          endif  ! re biorthog condition

! 
          term1 = dot_product(qp,pp) ! following netlib Alg 7.13

!  A measure of the biorthogonality
          zborth(k) = term1

          term = dot_product(qp,qp) ! following netlib Alg 7.13
!          write(6,*) ""
!          write(6,*) "(r^_(j+1),r^_(j+1))=",term
!          write(6,*) "norm=",sqrt(term) 
!          write(6,*) ""
          term = dot_product(pp,pp) ! following netlib Alg 7.13
!          write(6,*) ""
!          write(6,*) "(s^_(j+1),s^_(j+1))=",term
!          write(6,*) "norm=",sqrt(term) 
!          write(6,*) ""
!          write(6,*) ""
!          write(6,*) "(r^_(j+1),s^_(j+1))=",term1 
!          write(6,*) ""
          lnorm1 = sqrt(abs(term1))  ! following yambo paper
          lnorm2 = conjg(term1)/lnorm1
!          write(6,*) "|(r^_(j+1),s^_(j+1))|=",lnorm1 
!          write(6,*) "(v^_(j+1),w^_(j+1))=",lnorm1
!          write(6,*) "lnorm2 =",lnorm2
!          write(6,*) ""
          zbeta(k+1)= lnorm1 ! following yambo paper
          zgamma(k+1)=lnorm2 ! following Y. Saad's paper
!            write(6,*) ""
!            write(6,*) "for iteration ",k,"beta =",lnorm
!            write(6,*) ""


!  Updating v_(j+1), to be used at start of (j+1) iteration, can use the tol here also
          if (k .ne. mlanc) then
            q=(1.0/lnorm1)*qp! v_(j+1)=v^_(j+1)/delta_(j+1)
            p=(1.0/conjg(lnorm2))*pp ! w_(j+1)=w^_(j+1)/beta_(j+1)

!  Add new vectors 
            Planc(:,k+1) = p
            Qlanc(:,k+1) = q

! Check for orthog after partial rebiorthog
            call zcheck_orthog_frob(Planc(:,1:k+1),Qlanc(:,1:k+1),& 
                                       m,k+1,frob)
            write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(6,*) " "
            write(6,*) "Norm after adding new Lanc vecs = ",frob
            write(6,*) " "
            write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(6,*) " "
          endif

! update total orthogs and save to file 
          orthtot = orthtot + 2*oterm
          write(222,*)k,oterm,2*oterm,orthfreq,orthtot
        enddo ! END OF LANCZOS LOOP
!  Print out final total orthogonalizations
        write(6,*) " "
        write(6,*) "****** TOTAL PARTIAL ORTHOGS =",orthtot," VS FULL RE ORTHOGS = ",(mlanc+1)*mlanc,"******"
        write(6,*) " "
        write(6,*) " "
        return
        end subroutine lanczos_non_herm_minprbo


        subroutine lanczos_non_herm_prbo(A,v0,m,mlanc,zalpha,zbeta,zgamma,& 
                                   zborth,fnorm,paige,reorth,eps1,prbo)
        implicit none
!--------------------------------------------------------------------------------
! Carries out the non-hermitian Lanczos algorithm on a complex matrix. Does a
! partial re biorthogonalization strategy (PRBO) as originally formulated,
! where once its deemed necessary to rebiorthogonalize, a full rebiorthogonalization of the 
! candidate vectors is carried out, followed by a second full rebiorthogonalization step in the 
! subsequent iteration, regardless of whether semi orthogonality is lost or not..
! Follows the non-hermitian eigen value problem algorithm from 
! https://netlib.org/utk/people/JackDongarra/etemplates/node245.html#bd:lanalg
! Actual algorithm only uses two recurrence terms.
!
! This is what's implemented in "ABLE: An adaptive block-Lanczos method for non-Hermitian
! eigenvalue problems."
! Does  partial reorthogonalization on the unnormalized pp and qp
!-------------------------------------------------------------------------------
        integer :: i,j,k,l,m,n,f,mlanc,mhalf,iterm,jterm,F_j,prbo,orthtot,orthfreq
        real(kind=8) :: b0,temp1,tol,eps1,term_re,term_im
        real(kind=8) :: u0,csgn,rand_num(m),fnorm(mlanc),fterm,frob
        real,parameter :: pi=3.141592654d0,eps=1.11E-16
        complex,parameter :: ii=complex(0.d0,1.d0)
        integer,parameter :: sedsz=33
        integer :: seed(sedsz),poterm,qoterm,chk,oterm 
        complex(kind=8) ::zalpha(mlanc),zbeta(mlanc+1),zgamma(mlanc+1), &
                     &  lnorm1,lnorm2,term,term1,zborth(mlanc),pterm,qterm
        complex(kind=8):: v0(m),pm1(m),qm1(m),p(m),q(m),s(m),tq(m),tp(m), &
                      & pp(m),qp(m),qo(m),po(m),vp(m),vq(m),wp(m),wq(m),pp_j(m),qp_j(m),p_k(m),q_k(m) ! u is now vp/ v^_(j+1) ,wp is
                                    !  w^_(j+1)
        
        complex(kind=8) :: A(m,m)
        logical :: paige, reorth  ! to decide whether or not to do a Paige like step
        character (10) :: mult_type ! methods of doing the matrix-matrix-vector operation

!  Matrices for calculating biorthogonality
        integer, allocatable :: porthopt(:),qorthopt(:),porthindx(:),&
                                qorthindx(:)
        complex(kind=8) :: Planc(m,mlanc), Qlanc(m,mlanc)
        complex(kind=8), allocatable :: iden(:,:),Pj(:,:),Qj(:,:), &
                        iden1(:,:),iden2(:,:),&
                        pw_jk(:),qw_jk(:)

  ! Lapack routines
        complex( kind = kind( 1.0d0 ) ):: ZDOTC
                                      

!   Allocating the arrays needed: naming convention follows Y Saad
!    paper, "The Lanczos Biorthogonalization Algorithm and other Oblique
!    Projection Methods For Solving Large Unsymmetric Systems. " 
!   vqm1 = q_(j-1),v=v_j, tq=A*q_j, tp=A^H*p_j, qp=p^_(j+1)=A*q_j - gamma_j*q_(j-1)
!   pp = p^_(j+1)=A^H*p_j-conjg(alpha_j)*p_j-conjg(beta_j)*p_(j-1).
 
!   Initializing arrays to 0
        zalpha=0.0d0
        zbeta=0.0d0
        zgamma=0.0d0
        zborth = 0.0d0
        fnorm = 0.0d0
        pm1=0.0d0
        qm1=0.0d0
        p=0.0d0
        q=0.0d0
        s=0.0d0
        tp=0.0d0
        tq=0.0d0
        pp=0.0d0
        qp=0.0d0
        po = 0.0d0
        qo = 0.0d0


!  Normalizing initial vector , v_0
        u0=dot_product(v0,v0)
        u0=sqrt(u0)
!        write(6,*) "--- The Euclidean norm =",u0
        p=v0/u0
        q=p

!  Collecting q_1 and p_1 to add to P and Q
        Planc = 0.0d0
        Qlanc = 0.0d0
        Planc(:,1) = p
        Qlanc(:,1) = q
 
        write(6,*) " "
        write(6,*) " ****** eps1 =",eps1,"*******"
        write(6,*) " "

! total # of orthogonalizations
        orthtot = 0
        orthfreq = 0
        F_j = 0 ! integer flag for second full rebiorthogonalization

!----------  Starting the Lanczos iteration ------------------------
        do k=1,mlanc
          if (k==1) then
             zbeta(k)=0.0d0
             zgamma(k)=0.0d0
          end if

          if (k .ne. 1)then
            pm1 = Planc(:,k-1)
            qm1 = Qlanc(:,k-1)
          else
            pm1 = 0.0d0
            qm1 = 0.0d0
          endif

          write(6,*) "******** Lanczos iteration =",k,"***********"

!  Checking biorthogonality with the frobenius norm
!          allocate(iden(k,k),Pj(m,k),Qj(m,k))
          allocate(Pj(m,k),Qj(m,k))
!          iden = 0.0
          Pj = 0.0d0
          Qj = 0.0d0
          Pj = Planc(:,1:k)
          Qj = Qlanc(:,1:k)   
          call zcheck_orthog_frob(Pj,Qj,m,k,frob)
          write(6,*) " "
          write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^" 
          write(6,*) "Loss of biorthog. =",frob !,"alt="
          write(6,*) "^^^^^^^^^^^^^^^^^^^^^^^^^"
          write(6,*) " "
          fnorm(k) = frob
          deallocate(Pj,Qj)

          tq=matmul(A,q) !doing A*v_j
          tp=matmul(transpose(conjg(A)),p) ! doing A^H*w_j
          term=dot_product(p,tq) ! Following the SLEPc paper for Hermitian case
!          write(6,*) ""
!          write(6,*) "====== Before Lanczos, overlap ========"
!          write(6,*) " <p|A|q> =", term
!          write(6,*) ""

!   Including this Paige-like regularization, from hermitian lanczos, for better convergence
          qp = tq - zgamma(k)*qm1 ! v^_(j+1)
          pp = tp - conjg(zbeta(k))*pm1 ! v^_(j+1)

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
          if (paige) then
            qp = qp - zalpha(k)*q ! u_(j+1)  Paige-like regularization
            pp = pp - conjg(zalpha(k))*p ! u_(j+1), Paige-like reg
          else
            qp = tq - zalpha(k)*q ! u_(j+1) , netlib Alg 7.42
            pp = tp - conjg(zalpha(k))*p ! u_(j+1), netlib Alg 7.42
          endif

          ! Doing the second rebiorthogonalization if one is needed in the previous iter
          ! This is PRBO-2
          if (prbo .eq. 2) then
             if (k .eq. F_j) then
                print *
                print *, " Doing a second full rebiorthogonalization "
                print *
                print *," iter = ",k," F_j =",F_j

! Method 1 : Preferred/most robust way to do y to do Gram-Schmidt with (I-Q_k*(P_k)^H)
! Uses matrices and does the matrix-matrix-vector product
                    allocate(iden1(m,m),iden2(m,m),Pj(m,k),Qj(m,k))
                    iden1 = 0.0d0
                    iden2 = 0.0d0
                    Pj = 0.0d0
                    Qj = 0.0d0
                    Pj = Planc(:,1:k)
                    Qj = Qlanc(:,1:k)   
        !          write(6,*) "Pj=",Pj
        !          write(6,*) " "
        !          write(6,*) "Qj=",Qj
        !          write(6,*) " "  
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
             endif        
          endif

! Checking for and doing partial re biorthogonality
          if (k .ne. mlanc) then
            allocate(Pj(m,k),Qj(m,k),pw_jk(k),qw_jk(k))
            Pj = 0.0d0
            Qj = 0.0d0
            pw_jk = 0.0d0
            qw_jk = 0.0d0
            Pj = Planc(:,1:k)
            Qj = Qlanc(:,1:k)
!            write(6,*) "shape Pj/Qj =",shape(Pj)
!            write(6,*) " "

!!  Calculate local orthog
            vp = 0.0d0
            vq = 0.0d0
            term1 = dot_product(qp,pp) ! following netlib Alg 7.13
            lnorm1 = sqrt(abs(term1))  ! following yambo paper
            lnorm2 = conjg(term1)/lnorm1
            vp=(1.0d0/conjg(lnorm2))*pp ! w_(j+1)=w^_(j+1)/beta_(j+1)
            vq=(1.0d0/lnorm1)*qp! v_(j+1)=v^_(j+1)/delta_(j+1)
            pw_jk = matmul(conjg(transpose(Qj)),vp)
            qw_jk = matmul(conjg(transpose(Pj)),vq)
!            write(6,*) " "
!            write(6,*) "orthog, <p_j+1,q_k>,max",maxval(abs(pw_jk))
!            write(6,*) " "
!            write(6,*) pw_jk
!            write(6,*) " "
!            write(6,*) "orthog, <q_j+1,p_k>,max",maxval(abs(qw_jk))
!            write(6,*) " "
!            write(6,*) qw_jk
!            write(6,*) " "

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
!                write(6,*) " "
!                write(6,*) " vector",i,"lost orth, pw_jk=",term,eps1
!                write(6,*) " "
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
!                write(6,*) " "
!                write(6,*) " vector",i,"lost orth, qw_jk=",term,eps1
!                write(6,*) " "
                qorthopt(i) = 1
              endif
            enddo
            qoterm = sum(qorthopt)
            write(6,*) "---# vecs that lost biortho",poterm,qoterm,"---"
            write(6,*) " "
!
            deallocate(porthopt,qorthopt)


               ! update orthfreq
             oterm = max(poterm,qoterm)  
             if (oterm .gt. 0) then
               if (prbo .eq. 2) then    
                 orthfreq = orthfreq + 2
               else
                 orthfreq = orthfreq + 1
               endif
             endif 

! Construct the Pj and Qj
            if ((poterm .ge. 1) .or. (qoterm .ge. 1)) then
               if (prbo .eq. 2) then    
                 iterm = 2*k + 2*(k+1)
               else
                 iterm = 2*k
               endif

               write(6,*) " doing a full 2-sided GS as biorthog failed"
               print *, " # rebiorthogonalizations = ",iterm
               print *
               
! Method 1 : Preferred/most robust way to do y to do Gram-Schmidt with (I-Q_k*(P_k)^H)
! Uses matrices and does the matrix-matrix-vector product
                    allocate(iden1(m,m),iden2(m,m),Pj(m,k),Qj(m,k))
                    iden1 = 0.0d0
                    iden2 = 0.0d0
                    Pj = 0.0d0
                    Qj = 0.0d0
                    Pj = Planc(:,1:k)
                    Qj = Qlanc(:,1:k)   
        !          write(6,*) "Pj=",Pj
        !          write(6,*) " "
        !          write(6,*) "Qj=",Qj
        !          write(6,*) " "  
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

              ! Set F_j to iter+1 to do the second full reorthog in the subsequent iteration
               F_j = k + 1
            else
               iterm = 0
            endif   !actual rebiorthog

            deallocate(pw_jk,qw_jk)
          endif  ! re biorthog condition

! 
          term1 = dot_product(qp,pp) ! following netlib Alg 7.13

!  A measure of the biorthogonality
          zborth(k) = term1

          term = dot_product(qp,qp) ! following netlib Alg 7.13
!          write(6,*) ""
!          write(6,*) "(r^_(j+1),r^_(j+1))=",term
!          write(6,*) "norm=",sqrt(term) 
!          write(6,*) ""
          term = dot_product(pp,pp) ! following netlib Alg 7.13
!          write(6,*) ""
!          write(6,*) "(s^_(j+1),s^_(j+1))=",term
!          write(6,*) "norm=",sqrt(term) 
!          write(6,*) ""
!          write(6,*) ""
!          write(6,*) "(r^_(j+1),s^_(j+1))=",term1 
!          write(6,*) ""
          lnorm1 = sqrt(abs(term1))  ! following yambo paper
          lnorm2 = conjg(term1)/lnorm1
!          write(6,*) "|(r^_(j+1),s^_(j+1))|=",lnorm1 
!          write(6,*) "(v^_(j+1),w^_(j+1))=",lnorm1
!          write(6,*) "lnorm2 =",lnorm2
!          write(6,*) ""
          zbeta(k+1)= lnorm1 ! following yambo paper
          zgamma(k+1)=lnorm2 ! following Y. Saad's paper
!            write(6,*) ""
!            write(6,*) "for iteration ",k,"beta =",lnorm
!            write(6,*) ""


!  Updating v_(j+1), to be used at start of (j+1) iteration, can use the tol here also
          if (k .ne. mlanc) then
            q=(1.0d0/lnorm1)*qp! v_(j+1)=v^_(j+1)/delta_(j+1)
            p=(1.0d0/conjg(lnorm2))*pp ! w_(j+1)=w^_(j+1)/beta_(j+1)

!  Add new vectors 
            Planc(:,k+1) = p
            Qlanc(:,k+1) = q

! Check for orthog after partial rebiorthog
            call zcheck_orthog_frob(Planc(:,1:k+1),Qlanc(:,1:k+1),& 
                                       m,k+1,frob)
            write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(6,*) " "
            write(6,*) "Norm after adding new Lanc vecs = ",frob
            write(6,*) " "
            write(6,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(6,*) " "
          endif

! update total orthogs and write to file
          
          orthtot = orthtot + iterm
          write(222,*)k,oterm,2*oterm,orthfreq,orthtot
        enddo ! END OF LANCZOS LOOP
!  Print out final total orthogonalizations
        write(6,*) " "
        write(6,*) "****** TOTAL PARTIAL ORTHOGS =",orthtot," VS FULL RE ORTHOGS = ",(mlanc+1)*mlanc,"******"
        write(6,*) " "
        write(6,*) " "
        return
        end subroutine lanczos_non_herm_prbo
