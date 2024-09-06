!subroutine act_1e( psi, eband, nk, nkx, nky, nkz, nb, nbc, nbv )
!  implicit none
!  !
!  integer :: nk, nkx, nky, nkz, nb, nbc, nbv
!  !
!  complex( kind = kind( 1.0d0 ) ) :: psi( nbv, nbc, nk )
!  real( kind = kind( 1.0d0 ) ) :: eband( nb, nk )
!  !
!  integer :: i, j, ik
!  real( kind = kind( 1.0d0 ) ) :: ec, ev
!  ! 
!  do i = 1, nbv
!     do j = 1, nbc
!        do ik = 1, nk
!            ec = eband( j + nbv, ik )
!            ev = eband( i, ik )
!            psi( i, j, ik ) = psi( i, j, ik ) * ( ec - ev )
!        end do
!     end do
!  end do
!  !
!  return
!end subroutine act_1e

subroutine act_1e( psi, eband, nk, nkx, nky, nkz, nb, nbc, nbv )
  implicit none
  !
  integer :: nk, nkx, nky, nkz, nb, nbc, nbv
  !
  complex( kind = kind( 1.0d0 ) ) :: psi( nbv, nbc, nk )
  real( kind = kind( 1.0d0 ) ) :: eband( nb, nk )
  !
  integer :: i, j, ik
  real( kind = kind( 1.0d0 ) ) :: ec, ev,rtemp,rtemp1,scal
  complex( kind = kind( 1.0d0 ) ) :: rm1,ztemp
  ! 
  rm1 = -1
  rm1 = sqrt(rm1)

  ! Get value of scal for complex energies
  open(unit=99, file='scal_ze', form='formatted', status='unknown')
  rewind 99
  read(99,*) scal
  print *
  print *, "scal = ",scal
  print *
  close(unit=99)

  do i = 1, nbv
     do j = 1, nbc
        do ik = 1, nk
            ec = eband( j + nbv, ik )
            ev = eband( i, ik )
            rtemp1 = ec - ev
            ztemp = rtemp1 - scal*rtemp1*rm1

!            psi( i, j, ik ) = psi( i, j, ik ) * ( ec - ev )
            psi( i, j, ik ) = psi( i, j, ik ) * ztemp
        end do
     end do
  end do
  !
  return
end subroutine act_1e

!subroutine act_C( psi, sfact, nk, nkx, nky, nkz, nb, nbc, nbv )
!  ! For generalized Lanczos algorithm, eigenvalue problem, C*x = lambda*B*x
!  ! With H_BSE = A = C*B , C is given by act_C. It will give the diagonal blocks in the 
!  ! full BSE transition space Hamiltonian like the single particle energy difference operator. 
!  implicit none
!  !
!  integer :: nk, nkx, nky, nkz, nb, nbc, nbv
!  !
!  complex( kind = kind( 1.0d0 ) ) :: psi( nbv, nbc, nk ), sfact(nbv,nbc,nk)
!  real( kind = kind( 1.0d0 ) ) :: eband( nb, nk )
!  !
!  integer :: i, j, ik
!  real( kind = kind( 1.0d0 ) ) :: ec, ev
!  !
!  psi(:,:,:) = psi(:,:,:)*sfact(:,:,:)
!!  do i = 1, nbv
!!     do j = 1, nbc
!!        do ik = 1, nk
!!!            ec = eband( j + nbv, ik )
!!!            ev = eband( i, ik )
!!            psi( i, j, ik ) = psi( i, j, ik ) ! * sfact(i,j,ik)
!!        end do
!!     end do
!!  end do
!  !
!  return
!end subroutine act_C

subroutine create_iden( sfact, eband, nk, nkx, nky, nkz, nb, nbc, nbv )
  ! For generalized Lanczos algorithm, eigenvalue problem, C*x = lambda*B*x
  ! With H_BSE = A = C*B , C is given by act_1e_iden. It will give the diagonal blocks in the 
  ! full BSE transition space Hamiltonian like the single particle energy difference operator. 
  ! By scaling each element with the energy difference, in the b=0,l=0 case, it will give the identity
  ! submatrices once we multiply by sfact.

  implicit none
  !
  integer :: nk, nkx, nky, nkz, nb, nbc, nbv
  !
  complex( kind = kind( 1.0d0 ) ) :: psi( nbv, nbc, nk )
  real( kind = kind( 1.0d0 ) ) :: eband( nb, nk ), sfact( nbv, nbc, nk )
  !
  integer :: i, j, ik
  real( kind = kind( 1.0d0 ) ) :: ec, ev, rterm, iden(nbv,nbc,nk)
  ! 
  iden = 0.0d0
  open(unit=99, file='iden_matrix', form='formatted', status='unknown')
  rewind 99

  do i = 1, nbv
     do j = 1, nbc
        do ik = 1, nk
            ec = eband( j + nbv, ik )
            ev = eband( i, ik )
            rterm = 1.0d0/( ec -ev )
            iden( i, j, ik ) = rterm * sfact(i, j, ik )
            write(99,*) iden( i, j, ik )
        end do
     end do
  end do
  !
  close(unit=99)
  return
end subroutine create_iden

subroutine read_iden( sfact, eband, nk, nkx, nky, nkz, nb, nbc, nbv, iden )
  ! For generalized Lanczos algorithm, eigenvalue problem, C*x = lambda*B*x
  ! With H_BSE = A = C*B , C is given by act_1e_iden. It will give the diagonal blocks in the 
  ! full BSE transition space Hamiltonian like the single particle energy difference operator. 
  ! By scaling each element with the energy difference, in the b=0,l=0 case, it will give the identity
  ! submatrices once we multiply by sfact.

  ! Reads this identity matrix from the saved file to an array

  implicit none
  !
  integer :: nk, nkx, nky, nkz, nb, nbc, nbv
  !
  complex( kind = kind( 1.0d0 ) ) :: psi( nbv, nbc, nk )
  real( kind = kind( 1.0d0 ) ) :: eband( nb, nk ), sfact( nbv, nbc, nk )
  !
  integer :: i, j, ik
  real( kind = kind( 1.0d0 ) ) :: ec, ev, rterm, iden(nbv,nbc,nk)
  ! 
  iden = 0.0d0
  open(unit=99, file='iden_matrix', form='formatted', status='unknown')
  rewind 99

  do i = 1, nbv
     do j = 1, nbc
        do ik = 1, nk
!            ec = eband( j + nbv, ik )
!            ev = eband( i, ik )
!            rterm = 1.0d0/( ec -ev )
!            iden( i, j, ik ) = rterm * sfact(i, j, ik )
            read(99,*) iden( i, j, ik )
        end do
     end do
  end do
  !
  close(unit=99)
  return
end subroutine read_iden

subroutine act_cmplx_1e( psi, eband, nk, nkx, nky, nkz, nb, nbc, nbv )
! This routine includes complex energies to the TDA Hamiltonian.
! For now, random numbers will be used to determine the complex energies, which
! will capture the effect of phonons as in the paper "Ab initio finite temperature excitons".
! In the future, one will look into using DFPT to get the accurate phonon shifts
  implicit none
  !
  integer :: nk, nkx, nky, nkz, nb, nbc, nbv
  !
  complex( kind = kind( 1.0d0 ) ) :: psi( nbv, nbc, nk )
  real( kind = kind( 1.0d0 ) ) :: eband( nb, nk )
  !
  integer :: i, j, ik
  real( kind = kind( 1.0d0 ) ) :: ec, ev

  ! Generating random numbers
  integer, parameter :: sedsz=33
  integer :: seed(sedsz)
  real( kind = kind( 1.0d0 ) ) :: rand_e( nb, nk ),rand_ec,rand_ev,rtemp,rtemp1,scal
  complex( kind = kind( 1.0d0 ) ) :: rm1, ztemp
  ! 
  rm1 = -1
  rm1 = sqrt(rm1)

  ! Generating the array of random energies
  open( unit=99, file='scal_ze', form='formatted', status='unknown' )
  rewind 99
  read ( 99,* ) scal
  close( unit=99 )
!  scal = 1.0E-2
  rand_e = 0.0
  call random_seed(put=seed(1:sedsz))
  call random_number(rand_e)
!  write(6,*) " "
!  write(6,*) "Writing out the random numbers for k = 1"
!  write(6,*) " "
!  write(6,*) rand_e(:,1)
!  write(6,*) " "

  do i = 1, nbv
     do j = 1, nbc
        do ik = 1, nk
            ec = eband( j + nbv, ik )
            ev = eband( i, ik )
            rtemp1 = ec - ev
            rand_ec = rand_e( j + nbv, ik )
            rand_ev = rand_e( i, ik )
            rtemp = rand_ec - rand_ev
            ztemp = rtemp1 - scal*rtemp1*rm1
!            ztemp = rtemp1 + rtemp*rm1
!            psi( i, j, ik ) = psi( i, j, ik ) * ( ec - ev )
            psi( i, j, ik ) = psi( i, j, ik ) * ztemp
        end do
     end do
  end do
  !
  return
end subroutine act_cmplx_1e

subroutine act_cmplx_dag_1e( psi, eband, nk, nkx, nky, nkz, nb, nbc, nbv )
! This routine includes complex energies to the TDA Hamiltonian.
! For now, random numbers will be used to determine the complex energies, which
! will capture the effect of phonons as in the paper "Ab initio finite temperature excitons".
! In the future, one will look into using DFPT to get the accurate phonon shifts
! Does the Hermitian conjugate/dagger

  implicit none
  !
  integer :: nk, nkx, nky, nkz, nb, nbc, nbv
  !
  complex( kind = kind( 1.0d0 ) ) :: psi( nbv, nbc, nk )
  real( kind = kind( 1.0d0 ) ) :: eband( nb, nk )
  !
  integer :: i, j, ik
  real( kind = kind( 1.0d0 ) ) :: ec, ev

  ! Generating random numbers
  integer, parameter :: sedsz=33
  integer :: seed(sedsz)
  real( kind = kind( 1.0d0 ) ) :: rand_e( nb, nk ),rand_ec,rand_ev,rtemp,rtemp1,scal
  complex( kind = kind( 1.0d0 ) ) :: rm1,ztemp
  ! 
  rm1 = -1
  rm1 = sqrt(rm1)

  ! Generating the array of random energies

  open( unit=99, file='scal_ze', form='formatted', status='unknown' )
  rewind 99
  read ( 99,* ) scal
  close( unit=99 )
!  scal = 1.0E-1
  rand_e = 0.0
  call random_seed(put=seed(1:sedsz))
  call random_number(rand_e)
!  write(6,*) " "
!  write(6,*) "Writing out the random numbers for k = 1"
!  write(6,*) " "
!  write(6,*) rand_e(:,1)
!  write(6,*) " "

  do i = 1, nbv
     do j = 1, nbc
        do ik = 1, nk
            ec = eband( j + nbv, ik )
            ev = eband( i, ik )
            rand_ec = rand_e( j + nbv, ik )
            rand_ev = rand_e( i, ik )
            rtemp = rand_ec - rand_ev
            rtemp1 = ec - ev
            ztemp = rtemp1 - scal*rtemp1*rm1
!            ztemp = rtemp1 + rtemp*rm1
!            psi( i, j, ik ) = psi( i, j, ik ) * ( ec - ev )
            psi( i, j, ik ) = psi( i, j, ik ) * conjg(ztemp)
        end do
     end do
  end do
  !
  return
end subroutine act_cmplx_dag_1e
