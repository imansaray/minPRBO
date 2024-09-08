subroutine act( psi, u, eband, kk, nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv, xk, nkret, kret, mr, &
                lflag, bflag, backf, aldaf, allow, sfact, bande )
  implicit none
  !
  integer :: nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv
  integer :: lflag, bflag, backf, aldaf, nkret, bande
  !
  integer :: kk( nk, 3 ), kret( nkret )
  complex( kind = kind( 1.0d0 ) ) :: psi( nbv, nbc, nk ), u( nb, nk, ng )
  real( kind = kind( 1.0d0 ) ) :: eband( nb, nk ), xk( nk, 3 )
  real( kind = kind( 1.0d0 ) ) :: allow( nbv, nbc, nk ), sfact( nbv, nbc, nk )
  real( kind = kind( 1.0d0 ) ) :: mr( nkret, ng, ng )
  !
  complex( kind = kind( 1.0d0 ) ), allocatable :: c_l(:,:,:), c_b(:,:,:), ps2(:,:,:)
  complex( kind = kind( 1.0d0 ) ), allocatable :: c_a(:,:,:)
  !
  real( kind = kind( 1.0d0 ) ) :: zr, zi
  complex( kind = kind( 1.0d0 ) ) :: z, rm1
  complex( kind = kind( 1.0d0 ) ), external :: wvfdot
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  allocate( c_l( nbv, nbc, nk ), c_b( nbv, nbc, nk ), c_a( nbv, nbc, nk ) )
  allocate( ps2( nbv, nbc, nk ) )
  !
!  psi = 1.0d0 !IM
  psi( :, :, : ) = psi( :, :, : ) * allow( :, :, : )
  ps2( :, :, : ) = psi( :, :, : )
  c_l( :, :, : ) = psi( :, :, : )
  c_b( :, :, : ) = psi( :, :, : )
  c_a( :, :, : ) = psi( :, :, : )
  !
  write(6,*) 'begin act'
  if ( bande .eq. 1 ) then
     call act_1e( psi, eband, nk, nkx, nky, nkz, nb, nbc, nbv )
     if ( backf .eq. 1 ) call act_1e( psi, eband, nk, nkx, nky, nkz, nb, nbc, nbv )
     psi( :, :, : ) = psi( :, :, : ) * sfact( :, :, : )
  else
     psi = 0
     write ( 6, * ) 'bande = 0'
  end if
  !
  z = wvfdot( ps2, psi, allow, nk * nbc * nbv )
  zr = z
  zi = -rm1 * z
  write ( 6, '(1a6,2x,2(1x,1e22.15))' ) 'one-el', zr, zi
  !
  if ( bflag .eq. 1 ) then
     if ( backf .eq. 1 ) call bfnorm( c_b, eband, nk, nb, nbc, nbv )
     call smbubble( c_b, u, nk, ng, ngx, ngy, ngz, nb, nbc, nbv )
     if ( backf .eq. 1 ) call bfnorm( c_b, eband, nk, nb, nbc, nbv )
     z = wvfdot( ps2, c_b, allow, nk * nbc * nbv )
     zr = z
     zi = -rm1 * z
     write ( 6, '(1a6,2x,2(1x,1e22.15))' ) 'bubble', 2 * zr, 2 * zi
     psi( :, :, : ) = psi( :, :, : ) + 2.0d0 * c_b( :, :, : )
  end if
  deallocate( c_b )
  !
  write(6,*) 'p. 1'
  if ( aldaf .eq. 1 ) then
     if ( backf .eq. 1 ) call bfnorm( c_a, eband, nk, nb, nbc, nbv )
     call smkxc( c_a, u, nk, ng, ngx, ngy, ngz, nb, nbc, nbv )
     if ( backf .eq. 1 ) call bfnorm( c_a, eband, nk, nb, nbc, nbv )
     z = wvfdot( ps2, c_a, allow, nk * nbc * nbv )
     zr = z
     zi = -rm1 * z
     write ( 6, '(1a6,2x,2(1x,1e22.15))' ) 'ldakxc', 2 * zr, 2 * zi
     psi( :, :, : ) = psi( :, :, : ) + 2.0d0 * c_a( :, :, : )
  end if
  deallocate( c_a )
  !
  if ( lflag .eq. 1 ) then
  write(6,*) 'p. 2'
     call smlsauber( c_l, u, kk, nk, nkx, nky, nkz, ng, nb, nbc, nbv, nkret, kret, mr )
  write(6,*) 'p. 3'
     z = wvfdot( ps2, c_l, allow, nk * nbc * nbv )
     zr = z
     zi = -rm1 * z
     write ( 6, '(1a6,2x,2(1x,1e22.15))' ) 'ladder', -zr, -zi
     psi( :, :, : ) = psi( :, :, : ) - c_l( :, :, : )
  end if

  deallocate( c_l, ps2 )
   write(6,*) 'end act' !
  return
end subroutine act


subroutine act_cmplx( psi, u, eband, kk, nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv, xk, nkret, kret, mr, &
                lflag, bflag, backf, aldaf, allow, sfact, bande )

  ! For now, this only does the action of the complex energies in act_cmplx_1e.      
  implicit none
  !
  integer :: nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv
  integer :: lflag, bflag, backf, aldaf, nkret, bande
  !
  integer :: kk( nk, 3 ), kret( nkret )
  complex( kind = kind( 1.0d0 ) ) :: psi( nbv, nbc, nk ), u( nb, nk, ng )
  real( kind = kind( 1.0d0 ) ) :: eband( nb, nk ), xk( nk, 3 )
  real( kind = kind( 1.0d0 ) ) :: allow( nbv, nbc, nk ), sfact( nbv, nbc, nk )
  real( kind = kind( 1.0d0 ) ) :: mr( nkret, ng, ng )
  !
  complex( kind = kind( 1.0d0 ) ), allocatable :: c_l(:,:,:), c_b(:,:,:), ps2(:,:,:)
  complex( kind = kind( 1.0d0 ) ), allocatable :: c_a(:,:,:)
  !
  real( kind = kind( 1.0d0 ) ) :: zr, zi
  complex( kind = kind( 1.0d0 ) ) :: z, rm1
  complex( kind = kind( 1.0d0 ) ), external :: wvfdot
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  allocate( c_l( nbv, nbc, nk ), c_b( nbv, nbc, nk ), c_a( nbv, nbc, nk ) )
  allocate( ps2( nbv, nbc, nk ) )
  !
  psi( :, :, : ) = psi( :, :, : ) * allow( :, :, : )
  ps2( :, :, : ) = psi( :, :, : )
  c_l( :, :, : ) = psi( :, :, : )
  c_b( :, :, : ) = psi( :, :, : )
  c_a( :, :, : ) = psi( :, :, : )
  !
  write(6,*) 'begin act'
  if ( bande .eq. 1 ) then
     call act_cmplx_1e( psi, eband, nk, nkx, nky, nkz, nb, nbc, nbv )
     if ( backf .eq. 1 ) call act_cmplx_1e( psi, eband, nk, nkx, nky, nkz, nb, nbc, nbv )
     psi( :, :, : ) = psi( :, :, : ) * sfact( :, :, : )
  else
     psi = 0
     write ( 6, * ) 'bande = 0'
  end if
  !
  z = wvfdot( ps2, psi, allow, nk * nbc * nbv )
  zr = z
  zi = -rm1 * z
  write ( 6, '(1a6,2x,2(1x,1e22.15))' ) 'one-el', zr, zi
  !
  if ( bflag .eq. 1 ) then
     if ( backf .eq. 1 ) call bfnorm( c_b, eband, nk, nb, nbc, nbv )
     call smbubble( c_b, u, nk, ng, ngx, ngy, ngz, nb, nbc, nbv )
     if ( backf .eq. 1 ) call bfnorm( c_b, eband, nk, nb, nbc, nbv )
     z = wvfdot( ps2, c_b, allow, nk * nbc * nbv )
     zr = z
     zi = -rm1 * z
     write ( 6, '(1a6,2x,2(1x,1e22.15))' ) 'bubble', 2 * zr, 2 * zi
     psi( :, :, : ) = psi( :, :, : ) + 2.0d0 * c_b( :, :, : )
  end if
  deallocate( c_b )
  !
  write(6,*) 'p. 1'
  if ( aldaf .eq. 1 ) then
     if ( backf .eq. 1 ) call bfnorm( c_a, eband, nk, nb, nbc, nbv )
     call smkxc( c_a, u, nk, ng, ngx, ngy, ngz, nb, nbc, nbv )
     if ( backf .eq. 1 ) call bfnorm( c_a, eband, nk, nb, nbc, nbv )
     z = wvfdot( ps2, c_a, allow, nk * nbc * nbv )
     zr = z
     zi = -rm1 * z
     write ( 6, '(1a6,2x,2(1x,1e22.15))' ) 'ldakxc', 2 * zr, 2 * zi
     psi( :, :, : ) = psi( :, :, : ) + 2.0d0 * c_a( :, :, : )
  end if
  deallocate( c_a )
  !
  if ( lflag .eq. 1 ) then
  write(6,*) 'p. 2'
     call smlsauber( c_l, u, kk, nk, nkx, nky, nkz, ng, nb, nbc, nbv, nkret, kret, mr )
  write(6,*) 'p. 3'
     z = wvfdot( ps2, c_l, allow, nk * nbc * nbv )
     zr = z
     zi = -rm1 * z
     write ( 6, '(1a6,2x,2(1x,1e22.15))' ) 'ladder', -zr, -zi
     psi( :, :, : ) = psi( :, :, : ) - c_l( :, :, : )
  end if

  deallocate( c_l, ps2 )
   write(6,*) 'end act' !
  return
end subroutine act_cmplx

subroutine act_cmplx_dag( psi, u, eband, kk, nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv, xk, nkret, kret, mr, &
                lflag, bflag, backf, aldaf, allow, sfact, bande )

  ! For now, this only does the action of the complex energies in act_cmplx_dag_1e.      
  implicit none
  !
  integer :: nk, nkx, nky, nkz, ng, ngx, ngy, ngz, nb, nbc, nbv
  integer :: lflag, bflag, backf, aldaf, nkret, bande
  !
  integer :: kk( nk, 3 ), kret( nkret )
  complex( kind = kind( 1.0d0 ) ) :: psi( nbv, nbc, nk ), u( nb, nk, ng )
  real( kind = kind( 1.0d0 ) ) :: eband( nb, nk ), xk( nk, 3 )
  real( kind = kind( 1.0d0 ) ) :: allow( nbv, nbc, nk ), sfact( nbv, nbc, nk )
  real( kind = kind( 1.0d0 ) ) :: mr( nkret, ng, ng )
  !
  complex( kind = kind( 1.0d0 ) ), allocatable :: c_l(:,:,:), c_b(:,:,:), ps2(:,:,:)
  complex( kind = kind( 1.0d0 ) ), allocatable :: c_a(:,:,:)
  !
  real( kind = kind( 1.0d0 ) ) :: zr, zi
  complex( kind = kind( 1.0d0 ) ) :: z, rm1
  complex( kind = kind( 1.0d0 ) ), external :: wvfdot
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  allocate( c_l( nbv, nbc, nk ), c_b( nbv, nbc, nk ), c_a( nbv, nbc, nk ) )
  allocate( ps2( nbv, nbc, nk ) )
  !
  psi( :, :, : ) = psi( :, :, : ) * allow( :, :, : )
  ps2( :, :, : ) = psi( :, :, : )
  c_l( :, :, : ) = psi( :, :, : )
  c_b( :, :, : ) = psi( :, :, : )
  c_a( :, :, : ) = psi( :, :, : )
  !
  write(6,*) 'begin act'
  if ( bande .eq. 1 ) then
     call act_cmplx_dag_1e( psi, eband, nk, nkx, nky, nkz, nb, nbc, nbv )
     if ( backf .eq. 1 ) call act_cmplx_dag_1e( psi, eband, nk, nkx, nky, nkz, nb, nbc, nbv )
     psi( :, :, : ) = psi( :, :, : ) * sfact( :, :, : )
  else
     psi = 0
     write ( 6, * ) 'bande = 0'
  end if
  !
  z = wvfdot( ps2, psi, allow, nk * nbc * nbv )
  zr = z
  zi = -rm1 * z
  write ( 6, '(1a6,2x,2(1x,1e22.15))' ) 'one-el', zr, zi
  !
  if ( bflag .eq. 1 ) then
     if ( backf .eq. 1 ) call bfnorm( c_b, eband, nk, nb, nbc, nbv )
     call smbubble( c_b, u, nk, ng, ngx, ngy, ngz, nb, nbc, nbv )
     if ( backf .eq. 1 ) call bfnorm( c_b, eband, nk, nb, nbc, nbv )
     z = wvfdot( ps2, c_b, allow, nk * nbc * nbv )
     zr = z
     zi = -rm1 * z
     write ( 6, '(1a6,2x,2(1x,1e22.15))' ) 'bubble', 2 * zr, 2 * zi
     psi( :, :, : ) = psi( :, :, : ) + 2.0d0 * c_b( :, :, : )
  end if
  deallocate( c_b )
  !
  write(6,*) 'p. 1'
  if ( aldaf .eq. 1 ) then
     if ( backf .eq. 1 ) call bfnorm( c_a, eband, nk, nb, nbc, nbv )
     call smkxc( c_a, u, nk, ng, ngx, ngy, ngz, nb, nbc, nbv )
     if ( backf .eq. 1 ) call bfnorm( c_a, eband, nk, nb, nbc, nbv )
     z = wvfdot( ps2, c_a, allow, nk * nbc * nbv )
     zr = z
     zi = -rm1 * z
     write ( 6, '(1a6,2x,2(1x,1e22.15))' ) 'ldakxc', 2 * zr, 2 * zi
     psi( :, :, : ) = psi( :, :, : ) + 2.0d0 * c_a( :, :, : )
  end if
  deallocate( c_a )
  !
  if ( lflag .eq. 1 ) then
  write(6,*) 'p. 2'
     call smlsauber( c_l, u, kk, nk, nkx, nky, nkz, ng, nb, nbc, nbv, nkret, kret, mr )
  write(6,*) 'p. 3'
     z = wvfdot( ps2, c_l, allow, nk * nbc * nbv )
     zr = z
     zi = -rm1 * z
     write ( 6, '(1a6,2x,2(1x,1e22.15))' ) 'ladder', -zr, -zi
     psi( :, :, : ) = psi( :, :, : ) - c_l( :, :, : )
  end if

  deallocate( c_l, ps2 )
   write(6,*) 'end act' !
  return
end subroutine act_cmplx_dag
