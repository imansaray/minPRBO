
!/local/hadley/cvs_ai2nbse/AI2NBSE/Src/NBSE/zbridge/zf77/spectrum.f
! This program takes in the Lanczos coeffs a, b, c, where these coeffs
! are assumed complex numbers from a.dat, b.dat, c.dat

      program spectrum
      implicit none
      integer, parameter :: stdin = 5
      integer, parameter :: eps1 = 41, eps2 = 42
      integer, parameter :: loss = 43, refl = 44, inds = 45
      integer, parameter :: adat = 55, bdat = 56, cdat = 57
      character * 9, parameter :: f9 = 'formatted'
      character * 7, parameter :: u7 = 'unknown'
      integer i, n, ne, ie, nk, idum, nonherm
      double precision gam, el, eh, vol, pnorm, fact, de, ere, lossf
      double precision reeps, imeps, indref, indabs, rad
      double precision ref, theta, pi, term1,term2,term3
      double complex e, arg, r, rm, rp, al, be, eps
!      double precision, allocatable :: a(:), b(:)
      double complex, allocatable :: a(:), b(:), c(:)

      !mpp August 2008
      integer, parameter :: opcons = 46
      double complex refrac
      double precision eloss, reflct,mu, omega, hart, bohr, ryd
      character * 512 slog
      parameter ( ryd  = 13.605 698d0) ! Rydberg to eV
      parameter (hart = 2 * ryd)       ! Hartree to Rydberg
      double precision, parameter :: ha2ev = 27.211386245988 ! Hartree to eV
      integer :: niter, bflag, lflag, backf, aldaf, qpflg, bwflg,bande
      !
      pi = 4.0d0 * atan( 1.0d0 )
      open( unit=99, file='omega.h', form=f9, status=u7 )
      call rgetval( vol, 'volbo' )
      close( unit=99 )
      open( unit=99, file='ndone.h', form=f9, status=u7 )
      call igetval( niter, 'ndone' )
      close( unit=99 )
C Modified by fer
C     read ( stdin, * ) gam, el, eh, ne
C     read ( stdin, * ) n
      open( unit=99, file='spect.ipt', form=f9, status=u7 )
      read ( 99, * ) gam, el, eh, ne, n
      write(6,*) " "
      write(6,*) " gamma = ",gam," eV"
      write(6,*) " min energy = ",el, "eV"
      write(6,*) " max energy  = ",eh," eV"
      write(6,*) " # of energy points = ",ne
      write(6,*) " # of Lanczos iterations = ",n
      write(6,*) " "
      close( unit=99 )

      open( unit=99, file='niter.h', form='formatted', status='unknown')
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

c      niter = n
c      if ( ( n .lt. 0 ) .or. ( n .gt. niter ) ) n = niter
      allocate( a( 0 : niter ), b( 0 : niter + 1 ), c(0:niter+1) )
      write(6,*) " "
      print *, " nonherm =",nonherm, "niter =",niter
      write(6,*) " Reading in the Lanczos coeficients a, b, possibly c"
      write(6,*) " "
      a = 0
      b = 0
      c = 0

      open( unit=adat, file='a.dat', form=f9, status=u7 )
      rewind adat
      read ( adat, * ) pnorm, nk
      do i = 0, n
        read( adat, * ) idum, a( i )
      end do
      close( unit=adat )

      open( unit=bdat, file='b.dat', form=f9, status=u7 )
      rewind bdat
      read ( bdat, * ) pnorm, nk
      do i = 0, n + 1
        read( bdat, * ) idum, b( i )
      end do
      close( unit=bdat )
     
c      nonherm = 0 
      if (nonherm .eq. 1) then
         write(6,*) " Non Hermitian case, reading c.dat "
         write(6,*) " "
         open( unit=cdat, file='c.dat', form=f9, status=u7 )
         read ( cdat, * ) pnorm, nk
         do i = 0, n + 1
           read( cdat, * ) idum, c( i )
         end do
         close( unit=cdat )
      end if

!      fact = 8.0d0 * pi ** 2 * 27.2114d0 * pnorm ** 2 /
!     &       ( dble( nk ) * vol ) / pi

      fact = 8.0d0 * pi ** 2 * ha2ev * pnorm ** 2 /
     &       ( dble( nk ) * vol ) / pi
      open( unit=eps1, file='eps1', form=f9, status=u7 )
      open( unit=eps2, file='eps2', form=f9, status=u7 )
      open( unit=loss, file='loss', form=f9, status=u7 )
      open( unit=refl, file='refl', form=f9, status=u7 )
      open( unit=inds, file='inds', form=f9, status=u7 )

      !mpp Aug. 2008
      !open opcons.dat and write a short header
      open( unit=opcons, file='opcons_nherm.dat', form=f9, status=u7 )
      slog="#   omega (eV)      epsilon_1       epsilon_2       n"//
     &"               kappa           mu (cm^(-1))    R"//
     &"               epsinv"
      write(opcons,fmt="(a)") slog(1:125)


      rewind eps1
      rewind eps2
      rewind loss
      rewind refl
      rewind inds

      ! for non hermitian case
      if (nonherm .eq. 0) then
        c = b
      end if

!  With complex energies, set gamma = 0.0, switched c and b, since I
!  switched them in the Lanczos routines
!      gam = 0.0      
      de = ( eh - el ) / dble( ne )
      do ie = 0, ne
        ere = el + de * dble( ie )
        e = dcmplx( ere, gam )
        if (bwflg .eq. 0) then
!           arg = ( ere - real(a( n )) )**2-4.d0*real(b( n + 1 ))** 2
!           arg = cdsqrt( arg )
!           rp = 0.5d0 * ( ere - real(a( n )) + arg )
!           rm = 0.5d0 * ( ere - real(a( n )) - arg )
           arg = ( ere - a( n ) )**2 - 4.d0*c( n + 1 )*conjg(b(n+1))
           arg = cdsqrt( arg )
           rp = 0.5d0 * ( ere - a( n ) + arg )
           rm = 0.5d0 * ( ere - a( n ) - arg )
        else
           arg = ( ere - a( n ) )**2 - 4.d0*c( n + 1 )*conjg(b(n+1))
           arg = cdsqrt( arg )
           rp = 0.5d0 * ( ere - a( n ) + arg )
           rm = 0.5d0 * ( ere - a( n ) - arg )

        endif
        if ( dimag( rp ) .lt. 0.d0 ) then
          r = rp
        else
          r = rm
        end if
        al =    e - a( n ) - r
        be  = - e - a( n ) - r
        do i = n, 1, - 1 ! initially from 0, n
          al =    e - a( i ) - c( i + 1 )*conjg(b(i+1)) / al
          be  = - e - a( i ) - c( i + 1 )*conjg(b(i+1)) / be
        end do
        if (bwflg .eq. 0) then
           eps = 1.d0 - fact / al  !- fact / be ! Leave the negative enrgy part out for comparison to GMRES pos energy soln
        else
           eps = 1.d0 - fact / al
        endif
        reeps = dble( eps )
        imeps = dimag( eps )
        rad = sqrt( reeps ** 2 + imeps ** 2 )
        theta = acos( reeps / rad ) / 2
        indref = sqrt( rad ) * cos( theta )
        indabs = sqrt( rad ) * sin( theta )
        ref = ( ( indref - 1 ) ** 2 + indabs ** 2 ) /
     &        ( ( indref + 1 ) ** 2 + indabs ** 2 )
        lossf = imeps / ( reeps ** 2 + imeps ** 2 )
        write ( eps1, '(2x,1f10.5,1f30.20)' ) ere, reeps
        write ( eps2, '(2x,1f10.5,1f30.20)' ) ere, imeps
        write ( loss, * ) ere, lossf
        write ( inds, '(5(2x,1e15.8))' ) ere, indref, indabs,
     &                                   indref ** 2 - indabs ** 2,
     &                                   2 * indref * indabs
        write ( refl, * ) ere, ref

        !mpp August 2008
        omega=ere/hart
        call eps2opt(eps-1.0,omega,refrac,mu,reflct,eloss)
        write (unit=opcons,fmt="(19e16.6)")
     &   ere,eps,refrac-1.0,mu,reflct,eloss

      end do
      close( unit=eps1 )
      close( unit=eps2 )
      close( unit=inds )
      close( unit=loss )
      close( unit=refl )
      !mpp August 2008
      close( unit=opcons )
!
      write(6,*) " "
      write(6,*) " All done with calculating the spectrum. "
      end


      subroutine eps2opt(eps,omega,refrac,mu,reflct,eloss)
      !given complex dielectric contsant minus 1 in eps and frequency in
      !omega, computes the complex index of refraction (refrac), normal
      !incidence reflectivity (reflct), absorption coefficient in
      !inverse angstroms (mu), and energy loss function (eloss).
      implicit none
      double complex refrac,eps
      double precision eloss, reflct,mu, omega
      double precision alpinv ,alphfs , bohr, ryd
      parameter (alpinv = 137.035 989 56d0)
      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
      !index of refraction N=n+ik=(epsilon)^(1/2)
      refrac=sqrt(eps+1)
      !normal incidence reflectance
      reflct=abs((refrac-1)/(refrac+1))**2
      !absorption coefficient in inverse angstroms
      mu=2*omega*alphfs*aimag(refrac)/bohr*1000
      eloss=-1.0*aimag((eps+1)**(-1))

      end !subroutine eps2opt()


