
!/local/hadley/cvs_ai2nbse/AI2NBSE/Src/NBSE/zbridge/zf77/spectrum.f

      program spectrum
      implicit none
      integer, parameter :: stdin = 5
      integer, parameter :: eps1 = 41, eps2 = 42
      integer, parameter :: loss = 43, refl = 44, inds = 45
      integer, parameter :: adat = 55, bdat = 56
      character * 9, parameter :: f9 = 'formatted'
      character * 7, parameter :: u7 = 'unknown'
      integer niter, i, n, ne, ie, nk, idum,iterm
      double precision gam, el, eh, vol, pnorm, fact, de, ere, lossf
      double precision reeps, imeps, indref, indabs, rad
      double precision ref, theta, pi
      double complex e, arg, r, rm, rp, al, be, eps
      double precision, allocatable :: a(:), b(:)

      !mpp August 2008
      integer, parameter :: opcons = 46
      double complex refrac
      double precision eloss, reflct,mu, omega, hart, bohr, ryd
      character * 512 slog
      parameter ( ryd  = 13.605 698d0) ! Rydberg to eV
      parameter (hart = 2 * ryd)       ! Hartree to Ry
      double precision, parameter :: ha2ev = 27.211386245988 ! Hartree to ev

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
c      read ( stdin, * ) gam, el, eh, ne, n
      
      open( unit=99,file='spect.ipt',form=f9, status=u7 )
      read ( 99, * ) gam, el, eh, ne, iterm
      write(6,*) " "
      write(6,*) " gamma = ",gam," eV "
      write(6,*) " min energy = ",el," eV "
      write(6,*) " max energy = ",eh," eV "
      write(6,*) " # of energy points = ",ne
      write(6,*) " # of Lanczos iterations = ",niter
      write(6,*) " "
      close( unit=99)

      n = niter
c      if ( ( n .lt. 0 ) .or. ( n .gt. niter ) ) n = niter
      allocate( a( 0 : niter ), b( 0 : niter + 1 ) )
      write(6,*) " "
      write(6,*) " Reading in the Lanczos coefficients a, b"
      write(6,*) " "

      open( unit=adat, file='a.dat', form=f9, status=u7 )
      rewind adat
      read ( adat, * ) pnorm, nk
      do i = 0, n
        read( adat, * ) idum, a( i )
      end do
      close( unit=adat )
      open( unit=bdat, file='b.dat', form=f9, status=u7 )
      read ( bdat, * ) pnorm, nk
      do i = 0, n + 1
        read( bdat, * ) idum, b( i )
      end do
      close( unit=bdat )
c      fact = 8.0d0 * pi ** 2 * 27.2114d0 * pnorm ** 2 /
c     &       ( dble( nk ) * vol ) / pi


      fact = 8.0d0 * pi ** 2 * ha2ev * pnorm ** 2 /
     &       ( dble( nk ) * vol ) / pi

      write(6,*) " "
      write(6,*) "fact= ",fact
      write(6,*) "fact(Ha) = ",fact/hart,"vol =",vol,"pnorm=",pnorm
      write(6,*) " "
      open( unit=eps1, file='eps1', form=f9, status=u7 )
      open( unit=eps2, file='eps2', form=f9, status=u7 )
      open( unit=loss, file='loss', form=f9, status=u7 )
      open( unit=refl, file='refl', form=f9, status=u7 )
      open( unit=inds, file='inds', form=f9, status=u7 )

      !mpp Aug. 2008
      !open opcons.dat and write a short header
      open( unit=opcons, file='opcons_herm.dat', form=f9, status=u7 )
      slog="#   omega (eV)      epsilon_1       epsilon_2       n"//
     &"               kappa           mu (cm^(-1))    R"//
     &"               epsinv"
      write(opcons,fmt="(a)") slog(1:125)


      rewind eps1
      rewind eps2
      rewind loss
      rewind refl
      rewind inds
      de = ( eh - el ) / dble( ne )
      do ie = 0, ne
        ere = el + de * dble( ie )
        e = dcmplx( ere, gam )
        arg = ( ere - a( n ) ) ** 2 - 4.d0 * b( n + 1 ) ** 2
        arg = cdsqrt( arg )
        rp = 0.5d0 * ( ere - a( n ) + arg )
        rm = 0.5d0 * ( ere - a( n ) - arg )
        if ( dimag( rp ) .lt. 0.d0 ) then
          r = rp
        else
          r = rm
        end if
        al =    e - a( n ) - r
        be  = - e - a( n ) - r
        do i = n, 0, - 1
c im        do i = n-0, 1, - 1
          al =    e - a( i ) - b( i + 1 ) ** 2 / al
          be  = - e - a( i ) - b( i + 1 ) ** 2 / be
        end do
c        print *
c        print *," iter =",ie,"energy =",ere
c        print *
c        print *, " al =",al, "fact/al =",fact/al
c        print *
c        print *, " bl =",be, "fact/bl =",fact/be
c        print *

        eps = 1.d0 - fact / al - fact / be
c        print *, " eps = ",eps
c        print *
        reeps = dble( eps )
        imeps = dimag( eps )
c        print *, " e2 =",imeps
c        print *
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
      write(6,*) " Completed calculation of the optical constants"
      write(6,*) " "
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


