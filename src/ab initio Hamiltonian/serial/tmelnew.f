      program tmelnew 
c---------------------------------------------
c
c     purpose: writes the file tmels, 
c
c     first coded by h. lawler spr 08
c     Edited by I. Mansaray, to use ABINIT 9.8.3
c
c---------------------------------------------

      include 'license.h'
      integer npw,i,b1,b2,ik,nk,tmpo,tmpu,nspinor,nband
      integer minoband,maxoband,minuband,maxuband,totalbands,diffu
      integer, allocatable :: gvec(:,:)
      double precision orthcr, orthci,qval,q1,q2,q3
      double precision, allocatable :: zr(:,:),zi(:,:),zlin(:)
      complex(kind = kind( 1.0d0 ) ), allocatable :: tmels(:,:)
      double precision bv1(3),bv2(3),bv3(3)
      complex :: tmels_tmp
c   Number of kpoints, nkpt from ABINIT      
      open(unit=75,file='nkpt.ipt',form='formatted',
     &   status='unknown')
      read(75,*) nk
      write(6,*) " nkpt =",nk
      close(75)

c     Getting the band ranges for valence and conduction      
      open(unit=75,file='brange.ipt',form='formatted',status='unknown')
      read(75,*) minoband,maxoband
      read(75,*) minuband,maxuband
      close(75)
      write(6,*) " number of bands = ",maxuband
      totalbands=maxoband-minoband+maxuband-minuband+2
      tmpo = maxoband-minoband+1
      tmpu = maxuband-minuband+1
      write(6,*) "tmpo =",tmpo,"tmpu = ",tmpu
      write(6,*) " "
c      allocate(tmels(nk*tmpo,nk*tmpu))

c    Getting q in units of b vectors to calculate transition matrix elements
c    q = beta_1*b_1 + beta_2*b_2 + beta_3*b_3      
      open(unit=75,file='qinunitsofbvectors.ipt',form='formatted',
     &   status='unknown')
      read(75,*) q1,q2,q3
      close(75)

c    Reciprocal lattice vectors, b      
      open(unit=75,file='bvecs',form='formatted',
     &   status='unknown')
      read(75,*) bv1(1),bv1(2),bv1(3)
      read(75,*) bv2(1),bv2(2),bv2(3)
      read(75,*) bv3(1),bv3(2),bv3(3)
      close(75)


c    Magnitude of q vector
      qval =dsqrt( (q1*bv1(1)+q2*bv2(1)+q3*bv3(1))**2
     &            +(q1*bv1(2)+q2*bv2(2)+q3*bv3(2))**2
     &            +(q1*bv1(3)+q2*bv2(3)+q3*bv3(3))**2)
      write(6,*) "qval = ",qval
      write(6,*) "inv qval = ",1.0/qval

c    Opening up the wavefunction coeff file, and tmels file to store mat
c    elements      
       open(unit=75,file='u1.dat',form='unformatted',status='unknown')
       rewind 75
       open(unit=18,file='tmels',form='formatted',status='unknown')
       rewind 18
       open(unit=76,file='tmels2',form='formatted',status='unknown')
       rewind 76

       do ik = 1, nk
         read(75) npw,nspinor,nband
         write(6,*) " ik = ",ik
         write(6,*) " "
         write(6,*) "npw =",npw,"nspinor=",nspinor,"nband=",nband
         allocate (gvec(npw,3))
         allocate (zr(npw,nband),zi(npw,nband))
         allocate(zlin(2*npw*nspinor))
!         allocate(tmels(tmpo,tmpu))
!         tmels = 0.0
         read(75) gvec    
         write(6,*) "size gvec array =",size(gvec),"shape=",shape(gvec) 
         write(6,*) "shape zr, zi,zlin",shape(zr),shape(zi),shape(zlin)
         write(6,*) " "
         write(6,*) " "
         do iband = 1, nband
           read(75) zlin
           zr(:,iband) = zlin(1:npw)
           zi(:,iband) = zlin(npw+1:2*npw*nspinor)
         enddo
         diffu = maxuband - minuband
         do b2=minuband,maxuband
           do b1=minoband,maxoband
             orthcr = 0.d0
             orthci = 0.d0
             do i=1,npw
               orthcr = orthcr + zr(i,b1)*zr(i,b2) + zi(i,b1)*zi(i,b2)
               orthci = orthci + zr(i,b1)*zi(i,b2) - zi(i,b1)*zr(i,b2)
             enddo

! Original formatting
            write(18,'(8(1x,1e22.15))')orthcr/qval,orthci/qval,0.0,0.0, 
     &                                  0.0,0.0,0.0,0.0
             tmels_tmp = (1.0/qval)*cmplx(orthcr,orthci)
!             write(6,*) "occ,b1=",b1,"unocc, b2=",b2
!             write(6,*) "tmels(b1,b2)=",tmels_tmp
!             write(6,*) " "
c             write(18,*) tmels_tmp
            
!            write(76,*)orthcr/qval,orthci/qval
            write(76,*)orthcr,orthci
c  Do this to avoid going out of bounds, only needed when bwflag is not
c  set.
!             tmels(b1,b2-diffu) = (1.0/qval)*cmplx(orthcr,orthci)
           enddo
         enddo
!         write(6,*) tmels(:,:)
!         write(18,*) tmels
         deallocate (gvec,zr,zi,zlin)
!         deallocate (gvec,zr,zi,zlin,tmels)
        enddo  
      close(18)
      close(75)
      close(76)

      end
