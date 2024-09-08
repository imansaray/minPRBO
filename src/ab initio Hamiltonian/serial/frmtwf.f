      program frmtwf
c-------------------------------------------
c
c     purpose: writes wfout, the binary form of the wavefunctions
c
c     first coded by h. lawler spr 08
c
c-------------------------------------------

      include 'license.h'
      integer npw,i,j,ii,bandi,nband
      integer, allocatable :: gVec(:,:)
      double precision, allocatable :: zr(:,:),zi(:,:)
      open(unit=74,file='wfraw',form='formatted',status='unknown')
      read(74,*)nband
      read(74,*)npw
      allocate (gVec(npw,3))
      allocate (zr(npw,nband),zi(npw,nband))
      do i=1,npw

      read(74,*) (gvec(i,j),j=1,3)     
      enddo
      do ii=1,nband

       do i=1,npw
        read(74,*) bandi,zr(i,ii),zi(i,ii) 
       enddo
      enddo
      close(74)
      open(unit=75,file='wfout',form='unformatted',status='unknown')
       write(75) npw
       write(75) gvec     
       write(75) zr
       write(75)zi
      close(75)

      end
