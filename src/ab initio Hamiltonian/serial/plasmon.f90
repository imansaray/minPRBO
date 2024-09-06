!     Written by J. Vinson Dec. 08
!     The purpose of this is outlined in plasmon.tex hopefully nearby

      program plasmon

      implicit none

      integer :: gtot,occ_min,occ_max,unocc_min,unocc_max,bandtot
      integer, allocatable :: gvecs(:,:)
      double precision, allocatable :: kr(:,:),ki(:,:),G(:,:)
      double precision :: fermie,pi,energy,mult,invpi,broaden,          &
     &   broadensq,omega(3,3),bigA,gscale(3,3),kpt_un(3),kpt_sh(3),     &
     &   mult_total,kpt(3),vol

      integer :: filecounter,numfiles,dumint,bandcounter,i,j,gcounter
      character*11 :: filename
      pi=3.14159
      invpi=1/pi
!      broaden=.000007
!      broadensq=broaden**2
      omega(:,:) = 0.d0
      mult_total = 0.d0

! Open and read the band ranges since this data not in wf files
      open(unit=21,file='brange.ipt',form='formatted',status='old')
      read(21,*)occ_min,occ_max,unocc_min,unocc_max
      close(21)
!      write(6,'(4I5)')occ_min,occ_max,unocc_min,unocc_max
      bandtot=unocc_max-unocc_min+occ_max-occ_min+2

! Get the fermi energy
      open(22,file='efermiinrydberg.ipt',form='formatted',status='old')
      read(22,*) fermie
      close(22)

! Get number of wf files
      open(unit=22,file='masterwfile',form='formatted',status='old')
      read(22,*)numfiles
      close(22)

! Get the rescaling of the gvectors
      open(unit=22,file='bvecs',form='formatted',status='old')
      read(22,*)gscale(1,1),gscale(1,2),gscale(1,3)
      read(22,*)gscale(2,1),gscale(2,2),gscale(2,3)
      read(22,*)gscale(3,1),gscale(3,2),gscale(3,3)
      close(22)

! get the broaden param
      open(unit=22,file='pbroaden',form='formatted',status='old')
      read(22,*) broaden
      close(22)
      broadensq = broaden**2

! get the volume
      open(unit=22,file='vol',form='formatted',status='old')
      read(22,*) vol
      close(22)

! listwfile contains the names of the wf files
      open(unit=23,file='listwfile',form='formatted',status='old')

! enkfile has all the bandenergy info
      open(unit=24,file='enkfile',form='formatted',status='old')

! kfile contains a list of the kpoints
      open(unit=25,file='klist',form='formatted',status='old')

      do filecounter=1,numfiles
! Get the filename, openfile, allocate arrays, read in everything, close
       read(23,'(I6,2X,A11)')dumint,filename
       open(30,file=filename,form='unformatted',status='old')
        
       read(30)gtot
       allocate(gvecs(gtot,3))
       allocate(kr(gtot,bandtot),ki(gtot,bandtot))
       read(30) gvecs
       read(30) kr
       read(30) ki
       close(30)

! convert gvecs to G
       allocate(G(gtot,3))
       G(:,:) = 0.d0
       do gcounter=1,gtot
        do i=1,3
         do j=1,3
          G(gcounter,i) = G(gcounter,i) + gscale(i,j) * gvecs(gtot,j)
         enddo
        enddo
       enddo
       deallocate(gvecs)

! get k-point
      kpt_un(:) = 0.d0
       read(25,*)kpt(1),kpt(2),kpt(3)
       do i=1,3
        do j=1,3
         kpt_un(i) = kpt_un(i) + gscale(i,j) * kpt(j)
        enddo
       enddo
       read(25,*)kpt_sh(1),kpt_sh(2),kpt_sh(3)

! First run over all the occupied bands
       do bandcounter=1,occ_max

! get the energy and compute the multiplicative weight based on a Lorentz
! wrt the fermi energy and the broaden param. Extra pi from def. of omega
        read(24,*) energy
        mult = invpi*invpi*broaden / ((energy-fermie)**2 + broadensq) 
        mult_total = mult_total + mult

! This is equation (4) in tex file.
        do gcounter=1,gtot
         bigA=(kr(gcounter,bandcounter)**2+                             &
     &         ki(gcounter,bandcounter)**2)**2
        

         do i=1,3
          do j=1,3
           omega(i,j)=omega(i,j)+mult*bigA*(G(gcounter,i)+kpt_un(i))*   &
     &                (G(gcounter,j)+kpt_un(j))
          enddo
         enddo

        enddo ! gcounter=1,gtot

       enddo ! bandcounter=1,occ_max

       deallocate(G,kr,ki)

      enddo ! filecounter=1,numfiles
 
      omega(:,:) = sqrt(omega(:,:)) * 27.21138386 / vol
      write(6,*) mult_total,((pi*mult_total)**0.5)/(4*pi)
      write(6,*)omega(1,1),omega(1,2),omega(1,3)
      write(6,*)omega(2,1),omega(2,2),omega(2,3)
      write(6,*)omega(3,1),omega(3,2),omega(3,3)

      
      close(23)
      close(24)
      close(25)

      end
