!
! Written by J. Vinson Dec. 09
! Purpose: Read in enkfile, take the functional derivative of the energy
! wrt the bvecs to form the plasmon pole tensor.
!

      program plasmon

      implicit none

      integer :: occ_min,occ_max,unocc_min,unocc_max,numk(3),bandcount, &
     & totalk,occbandtot,i,kcounti,kcountj,kcountk,found,bandtot
      integer, allocatable :: below(:),above(:)
      double precision :: fermie,pi,omega,dumf,bvecs(3,3),delta(3),ucvol&
          ,rydberg,kvol
      double precision, allocatable :: energy(:,:,:,:)

      pi=3.14159
      rydberg = 13.60569193

! Open and read the band ranges since this data not in wf files
      open(unit=21,file='brange.ipt',form='formatted',status='old')
      read(21,*)occ_min,occ_max,unocc_min,unocc_max
      close(21)
      occbandtot=occ_max-occ_min+1
      bandtot=occ_max-occ_min+unocc_max-unocc_min+2

! Get the fermi energy
      open(21,file='efermiinrydberg.ipt',form='formatted',status='old')
      read(21,*) fermie
      close(21)

! Get the dimensions of the k-point grid
      open (21,file='nkpt',form='formatted',status='old')
      read(21,*) numk(1:3)
      close(21)
      totalk = numk(1)*numk(2)*numk(3)

! Get the volume 
      open(21,file='ucvol',form='formatted',status='old')
      read(21,*) ucvol
      close(21)

! Open up enkfile and populate energy(:,:)      
      open(21,file='enkfile',form='formatted',status='old')

      allocate( energy(numk(1),numk(2),numk(3),bandtot) )

      do kcounti=1,numk(1)
       do kcountj=1,numk(2)
        do kcountk=1,numk(3)
          read(21,*) energy(kcounti,kcountj,kcountk,1:bandtot)
!         enddo
!         do bandcount=unocc_min,unocc_max
!          read(21,*)dumf
!         enddo
        enddo
       enddo
      enddo
      close(21)

! Open up klist 
!      open(21,file='klist',form='formatted',status='old')
!      allocate(klist(totalk,3))
!      do counter=1,totalk
!       read(21,*) klist(counter,1:3)
!       read(21,*) dumf,dumf,dumf
!      enddo
!      close(21)

! need to read in the scaling of the divisions in each direction
      open(21,file='bvecs',form='formatted',status='old')
      read(21,*) bvecs(1:3,1:3)
      close(21)
      bvecs(:,:) = bvecs(:,:) * 2 * pi
      do i=1,3
       delta(i)=(bvecs(i,1)**2+bvecs(i,2)**2+bvecs(i,3)**2)/(numk(i)**2)
      enddo

! kvol is the infinitesimal for the integral approx
      kvol =bvecs(1,1) * (bvecs(2,2)*bvecs(3,3) - bvecs(3,2)*bvecs(2,3))&
     &    - bvecs(2,1) * (bvecs(1,2)*bvecs(3,3) - bvecs(3,2)*bvecs(1,3))&
     &    + bvecs(3,1) * (bvecs(1,2)*bvecs(2,3) - bvecs(2,2)*bvecs(1,3))

      kvol = kvol/(numk(1)*numk(2)*numk(3))
 

! Run over bands to find which ones might have crossings
      allocate(below(occbandtot),above(occbandtot))
      above(:) = 0
      below(:) = 0

      do bandcount=occ_min,occ_max
       do kcounti=1,numk(1)
        do kcountj=1,numk(2)
         do kcountk=1,numk(3)
          if (energy(kcounti,kcountj,kcountk,bandcount) .gt.            &
     &        fermie ) then
           above(bandcount) = above(bandcount) + 1
          else
           below(bandcount) = below(bandcount) + 1
          endif
         enddo
        enddo
       enddo
      enddo


      omega = 0.d0
      found = 0
      do bandcount=occ_min,occ_max
       if ((above(bandcount) .gt. 0 ) .and. (below(bandcount) .gt. 0 )) &
     &  then
       write(6,*) bandcount
        do kcounti=1,numk(1)
         do kcountj=1,numk(2)
          do kcountk=1,numk(3) 
           if (above(bandcount) .gt. below(bandcount) ) then
! there are more points above, so search starting with those below
! (this isn't needed, but will be slightly faster
            if (energy(kcounti,kcountj,kcountk,bandcount) .lt.          &
     &          fermie ) then ! 
             call neighbors(energy(1,1,1,bandcount),fermie,numk,kcounti,&
     &                      kcountj,kcountk,omega,delta,found,1)
            endif
           else  ! There are more points below
            if (energy(kcounti,kcountj,kcountk,bandcount) .gt.          &
                 &          fermie ) then
             call neighbors(energy(1,1,1,bandcount),fermie,numk,kcounti,&
     &                      kcountj,kcountk,omega,delta,found,-1)
            endif
           endif             
          enddo
         enddo
        enddo
       endif
      enddo

      deallocate(energy,below,above)
      write(6,*) "Found:",found, "out of ",totalk, "kpts"
! need to scale omega by the factors in front of the integral
      omega = 8 * pi *omega / (3 * ucvol )
! need to scale omega by the infinitesimal since we did sum not integral
      omega = omega * kvol
! our expression was for omega squared in hartree, so moving to eV
      omega = dsqrt(omega) * 2 * rydberg
      write(6,*) "omega: ",omega

      end

      subroutine neighbors (energy,fermie,numk,kcounti,kcountj,kcountk, &
     &    omega,delta,found,compare)
! if compare is +1, then the found point is below the fermi, looking for points above
! otherwise the point is above, looking for below
      integer :: numk(3),kcounti,kcountj,kcountk,i,j,kcompare(3),       &
     &           compare,found,kcount(3)
      double precision :: omega, energy(numk(1),numk(2),numk(3)),       &
     & fermie,mag,delta(3)

      kcount(1) = kcounti
      kcount(2) = kcountj
      kcount(3) = kcountk

! looking only at nearest k-point neighbors, this should probably be modified
! for highly non-cubic k-space grids
! i goes between the three directions
! j switches between the two neighbors in that direction
! then we test to see if we have wrapped around our grid (1,N), 
!   setting 0 -> N and N+1 -> 1
      do i=1,3 
       kcompare(:) = kcount(:)
       do j=-1,1,2
        if (j+kcount(i) .le. 0 ) then
         kcompare(i) = numk(i)
        else if (j+kcount(i) .gt. numk(i) ) then
         kcompare(i) = 1
        else
         kcompare(i) = j+kcount(i)
        endif

        if((energy(kcompare(1),kcompare(2),kcompare(3))-fermie)*compare &
     &       .gt. 0 ) then
         found = found + 1
         mag = 1 / delta(i)
! delta is already squared
         mag = mag *( ( energy(kcompare(1),kcompare(2),kcompare(3)) -   &
     &          energy(kcount(1),kcount(2),kcount(3)) )/2)**2 
! need to convert energies to hartree from rydbergs hence the 2 in previous exp.

         omega = omega + mag

        endif
       enddo
      enddo
      
      end subroutine neighbors


      
