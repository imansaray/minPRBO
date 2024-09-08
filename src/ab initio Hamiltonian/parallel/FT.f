
      program FT

c-------------------------------------------------------------
c
c     purpose: Uses a FT to switch from density as a function of 
c        real-space to reciprocal space.
c        Reads the file rhoofr. Writes the file rhoG.
c     
c     first coded by h. lawler spr 08
c     
c-------------------------------------------------------------

      INCLUDE 'license.h'
      integer j,ii,jj,l,m,n
      integer kk
      integer Ngvecsi
      character*2 ntxt
      double precision,allocatable::Gi(:,:),phi(:),rhor(:)
      double precision rho1,rho2,vol,bv1(3)
      double precision bv2(3),bv3(3),modG,dumf
      integer b(3),n1,n2,n3,irang,jrang,krang
 
      b(:) = 0
      rho1 = 0
      rho2 = 0

      open(unit=17,file='bvecs',form='formatted',status='unknown')  
       read(17,*) bv1(1),bv1(2),bv1(3)
       read(17,*) bv2(1),bv2(2),bv2(3)
       read(17,*) bv3(1),bv3(2),bv3(3)
      close(17)

      open(unit=17,file='rhoofr',form='formatted',status='unknown')   
      open(unit=27,file='Glist',form='formatted',status='unknown')

      read(17,*) ntxt

      open(unit=18,file='rho.head',form='formatted',status='unknown')
        read(18,*) n1,n2,n3,dumf
        read(18,*) vol
      close(18)

      Ngvecsi = n1*n2*n3
      allocate (Gi(3,Ngvecsi),rhor(Ngvecsi),phi(Ngvecsi))

      j=0
      do l=1,n1
       do m=1,n2
        do n=1,n3
         j = j+1
         read(17,*) Gi(1,j),Gi(2,j),Gi(3,j),rhor(j)
        enddo
       enddo
      enddo
      close(17)      

      irang = nint(n1/2.d0)
      jrang = nint(n2/2.d0)
      krang = nint(n3/2.d0)
      
      open(unit=17,file='rhoG',form='formatted',status='unknown')
      open(unit=19,file='rhoofg',form='formatted',status='unknown')
       write(19,*)(irang*2+1)*(jrang*2+1)*(krang*2+1)
      close(19)

      do ii=-irang,irang
       do jj=-jrang,jrang
        do kk=-krang,krang
         b(1) = ii
         b(2) = jj
         b(3) = kk

         j=0
         rho1 = 0
         rho2 = 0
         do l=1,n1
          do m=1,n2
           do n=1,n3
            j=j+1


            phi(j) = 8*datan(1.d0)*(b(1)*((Gi(1,j)-1)/dble(n1))
     &                   +          b(2)*((Gi(2,j)-1)/dble(n2))
     &                   +          b(3)*((Gi(3,j)-1)/dble(n3)))

            rho1 = rho1 + rhor(j)*dcos(phi(j))
            rho2 = rho2 + rhor(j)*dsin(phi(j))
           enddo
          enddo
         enddo
         modG = (bv1(1)*b(1) + bv2(1)*b(2) + bv3(1)*b(3))**2 +
     &          (bv1(2)*b(1) + bv2(2)*b(2) + bv3(2)*b(3))**2 +
     &          (bv1(3)*b(1) + bv2(3)*b(2) + bv3(3)*b(3))**2 
         write(17,'(3(i5,1x),2(f18.9,1x),f9.3)')b(1),b(2),b(3),
     &           rho1*vol/Ngvecsi,-1*rho2*vol/Ngvecsi,modG/1.d0

        enddo
       enddo
      enddo

      write(6,*)rho1*vol/Ngvecsi,rho2*vol/Ngvecsi
      do j=1,Ngvecsi
         write(27,*)Gi(1,j),Gi(2,j),Gi(3,j)
      enddo
      close(27)

      open(unit=37,file='Ng',form='formatted',status='unknown')
       write(37,*) Ngvecsi
      close(37)
      
      end
