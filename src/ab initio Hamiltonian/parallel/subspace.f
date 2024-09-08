      subroutine subspace(n,nvec,xvr,xvi,axr,axi,w)
      implicit none
c
      integer i,j,n,nvec,matz,ierr
      double precision prec
c
      double precision xvr(n,nvec),xvi(n,nvec),w(nvec)
      double precision axr(n,nvec),axi(n,nvec)
c
      double precision, allocatable :: assr(:,:), assi(:,:)
      double precision, allocatable :: evr(:,:), evi(:,:)
      double precision, allocatable :: fv1(:), fv2(:), fm1(:)
c
      allocate( assr(nvec,nvec), assi(nvec,nvec) )
      allocate( evr(nvec,nvec), evi(nvec,nvec) )
      allocate( fv1(nvec), fv2(nvec), fm1(2*nvec) )
c
      prec=0.000001d0
c
c  build up ss matrix
c
      do i=1,nvec
        do j=1,nvec
          call inner(assr(j,i),assi(j,i),xvr(1,j),xvi(1,j),
     &               axr(1,i),axi(1,i),n)
        end do
      end do
c
c  diagonalize
c
      matz=1
      ierr=0
      call elsch(nvec,nvec,assr,assi,w,matz,
     &           evr,evi,fv1,fv2,fm1,ierr)
      write (6,'(2x,6f10.5)') (w(i),i=1,nvec)
      open(unit=99,file='eigvals',form='formatted',
     &     status='unknown')
      rewind 99
      write (99,'(2x,6f10.5)') (w(i),i=1,nvec)
      close(unit=99)
      if (ierr.ne.0) stop 'error in the eispack'
c
c  get the new deals
c
      call ssdrot( n, nvec, xvr, xvi, evr, evi )
      call ssdrot( n, nvec, axr, axi, evr, evi )
c
      deallocate( assr, assi, evr, evi, fv1, fv2, fm1 )
c
      return
      end
c
c------------------------------------------------------------
c
      subroutine ssdrot( n, nvec, xr, xi, zr, zi )
      implicit none
c
      integer n, nvec
c
      double precision xr(n,nvec), xi(n,nvec)
      double precision zr(nvec,nvec), zi(nvec,nvec)
c
      integer i, j, k
c
      double precision tmpr, tmpi
      double precision, allocatable :: br(:),bi(:),cr(:),ci(:)
c
      allocate( br(nvec), bi(nvec), cr(nvec), ci(nvec) )
c
      do i=1,n
        do j=1,nvec
          br(j)=xr(i,j)
          bi(j)=xi(i,j)
        end do
        do k=1,nvec
          tmpr=0.d0
          tmpi=0.d0
          do j=1,nvec
            tmpr=tmpr+br(j)*zr(j,k)-bi(j)*zi(j,k)
            tmpi=tmpi+br(j)*zi(j,k)+bi(j)*zr(j,k)
          end do
          cr(k)=tmpr
          ci(k)=tmpi
        end do
        do j=1,nvec
          xr(i,j)=cr(j)
          xi(i,j)=ci(j)
        end do
      end do
      deallocate( br, bi, cr, ci )
      return
      end
