      subroutine inner(ar,ai,xr,xi,yr,yi,n)
      implicit none
      integer i,n
      double precision ar,ai,xr(n),xi(n),yr(n),yi(n)
      ar=0.d0
      ai=0.d0
      do i=1,n
        ar=ar+xr(i)*yr(i)+xi(i)*yi(i)
        ai=ai+xr(i)*yi(i)-xi(i)*yr(i)
      end do
      return
      end
