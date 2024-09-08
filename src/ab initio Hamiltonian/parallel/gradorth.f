      subroutine gradient(n,xr,xi,axr,axi,xx,xax,s,gr,gi)
      implicit none
c
      integer i,n
c
      double precision xr(n),xi(n),axr(n),axi(n)
      double precision gr(n),gi(n),xx,xax,s
      double precision tmpi
c
c  get (x,x), which we invert until the end
c
      gr = xr
      gi = xi
      call inner(xx,tmpi,gr,gi,xr,xi,n)
c
      xx = 1 / xx
c
c  get (x,Ax)
c
      call inner(xax,tmpi,gr,gi,axr,axi,n)
c
c  get Rayleight factor, s=(x,Ax) divided by (x,x)
c
      s=xax*xx
c
c  from whence we get the gradient
c
      do i=1,n
        gr(i)=-(axr(i)-s*xr(i))*xx
        gi(i)=-(axi(i)-s*xi(i))*xx
      end do
c
      xx = 1 / xx
c
      return
      end
c
c------------------------------------------------------------
c
      subroutine orthog(n,nvec,ivec,xr,xi,xvr,xvi)
      implicit none
c
      integer i,n,iv,nvec,ivec
c
      double precision dvr,dvi
      double precision xr(n),xi(n),xvr(n,nvec),xvi(n,nvec)
c
      do iv=1,nvec
        if (iv.ne.ivec) then
          call inner(dvr,dvi,xvr(1,iv),xvi(1,iv),xr,xi,n)
          do i=1,n
            xr(i)=xr(i)-(dvr*xvr(i,iv)-dvi*xvi(i,iv))
            xi(i)=xi(i)-(dvr*xvi(i,iv)+dvi*xvr(i,iv))
          end do
        end if
      end do
c
      return
      end
