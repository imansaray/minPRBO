      subroutine solfrn( v1, v2, av1, av2, n, get, nhave,
     &                   nloop, iloop, nrest, irest, nit, iit, ii,
     &                   vr, vi, avr, avi, w, flu, tol,
     &                   wrk, idwrk, cmd, need, point, itau,
     &                   vec, xxtrue, xx, xax )
      implicit none
      integer n, get, nhave, idwrk, need, point, itau, vec
      double precision xxtrue, xx, xax
      double precision v1( n ), v2( n ), av1( n ), av2( n )
      double precision vr( n * get ), vi( n * get )
      double precision avr( n * get ), avi( n * get )
      double precision w( get ), flu, tol
      double precision wrk( idwrk )
      character * 3 cmd
      integer xr, xi, axr, axi, gr, gi, hr, hi 
      integer dr, di, gor, goi, ahr, ahi
      integer nloop, iloop, nrest, irest, nit, iit, ii
c
        xr = 1
        xi = n + xr
       axr = n + xi
       axi = n + axr
        gr = n + axi 
        gi = n + gr
        hr = n + gi
        hi = n + hr
        dr = n + hi
        di = n + dr
       gor = n + di
       goi = n + gor
       ahr = n + goi
       ahi = n + ahr
      need = n + ahi - 1
      if ( need .gt. idwrk ) then
         if ( point .ne. 0 ) stop 'bad wrk change'
         cmd = 'all'
      else
         call xasbck( n, get, nhave, itau,
     &                nloop, iloop, nrest, irest, nit, iit, ii,
     &                point, cmd, v1, v2, av1, av2,
     &                vr, vi, avr, avi, w, flu, tol,
     &                vec, xxtrue, xx, xax,
     &                wrk( xr ), wrk( xi ), wrk( axr ), wrk( axi ),
     &                wrk( gr ), wrk( gi ), wrk( hr ), wrk( hi ),
     &                wrk( dr ), wrk( di ), wrk( gor ), wrk( goi ),
     &                wrk( ahr ), wrk( ahi ) )
      end if
      return
      end
c
c-------------------------------------------------------
c
      subroutine xasbck( n, get, nhave, itau,
     &                   nloop, iloop, nrest, irest, nit, iit, ii,
     &                   point, cmd, v1, v2, av1, av2,
     &                   vr, vi, avr, avi, w, flu, tol,
     &                   vec, xxtrue, xx, xax,
     &                   xr, xi, axr, axi,
     &                   gr, gi, hr, hi, dr, di, gor, goi, ahr, ahi )
      implicit none
c
      integer n, vec, get, nhave, itau, point
      integer nloop, nrest, nit, iloop, irest, iit, i, ii, itmp
c
      character * 3 cmd
c
      double precision xx, xxtrue, xax, gg, dd, hh
      double precision xhr, xhi, lamr, lami, tmpr, tmpi, gamr, gami
      double precision an, bnr, bni, bn2, cn, eps, norm
      double precision flu, s, tol
c
      double precision v1( n ), v2( n ), av1( n ), av2( n )
      double precision xr( n ), xi( n ), axr( n ), axi( n )
      double precision vr( n * get ), vi( n * get )
      double precision avr( n * get ), avi( n * get ), w( get )
      double precision gr( n ), gi( n ), hr( n ), hi( n )
      double precision dr( n ), di( n ), gor( n ), goi( n )
      double precision ahr( n ), ahi( n )
c
      point = abs( point )
      do while ( point .ge. 0 )
         select case ( point )
            case( 0 )
               iloop = 1
               tol = 0.001
               point = 5
            case( 5 )
               if ( iloop .le. nloop ) then
                  point = 10
               else
                  point = 510
               end if 
            case( 10 )
               vec = 1
               point = 15
            case( 15 )
               if ( vec .le. get ) then
                  point = 20
               else
                  point = 420
               end if
            case( 20 )
               ii = 1 + ( vec - 1 ) * n
               if ( vec .gt. nhave ) then
                  point = -25
                  cmd = 'sup'
               else    
                  point = 25
               end if
            case( 25 )
               flu = tol + 1
               irest = 0
               point = 30
            case( 30 )
               if ( ( flu .gt. tol ) .and. ( irest .lt. nrest ) ) then
                  point = 40
               else
                  point = 340 
               end if
            case( 40 )
               do i = 1, n
                  xr( i ) = vr( ii + i - 1 )
                  xi( i ) = vi( ii + i - 1 )
               end do
               call orthog( n, nhave, vec, xr, xi, vr, vi )
               point = -50
               cmd = 'act'
               v1 = xr
               v2 = xi
            case( 50 )
               axr = av1
               axi = av2
               call gradient(n,xr,xi,axr,axi,xx,xax,s,gr,gi)
               point= -55
               cmd = 'prc'
               v1 = gr
               v2 = gi
            case( 55 )
               hr = av1
               hi = av2
               dr = axr - s * xr
               di = axi - s * xi
               call inner( dd, tmpi, dr, di, dr, di, n )
               flu = dd / xx
               iit = 1
               point = 60
            case( 60 )
               if ( ( iit .le. nit ) .and. ( flu .ge. tol ) ) then
                  point = 70
               else
                  point = 270
               end if
            case( 70 )
               call orthog(n,nhave,vec,xr,xi,vr,vi)
               call orthog(n,nhave,vec,hr,hi,vr,vi)
               call inner(xhr,xhi,xr,xi,hr,hi,n)
               call inner(xxtrue,tmpi,xr,xi,xr,xi,n)
               xhr = xhr / xxtrue
               xhi = xhi / xxtrue 
               hr = hr - xr * xhr + xi * xhi
               hi = hi - xr * xhi - xi * xhr
               point = -80
               cmd = 'act'
               v1 = hr
               v2 = hi
            case( 80 )
               ahr = av1
               ahi = av2
               call inner(hh,tmpi,hr,hi,hr,hi,n)
               call inner(an,tmpi,hr,hi,ahr,ahi,n)
               an=an/hh
               call inner(bnr,bni,ahr,ahi,xr,xi,n)
               bnr = bnr / sqrt( hh * xxtrue )
               bni = bni / sqrt( hh * xxtrue )
               cn=xax/xx
               bn2 = bnr ** 2 + bni ** 2
               eps=((an+cn)-sqrt((an-cn)**2+4*bn2))/2
               lamr=bnr/(eps-an)*sqrt( xxtrue / hh )
               lami=bni/(eps-an)*sqrt( xxtrue / hh )
               xr = xr + lamr * hr - lami * hi
               xi = xi + lamr * hi + lami * hr
               axr = axr + lamr * ahr - lami * ahi
               axi = axi + lamr * ahi + lami * ahr
               call inner( gg, tmpi, gr, gi, gr, gi, n )
               gor = gr
               goi = gi
               call gradient(n,xr,xi,axr,axi,xx,xax,s,gr,gi)
               call inner(gamr,gami,gr,gi,gr,gi,n)
               call inner(tmpr,tmpi,gor,goi,gr,gi,n)
               gamr=(gamr-tmpr)/gg
               gami=(gami-tmpi)/gg
               hr=gr+gamr*hr-gami*hi
               hi=gi+gamr*hi+gami*hr
               dr = axr - s * xr
               di = axi - s * xi
               call inner(dd,tmpi,dr,di,dr,di,n)
               flu = dd / xx
               itmp=iit / itau
               if (iit.eq.itmp*itau) then
                  open(unit=12,form='formatted',status='unknown')
                  rewind 12
                  write (  6,'(3i8,2(1x,1e15.8))' ) n,iit,nit,s,flu
                  write ( 12,'(3i8,2(1x,1e15.8))' ) n,iit,nit,s,flu
                  close(unit=12)
               end if
               norm = 1 / sqrt( xx )
               xr = xr * norm
               xi = xi * norm
               axr = axr * norm
               axi = axi * norm
               gr = gr / norm
               gi = gi / norm
               hr = hr / norm
               hi = hi / norm
               gg = gg / norm ** 2
               iit = iit + 1
               point = 60
            case ( 270 ) 
               do i = 1, n
                  vr( ii + i - 1 ) = xr( i )
                  vi( ii + i - 1 ) = xi( i )
               end do
               irest = irest + 1
               point = 30
            case( 340 )
               call inner(norm,tmpi,vr(ii),vi(ii),vr(ii),vi(ii),n)
               norm = 1 / sqrt( norm )
               do i=0,n-1
                  vr(i+ii)=vr(i+ii)*norm
                  vi(i+ii)=vi(i+ii)*norm
               end do
               point = -350
               cmd = 'act'
               v1 = vr( ii : ii + n - 1 )
               v2 = vi( ii : ii + n - 1 )
            case( 350 )
               avr( ii : ii + n - 1 ) = av1
               avi( ii : ii + n - 1 ) = av2
               if ( vec .gt. nhave ) nhave = vec
               call subspace( n, nhave, vr, vi, avr, avi, w )
               vec = vec + 1
               point = 15
            case( 420 )
               tol = tol / 2
               iloop = iloop + 1
               point = 5
            case( 510 )
               point = -1
               cmd = 'end'
         end select
      end do
      return
      end
