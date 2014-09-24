      subroutine spentf(n2,kj,xl,alpha,beta,trig,y,x,ws)
c
      integer n2,kj
      real    xl,alpha,beta
      real    x(*),y(*)
      real    trig(*),ws(*)
c
c     this routine helps solving special pentadiagonal systems
c     that arise when solving certain biharmonic problems.
c
c     y has length n2  (input )
c     x has length n2  (output)
c     trig is trigonometric information.
c     ws is workspace of lenght at least n2.
c
c     operation count is n2*(2(div)+4(mult)+6(add)).
c     n2(div) can be turned into n2(mult) by storing
c     1/trig(j) as well.
c
c     one mult and one add can be saved by precomputing the
c     quantity c1 (requiring slightly more storage.)
c     this will be implemented in later versions.
c
c     local.
c
c     blas:          sdot
c
      integer j
      real    x1,c1,c2
      real    sdot
c
c     form the inverse of the diagonal matrix.
c
      do 20 j=1,n2
         ws(j) = trig(j)/((xl+trig(n2+j))*(xl+trig(n2+j)-alpha)+beta)
   20   continue
      c1 = sdot(n2,trig,1,ws,1)
      c2 = sdot(n2,y,1,ws,1)
      x1 = 4.0e0/(n2+kj-1.0e0)
      c2 = x1*c2/(1.0e0+x1*c1)
      do 30 j=1,n2
         x(j) = (y(j)/trig(j)-c2)*ws(j)
   30    continue
      return
      end
