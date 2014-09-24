      subroutine sftrny(m,n,f,idf,ws)
c
      integer m,n,idf
      real    f(idf,*)
      real    ws(*)
c
c     performs m sine transforms. the transform is unscaled.
c     note that this transform must be scaled up by a factor
c     two in order to correspond to the transform described
c     in chapter 4.2 of "numerical solution of the biharmonic
c     equation."
c     n should be odd. (when ssint is called.)
c     workspace ws must have at least int(3.5*n+16) elements.
c
c     local.
c
c     fourier:            ssinti, dsint   (swarztrauber version 3.)
c     blas:               scopy
c
      integer i
c
      call ssinti(n,ws(n+2))
      do 500 i   =1,m
         call scopy(n,f(i,1),idf,ws,1)
         call ssint(n,ws,ws(n+2))
         call scopy(n,ws,1,f(i,1),idf)
  500    continue
      return
      end
