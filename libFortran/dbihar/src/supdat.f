      subroutine supdat(nn,kj,maxu,iflag,km,tol,y,s,d,al,h,ws)
c
      integer nn,kj,maxu,iflag,km
      real    tol
      real    y(*),s(*)
      real    d(*),al(*),h(nn,*),ws(*)
c
c     updates the current preconditioning matrix
c     using a symmetric rank one qn update. the
c     vectors are kept in h (never multiplied
c     together.
c     maxu is the maximum number of updates.
c     km is current number of updates.
c     y and s are nn vectors that defines a new update.
c     d is a diagonal matrix of dimension nn used in shzero.
c     h is array of size h(nn,maxu)
c     al is array of size al(maxu) used to keep scaling of the
c        symmetric rank one updates.
c     ws is workspace of lenght at least 3*nn+18 if kj=0.
c     if kj=1 or kj=2 then ws is a dummy argument.
c
c     local.
c
c     biharmonic:         spreco
c     blas:               sdot, saxpy
c     fortran:            abs
c
      integer k
      real    sdot
c
      if(km.eq.maxu) km = km+1
      if(km.eq.maxu+1) return
      k = km+1
      call spreco(nn,kj,5,maxu,km,h(1,k),y,d,al,h,ws)
      call saxpy(nn,-1.0e0,s,1,h(1,k),1)
      al(k) = -sdot(nn,h(1,k),1,y,1)
      if(abs(al(k)).lt.tol*sdot(nn,h(1,k),1,h(1,k),1)) return
      al(k) = 1.0e0/al(k)
      km = km+1
      return
      end
