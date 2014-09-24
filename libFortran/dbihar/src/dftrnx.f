      subroutine dftrnx(m,n,f,idf,ws)
c
      integer m,n,idf
      double precision    f(idf,*)
      double precision    ws(*)
c
c     performs n sine transforms. the transform is unscaled.
c     note that this transform must be scaled up by a factor
c     two in order to correspond to the transform described
c     in chapter 4.2 of "numerical solution of the biharmonic
c     equation."
c     m must be odd.(when dsint is called.)
c     workspace ws must have at least int(2.5*m+15) elements.
c     note that routine dsint overwrites f(m+1,j).
c
c     local.
c
c     fourier:            dsinti, dsint   (swarztrauber version 3.)
c
      integer j
      double precision    x1
c
      call dsinti(m,ws)
      do 500 j=1,n
         x1       = f(m+1,j)
         call dsint(m,f(1,j),ws)
         f(m+1,j) = x1
  500    continue
      return
      end
