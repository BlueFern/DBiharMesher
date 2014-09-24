      subroutine dpreco(nn,kj,iflag,maxu,km,p,z,d,al,h,ws)
c
      integer nn,kj,iflag,km,maxu
      double precision    p(*),z(*)
      double precision    d(*),al(*),h(nn,*)
      double precision    ws(*)
c
c     preconditioning using a symmetric rank
c     one approximation to improve an initial
c     preconditioner defined in subroutine dhzero.
c     p=hz where h is the matrix defined above.
c     km rank one correctins are performed.
c     p and z are nn vectors.
c     d is diagonal of dimension nn used in shzero.
c     h is array of size h(nn,maxu)
c     al is array of size al(maxu) holding scaling factors
c     from the symmetric rank one update.
c     ws is workspace with at least 3*nn+18 elements if kj=0.
c     if kj=1 or kj=2 then ws is a dummy argument.
c
c     local.
c
c     biharmonic:         dhzero
c     blas:               ddot, daxpy
c     fortran:            min
c
      integer i,k
      double precision    x1
      double precision    ddot
c
      call dhzero(nn,kj,z,p,d,ws)
      if(iflag.le.2) return
      k = min(km,maxu)
      if(k.eq.0) return
      do 20 i = 1,k
         x1 = ddot(nn,h(1,i),1,z,1)*al(i)
         call daxpy(nn,x1,h(1,i),1,p,1)
   20    continue
      return
      end
