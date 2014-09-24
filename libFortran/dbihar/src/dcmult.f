      subroutine dcmult(mm,nn,ki,kj,del,alpha,beta,x,y,trig,ws)
c
      integer mm,nn,ki,kj
      double precision    del,alpha,beta
      double precision    x(*),y(*)
      double precision    trig(*),ws(*)
c
c     dcmult defines the capacitance-matrix.
c     y is c*x where c represents the matrix.
c
c     local.
c
c     biharmonic:         dpentf
c     blas:               dcopy,daxpy
c
      integer i,m,n,m2,n2,i1,i2
      double precision    x1,scal
c
c     this is the case where two ffts are used.
c
      m2   = mm
      n2   = nn
      m    = 2*(m2+ki-2)+1
      n    = 2*(n2+kj-2)+1
      i1   = (m+1)*(ki-1)
      i2   = 2*m+1+(n+1)*(kj-1)
      scal = 4.0d0*del*del/(m2+ki-1.0d0)
      call dcopy(n2,x,1,y,1)
      do 30 i=1,m2
         x1 = scal*trig(i1+i)*trig(i1+i)
         call dpentf(n2,kj,trig(i1+m2+i),alpha,beta,
     +              trig(i2),x,ws,ws(n2+1))
         call daxpy(n2,x1,ws,1,y,1)
   30 continue
      return
      end
