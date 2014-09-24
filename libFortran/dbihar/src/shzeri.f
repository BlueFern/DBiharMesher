      subroutine shzeri(m,n,kj,del,alpha,beta,diag,trig,ws)
c
      integer m,n,kj
      real    del,alpha,beta
      real    diag(*),trig(*),ws(*)
c
c     this routine computes (initializes) the diagonal
c     matrix diag that is needed in routine shzero.
c     diag has at least 2*n elements.
c     ws   has at least m/2+1 elements.
c
c     local.
c
      integer i,j,i1,i2,n2,m2,incr,kii,kjj,ip
      real    x1,x2,x3
c
      x1 = .1250e0/(n+1.0e0)
      x2 = 8.0e0*del*del/(m+1.0e0)
      incr=1
      if(kj.eq.0) incr=2
      ip = 0
      do 800 kjj=1,2
         n2 = n/2+2-kjj
         i2 = 2*m+(n+1)*(kjj-1)
         do 700 kii = 1,2
            m2 = m/2+2-kii
            i1 = (m+1)*(kii-1)
            if(kj.eq.0) ip = (kii-1)*n+kjj-2
            do 50 i = 1,m2
               ws(i) = trig(i1+i)**2
   50          continue
            do 200 j = 1,n2
               x3 = 0.0e0
               ip = ip+incr
               do 100 i = 1,m2
                  x3 = x3+ws(i)/((trig(i1+m2+i)+trig(i2+n2+j))*
     +                         (trig(i1+m2+i)+trig(i2+n2+j)-
     +                          alpha)+beta)
  100             continue
               x3 = x2*x3+1.0e0
               diag(ip) = x1/x3
  200          continue
  700       continue
  800    continue
      return
      end
