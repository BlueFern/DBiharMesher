      subroutine dmatge(m2,n2,ki,kj,del,alpha,beta,trig,cmat,w)
c
      integer m2,n2,ki,kj
      double precision    del,alpha,beta
      double precision    trig(*)
      double precision    cmat(*)
      double precision    w(*)
c
c     this routine computes the elements of the capacitance
c     matrix defined by ki and kj (ki=1 or 2, kj=1 or 2)
c     the matrix is stored in compact form in the array cmat.
c     this storage scheme is consistent with linpack.
c     operation count is m2*n2*(n2-1)/2 (mult+add)+4*m2*n2 (mult)
c                       +2*m2*n2 (div)+7*m2*n2 (add)  (plus order m2)
c
c     cmat must have at least n2*(n2+1)/2 elements. (this call.)
c     trig is assumed initialized by two calls to routine trigin.
c     workspace w must have at least n2 elements.
c
c     local.
c
c     blas:             ddot, daxpy
c
      integer i,j,k,km,ip,m,n,i1,i2
      double precision    scal,x1,x2,x3,x4
      double precision    ddot
c
      m    = 2*(m2+ki-2)+1
      n    = 2*(n2+kj-2)+1
      i1   = (m+1)*(ki-1)
      i2   = 2*m+(n+1)*(kj-1)
      x2   = 4.0d0/(n2+kj-1.0d0)
      scal = 4.0d0*del*del/(m2+ki-1.0d0)
c
c    the following loops (5 and 15) could be simplified
c    under fortran 77 assumptions. this is correct also
c    on fortran 66.
c
      ip   = 0
      do 15 k=1,n2
         if(k.eq.1) go to 10
         km = k-1
         do 5 j=1,km
            ip      = ip+1
            cmat(ip)= 0.0d0
    5       continue
   10    ip       = ip+1
         cmat(ip) = 1.0d0
   15    continue
      do 50 i=1,m2
         x1 = scal*trig(i1+i)*trig(i1+i)
         do 20 j=1,n2
            w(j) = trig(i2+j)/((trig(i1+m2+i)+trig(i2+n2+j))*
     +                         (trig(i1+m2+i)+trig(i2+n2+j)-alpha)+
     +                         beta)
   20       continue
         x3 = ddot(n2,trig(i2+1),1,w,1)
         ip = 0
         x3 = x1*x2/(1.0d0+x2*x3)
         do 40 k=1,n2
            x4 = -x3*w(k)
            call daxpy(k-1,x4,w,1,cmat(ip+1),1)
            ip = ip+k
            cmat(ip)=cmat(ip)+w(k)*(x1/trig(i2+k)+x4)
   40       continue
   50    continue
      return
      end
