      subroutine dbplrm(a,b,m,n,x0,x,idx,y0,y,idy,w)
c
      integer             m,n,idx,idy
      double precision    a,b,x0,y0
      double precision    x(idx,*),y(idy,*),w(*)
c
c******
c  july 1984.
c******
c
c     this routine multiplies a discrete biharmonic operator
c     in polar coordinates by a vector x, result in y. the
c     vector y can overwrite x.
c     this routine is the inverse of the routine dbiplr
c     if alfa = beta = 0.
c
c     if y = (dbplrm + alfa*dlplrm + beta)*x
c     then a call to dbiplr with y as the right hand side
c     function, will produce the vector x. )
c
c     note that dlplrm is called by this routine with
c     the special argument a = -1, in the case a=0.
c
c     if a is zero, then x0 and y0 is used to hold the
c     component of the vector corresponding to the center
c     point. y0 can overwrite x0.
c
c     a,b,m and n are the same as in dbiplr.
c     note that this routine must be called with f(2,1) when
c     f is the array used in dbiplr.
c
c     the workspace w must have at least 6*m+2*n elements.
c
      integer j
      double precision    a1,bm,at,xt,sum1,sum2
      double precision    dr,tpi,drh,zero(1)
c
      zero(1) = 0.0d0
      at  = a
      xt  = x0
      tpi = 8.0d0*atan(1.0d0)
      dr  = (b-a)/(m+1)
      drh = .50d0*dr
      a1  = 2.0d0*(a+drh)/((a+dr)*dr**4)
      bm  = 2.0d0*(b-drh)/((b-dr)*dr**4)
c
      call dcopy(2*n,zero,0,w,1)
      call daxpy(n,a1,x(1,1),idx,w,1)
c
      if(a .eq. 0.0d0) then
         at = -1.0d0
         sum1 = 0.0d0
         sum2 = 0.0d0
         do 10 j=1,n
            sum1 = sum1 + x(1,j)
            sum2 = sum2 + x(2,j)
   10       continue
         y0 = 16.0d0 * a1 * (x0 + (sum2 - 4.0d0 * sum1)/(3.0d0 * N))
      end if
c
      call daxpy(n,bm,x(m,1),idx,w(n+1),1)
      call dlplrm(at,b,m,n,x0,x,idx,y0,y,idy,w(2*n+1))
      call dlplrm(at,b,m,n,y0,y,idy,y0,y,idy,w(2*n+1))
      call daxpy(n,1.0d0,w(n+1),1,y(m,1),idy)
c
      if(a .eq. 0.0d0) then
         sum1 = a1 * (2.0d0*sum1/n - 3.0d0*xt)
         sum2 = 3.0d0*a1*xt/8
         do 20 j=1,n
            y(1,j) = y(1,j) + sum1
            y(2,j) = y(2,j) + sum2
   20       continue
      else
         call daxpy(n,1.0d0,w,1,y(1,1),idy)
      end if
c
      end
