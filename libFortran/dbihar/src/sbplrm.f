      subroutine sbplrm(a,b,m,n,x0,x,idx,y0,y,idy,w)
c
      integer             m,n,idx,idy
      real    a,b,x0,y0
      real    x(idx,*),y(idy,*),w(*)
c
c******
c  july 1984.
c******
c
c     this routine multiplies a discrete biharmonic operator
c     in polar coordinates by a vector x, result in y. the
c     vector y can overwrite x.
c     this routine is the inverse of the routine sbiplr
c     if alfa = beta = 0.
c
c     if y = (sbplrm + alfa*slplrm + beta)*x
c     then a call to sbiplr with y as the right hand side
c     function, will produce the vector x. )
c
c     note that slplrm is called by this routine with
c     the special argument a = -1, in the case a=0.
c
c     if a is zero, then x0 and y0 is used to hold the
c     component of the vector corresponding to the center
c     point. y0 can overwrite x0.
c
c     a,b,m and n are the same as in sbiplr.
c     note that this routine must be called with f(2,1) when
c     f is the array used in sbiplr.
c
c     the workspace w must have at least 6*m+2*n elements.
c
      integer j
      real    a1,bm,at,xt,sum1,sum2
      real    dr,tpi,drh,zero(1)
c
      zero(1) = 0.0e0
      at  = a
      xt  = x0
      tpi = 8.0e0*atan(1.0e0)
      dr  = (b-a)/(m+1)
      drh = .50e0*dr
      a1  = 2.0e0*(a+drh)/((a+dr)*dr**4)
      bm  = 2.0e0*(b-drh)/((b-dr)*dr**4)
c
      call scopy(2*n,zero,0,w,1)
      call saxpy(n,a1,x(1,1),idx,w,1)
c
      if(a .eq. 0.0e0) then
         at = -1.0e0
         sum1 = 0.0e0
         sum2 = 0.0e0
         do 10 j=1,n
            sum1 = sum1 + x(1,j)
            sum2 = sum2 + x(2,j)
   10       continue
         y0 = 16.0e0 * a1 * (x0 + (sum2 - 4.0e0 * sum1)/(3.0e0 * N))
      end if
c
      call saxpy(n,bm,x(m,1),idx,w(n+1),1)
      call slplrm(at,b,m,n,x0,x,idx,y0,y,idy,w(2*n+1))
      call slplrm(at,b,m,n,y0,y,idy,y0,y,idy,w(2*n+1))
      call saxpy(n,1.0e0,w(n+1),1,y(m,1),idy)
c
      if(a .eq. 0.0e0) then
         sum1 = a1 * (2.0e0*sum1/n - 3.0e0*xt)
         sum2 = 3.0e0*a1*xt/8
         do 20 j=1,n
            y(1,j) = y(1,j) + sum1
            y(2,j) = y(2,j) + sum2
   20       continue
      else
         call saxpy(n,1.0e0,w,1,y(1,1),idy)
      end if
c
      end
