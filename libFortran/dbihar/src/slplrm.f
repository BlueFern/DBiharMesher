      subroutine slplrm(a,b,m,n,x0,x,idx,y0,y,idy,w)
c
      integer m,n,idx,idy
      real a,b,x0,y0
      real x(idx,*),y(idy,*),w(*)
c
c     this routine multiplies a discrete approximation to the
c     laplace operator in polar coordinates by a vector x,
c     result in y. the vector y can overwrite x.
c     a,b,m,n should be the same as in the routine sbiplr.
c     the workspace w must have at least 6m elements.
c
c     if a is zero, then the center point is represented
c     by the variables x0 and y0, again y0 can overwrite x0.
c     x0 and y0 are not referenced if a is positive.
c
c     if a is -1, then a should be considered zero, but
c     the special cases for zero a in the code below
c     should not be executed. (this is the case when
c     this routine is called by sbplrm)
c
      integer i,j,nn,case
      real    dr,drh,dt,dr2,drdt,tpi
      real    t1,t2,t3
c
      if( a .eq. -1.0e0) then
         case = 1
         a    = 0.0e0
      else
         case = 0
      end if
c
      nn  = n-1
      tpi = 8.0e0*atan(1.0e0)
      dr  = (b-a)/(m+1)
      drh = .50e0*dr
      dt  = tpi/n
      dr2 = dr**2
      drdt= dr*dt
c
      t1  = a*dt
      t2  = 2.0e0*a*dr
      t3  = 1.0e0/dr2
      dr2 = 2.0e0*dr2
      do 50 i=1,m
         t1 = t1 + drdt
         t2 = t2 + dr2
         w(i) = 1.0e0/t1**2
         w(m+i) = t3-1.0e0/t2
         w(2*m+i) = t3+1.0e0/t2
   50    continue
c
c      case a=0 first part.
c
      if(a .eq. 0.0e0 .and. case .eq. 0) then
        t1 = 8.0e0/(n*dr2)
        t2 = x0/dr2
        t3 = 0.0e0
        do 60 j=1,n
           t3 = t3 + x(1,j)
   60      continue
        y0 = -8.0e0/dr2 * x0 + t1*t3
      end if
c
      call scopy(m,x(1,1),1,w(3*m+1),1)
      call scopy(m,x(1,n),1,w(5*m+1),1)
      do 100 j=1,nn
         call scopy(m,x(1,j),1,w(4*m+1),1)
         call spplrm(m,dr,w(m+1),w(2*m+1),x(1,j),y(1,j))
         do 80 i=1,m
            y(i,j) = y(i,j)+w(i)*(w(5*m+i)+x(i,j+1)-2.0e0*w(4*m+i))
   80       continue
         call scopy(m,w(4*m+1),1,w(5*m+1),1)
  100    continue
      call scopy(m,x(1,n),1,w(4*m+1),1)
      call spplrm(m,dr,w(m+1),w(2*m+1),x(1,n),y(1,n))
      do 110 i=1,m
         y(i,n) = y(i,n)+w(i)*(w(3*m+i)+w(5*m+i)-2.0e0*w(4*m+i))
  110    continue
c
c     finish the case when a is zero.
c
      if(a .eq. 0.0e0 .and. case .eq. 0) then
         do 120 j=1,n
            y(1,j) = y(1,j) + t2
  120       continue
      end if
c
      if(case .eq. 1) a = -1.0e0
c
      return
      end
