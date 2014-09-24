      subroutine dpplrm(m,dr,amat,bmat,x,y)
c
      integer m
      double precision    dr
      double precision    amat(*),bmat(*)
      double precision    x(*),y(*)
c
c     this routine multiplies a vector x by a primitive (unsymmetric)
c     matrix used when forming the discrete laplacian in polar
c     coordinates. the result is stored in y, y can overwrite x.
c     amat and bmat must be passed computed as in dlplrm. (w(m+1) and
c     w(2*m+1) in that routine.) (m and dr are also defined there.)
c
      integer i,mm
      double precision    t1,t2,t3
c
      mm = m-1
      t2 = x(1)
      t3 = -2.0d0/dr**2
      y(1) = t3*t2+bmat(1)*x(2)
      do 10 i=2,mm
         t1 = x(i)
         y(i) = t3*x(i)+amat(i)*t2+bmat(i)*x(i+1)
         t2=t1
   10    continue
      y(m)=t3*x(m)+amat(m)*t2
      return
      end
