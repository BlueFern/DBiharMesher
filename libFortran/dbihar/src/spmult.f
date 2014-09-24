      subroutine spmult(m,a,b,x,y)
c
      integer m
      real  a,b
      real  x(*),y(*)
c
c     y= tridiag(-b,a,-b) * x
c     for example a one-dimensional discrete laplace operator multiply.
c     x and y must have dimension at least m.
c     y can overwrite x if wanted.
c     m must be at least 3.
c
      integer i
      real  x1,x2
c
      x1  = x(1)
      y(1)= x1*a-x(2)*b
      do 10 i= 3,m
         x2    = x(i-1)
         y(i-1)= x2*a-b*(x1+x(i))
         x1    = x2
   10    continue
      y(m)= x(m)*a-x1*b
      return
      end
