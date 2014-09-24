      subroutine slmult(m,n,id1,f,id2,y,a,b,w)
c
      integer m,n,id1,id2
      real  a,b
      real  f(id1,*),y(id2,*)
      real  w(*)
c
c     y= blocktridiag(-i,p,-i) * f where p is tridiag(-b,a,-b).
c     ie minus 1 times some discrete laplace-operator times f.
c     x1 and x2 must have at least m elements.
c     y can overwrite f if wanted.
c     note the first dimension of f and y , id1 and id2.
c     f and y are treated as m by n in the routine.
c     m,n  must be at least 3.
c     workspace w must have at least 2*m elements.
c
      integer i,j,k
c
      k  = m
      call scopy(m,f(1,1),1,w,1)
      call spmult(m,a,b,f(1,1),y(1,1))
      call saxpy(m,-1.0e0,f(1,2),1,y,1)
c
      do 50 j=3,n
         call scopy(m,f(1,j-1),1,w(k+1),1)
         call spmult(m,a,b,f(1,j-1),y(1,j-1))
         call saxpy(m,-1.0e0,w(m-k+1),1,y(1,j-1),1)
         call saxpy(m,-1.0e0,f(1,j),1,y(1,j-1),1)
         k=m-k
   50    continue
c
      call spmult(m,a,b,f(1,n),y(1,n))
      call saxpy(m,-1.0e0,w(m-k+1),1,y(1,n),1)
      return
      end
