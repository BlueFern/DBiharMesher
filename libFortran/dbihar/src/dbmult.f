      subroutine dbmult(m,n,del,idf,f,idy,y,w)
c
      integer m,n,idf,idy
      double precision  del
      double precision  f(idf,*),y(idy,*)
      double precision  w(*)
c
c     multiplies the discrete 13-point biharmonic
c     operator (times hy**4) with vector f.
c     result in y , y can overwrite f if wanted.
c     note that both f and y are treated as m by n.
c     note the leading dimension of f and y ,idf and idy.
c     m and n must be at least 3.
c     del is (hy/hx)**2  (= (float(m+1)/float(n+1))**2 for square)
c     workspace w must have at least 4*m+2*n elements.
c
      integer i,j
      integer i1,i2,i3,i4,i5
      double precision  a,b,two
c
      i1  = 1
      i2  = i1+m
      i3  = i2+m
      i4  = i3+2*m
      i5  = i4+n
      two = 2.0d0
      a   = two*(del+1.0d0)
      b   = two*del*del
c
      call dcopy(m,f(1,1),1,w(i1),1)
      call dcopy(m,f(1,n),1,w(i2),1)
      call dcopy(n,f(1,1),idf,w(i4),1)
      call dcopy(n,f(m,1),idf,w(i5),1)
      call dlmult(m,n,idf,f,idy,y,a,del,w(i3))
      call dlmult(m,n,idy,y,idy,y,a,del,w(i3))
      call daxpy(m,two,w(i1),1,y(1,1),1)
      call daxpy(m,two,w(i2),1,y(1,n),1)
      call daxpy(n,b,w(i4),1,y(1,1),idy)
      call daxpy(n,b,w(i5),1,y(m,1),idy)
      return
      end
