      subroutine sbmult(m,n,del,idf,f,idy,y,w)
c
      integer m,n,idf,idy
      real  del
      real  f(idf,*),y(idy,*)
      real  w(*)
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
      real  a,b,two
c
      i1  = 1
      i2  = i1+m
      i3  = i2+m
      i4  = i3+2*m
      i5  = i4+n
      two = 2.0e0
      a   = two*(del+1.0e0)
      b   = two*del*del
c
      call scopy(m,f(1,1),1,w(i1),1)
      call scopy(m,f(1,n),1,w(i2),1)
      call scopy(n,f(1,1),idf,w(i4),1)
      call scopy(n,f(m,1),idf,w(i5),1)
      call slmult(m,n,idf,f,idy,y,a,del,w(i3))
      call slmult(m,n,idy,y,idy,y,a,del,w(i3))
      call saxpy(m,two,w(i1),1,y(1,1),1)
      call saxpy(m,two,w(i2),1,y(1,n),1)
      call saxpy(n,b,w(i4),1,y(1,1),idy)
      call saxpy(n,b,w(i5),1,y(m,1),idy)
      return
      end
