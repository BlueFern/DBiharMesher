      subroutine dstart(m,n,alpha,beta,f,idf,bda,bdb,bdc,bdd,dx,dy,del)
c
      integer m,n,idf
      double precision    dx,dy,del,alpha,beta
      double precision    bda(*),bdb(*),bdc(*),bdd(*)
      double precision    f(idf,*)
c
c     this routine computes the right hand side
c     of the discrete system.
c
c     input
c           m,n,alpha,beta,idf,f,bda,bdb,bdc,bdd,dx,dy,del
c     output
c           f
c
c     local.
c
c     blas:          dscal
c
      integer mp,np
      integer i,j
      double precision    twody,twodel,del2,dy4,d1,d2,d3
c
      mp    =m+1
      np    =n+1
      twody =2.0d0*dy
      twodel=2.0d0*del
      del2  =del*del
      dy4   =dy**4.0d0
      twody =2.0d0*dy
      twodel=2.0d0*del
      dy4   =dy*dy*dy*dy
      d1    =twodel+twodel+4.0d0-alpha
      d2    =del*d1
      d3    =2.0d0*dx*del2
c
c     scale right hand side.
c
      do 300 j=2,np
         call dscal(m,dy4,f(2,j),1)
  300    continue
c
c     add in contribution from the boundary.
c
      do 400 i=2,mp
         f(i,2)  =f(i,2)    +d1*f(i,1)  -twodel*
     +           (f(i+1,1)  +f(i-1,1))  -twody*bdc(i-1)
         f(i,3)  =f(i,3)    -f(i,1)
         f(i,n+1)=f(i,n+1)  +d1*f(i,n+2)-twodel*
     +           (f(i+1,n+2)+f(i-1,n+2))-twody*bdd(i-1)
         f(i,n)  =f(i,n)    -f(i,n+2)
  400    continue
      do 500 j=2,np
         f(2,j)  =f(2,j)    +d2*f(1,j)  -twodel*
     +           (f(1,j+1)  +f(1,j-1))  -d3*bda(j-1)
         f(3,j)  =f(3,j)    -del2*f(1,j)
         f(m+1,j)=f(m+1,j)  +d2*f(m+2,j)-twodel*
     +           (f(m+2,j+1)+f(m+2,j-1))-d3*bdb(j-1)
         f(m,j)  =f(m,j)    -del2*f(m+2,j)
  500        continue
      f(2,2)     =f(2,2)    +twodel*f(1,1)
      f(m+1,2)   =f(m+1,2)  +twodel*f(m+2,1)
      f(2,n+1)   =f(2,n+1)  +twodel*f(1,n+2)
      f(m+1,n+1) =f(m+1,n+1)+twodel*f(m+2,n+2)
      return
      end
