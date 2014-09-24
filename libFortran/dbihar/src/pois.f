      subroutine spois(m,n,idf,f,del,tm,tn,w)
c
      integer m,n,idf
      real    del
      real    f(idf,*),tm(*),tn(*),w(*)
c
c      this routine solves the second order five point discrete
c      poissons equation using fourier transforms.
c       the boundary condition is assumed to be dirichlet.
c
c       petter bjorstad mai 1980.
c
c      input.
c      m,n	number of interior gridpoints in
c      	x and y directions. m and n must both be odd.
c      del	meshratio (dy/dx).
c      idf 	declared rowdimension of array f.
c      tm	array containing sin(.5*i*pi/(m+1)) , i=1,2,..m.
c      tn	array containing sin(.5*j*pi/(n+1)) , j=1,2,..n.
c       f       array containing the appropriate right hand
c               side. (the discrete forcing function and
c               the contribution from the boundary condition.)
c               this can be defined by calling the routine spoirh.
c
c      output.
c      f	containing the solution.
c      	all other input variables are unchanged.
c
c      workspace.
c      w	having at least int(3.5*max(m,n)+16) elements.
c
c      note that the current version wastes space and time
c      by recomputing and double storing trigonometric
c      information. (in the fourier transform part)
c       this routine is not designed to be very fast , but it
c       does have an operationcount of order n*m*max(log(n),log(m)).
c
c       routines called:      sftrnx,sftrny
c       fortran:              float
c
c      local variables.
c
      integer i,j
      real    scal1,scal2
c
      call sftrnx(m,n,f,idf,w)
      call sftrny(m,n,f,idf,w)
      i     = 16*(m+1)*(n+1)
      scal1 = float(i)
      scal2 = del*del*scal1
      do 10 i=1,m
      	w(i) = scal2*tm(i)*tm(i)
   10      	continue
      do 20 j=1,n
      	w(m+j) = scal1*tn(j)*tn(j)
   20      	continue
      do 30 j=1,n
      do 30 i=1,m
      	f(i,j) = f(i,j)/(w(i)+w(m+j))
   30      	continue
      call sftrny(m,n,f,idf,w)
      call sftrnx(m,n,f,idf,w)
      return
      end
      subroutine spoirh(m,n,idf,f,dx,dy,xmin,ymin)
c
      integer m,n,idf,nr
      real    dx,dy,xmin,ymin
      real    f(idf,*)
c
c      sets up the right hand side (from data f and g)
c      for the discrete poissons equation. this routine
c      is often called before the routine spois.
c       the user must supply functions func(x,y) and g(x,y).
c       func(x,y) computes the forcingfunction at an interor
c                 gridpoint. (corresponding to (x,y)).
c       g(x,y)    computes the value of the boundary function
c                 at a given boundary gridpoint (x,y).
c
c      input.
c      m,n       	number of gridpoints in x and y direction.
c       idf             declared rowdimension of f.
c      dx,dy		gridspacing in x and y direction.
c       xmin,ymin       x at west edge,y at south edge.
c
c       output.
c
c       f               array containing the appropriate right
c                       hand side.
c
c       user supplied functions:    func,g
c
c
c       local variables.
c
      integer i,j
      real    x,y,dy2,del2,x1,y1
c
      dy2 = dy*dy
        del2= (dy/dx)**2
      x1  = xmin+(m+1)*dx
      y1  = ymin+(n+1)*dy
      do 10 j=1,n
      	y = ymin+dy*j
      	do 10 i=1,m
      		x = xmin+dx*i
      		f(i,j) = -dy2*func(x,y)
   10                   continue
      do 20 j=1,n
      	y = ymin+dy*j
      	f(1,j) = f(1,j)+g(xmin,y)*del2
      	f(m,j) = f(m,j)+g(x1,y)*del2
   20           continue
        do 30 i=1,m
      	x = xmin+dx*i
      	f(i,1) = f(i,1)+g(x,ymin)
      	f(i,n) = f(i,n)+g(x,y1)
   30           continue
      return
      end
      subroutine strig(n,tn)
c
      integer n
      real    tn(*)
c
c      this routine fills in the array tn(n) with
c      trigonometric information needed in spois.
c      (this is a very slow version.)
c
c      tn(i)=sin(.5*i*pi/(n+1))  i=1,2,.. n.
c
c       fortran:         atan,sin,cos,float
c
c       local variables.
c
      integer i
      real    pi,arg
c
      pi  = 4.e0*atan(1.e0)
      arg = .5e0*pi/float(n+1)
      do 10 i = 1,n
      	tn(i) = sin(arg*i)
   10           continue
      return
      end
      subroutine sftrnx(m,n,f,idf,ws)
c
      integer m,n,idf
      real    f(idf,*)
      real    ws(*)
c
c     performs n sine transforms. the transform is unscaled.
c     m must be odd.(when ssint is called.)
c     the workspace ws must have at least int(3.5*m+16) elements.
c     the call to scopy is performed since ssint uses a vector
c     of lenght m+1. (if f(m+1,n) is a legal address then one
c     can save and restore this value instead.)
c
c     fourier:            ssinti,ssint   (swarztrauber version 3.)
c     blas:               scopy
c
      integer j
      real    x1
c
      call ssinti(m,ws(m+2))
      do 500 j=1,n
         call scopy(m,f(1,j),1,ws,1)
         call ssint(m,ws,ws(m+2))
         call scopy(m,ws,1,f(1,j),1)
  500    continue
      return
      end
      subroutine sftrny(m,n,f,idf,ws)
c
      integer m,n,idf
      real    f(idf,*)
      real    ws(*)
c
c     performs m sine transforms. the transform is unscaled.
c     n should be odd. (when ssint is called.)
c     workspace ws must have at least int(3.5*n+16) elements.
c
c     fourier:            ssinti,ssint   (swarztrauber version 3.)
c     blas:               scopy
c
      integer i,j
c
      call ssinti(n,ws(n+2))
      do 500 i   =1,m
         call scopy(n,f(i,1),idf,ws,1)
         call ssint(n,ws,ws(n+2))
         call scopy(n,ws,1,f(i,1),idf)
  500    continue
      return
      end
