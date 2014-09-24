      subroutine sbiplr(a,b,m,bda,bdb,n,f,idf,
     +                  alpha,beta,iflag,w,lw)
c
      integer           m,n,idf,iflag,lw
      real  a,b,alpha,beta
      real  bda(*),bdb(*),f(idf,*),w(*)
c
c     this subroutine solves the generalized biharmonic
c     equation
c
c
c     1           1  2  1            1  2
c     - ( d r d + - d ) - ( d r d  + - d ) u(r,t)  +
c     r    r   r  r  t  r    r   r   r  t
c
c            1           1  2
c     alpha (- ( d r d + - d ) u(r,t)  + beta u(r,t)  = f(r,t)
c            r    r   r  r  t
c
c
c     in polar coordinates. the boundary conditions are of the first
c     kind, u and the first derivative of u with respect to r must
c     be specified on the boundary.
c     r is radial coordinate, while t denotes the theta coordinate.
c     the equation is solved in a region between two concentric
c     circles. the radius of the inner circle being a and that
c     of the outer circle being b.
c     the special case when the region
c     is a disk is accepted when a is set to zero.
c
c     the solution time is proportional to m*n*log(n).
c
c     this code is written by petter bjorstad, the author
c     acknowledges paul n swarztrauber for the use of his
c     fast fourier transform package and the linpack project
c     for the use of linear equation solvers.
c
c     this is version 2.0 , single precision. january 1984.
c
c     any comments, questions or suggestions are most welcome
c     and should be directed to :
c
c                  petter bjorstad
c
c                  veritas research
c                  p.o. box 300
c                  n-1322 hovik
c                  norway
c
c     description of input variables.
c
c       a,b          defines the domain where the equation is defined.
c                    see explanation above.
c                    note the special meaning if a is set to zero.
c
c       m            the number of interior gridpoints in the
c                    r direction. the interval a .le. x .le. b is divide
c                    in m+1 panels each of width (b-a)/(m+1).
c                    m must be at least 3.
c
c       n            the number of  gridpoints in the theta direction.
c                    the angle 2 pi (6.28...) is divided into n panels
c                    each of width 2*pi/n. n should be at least 3.
c                    the method is more efficient if n is a product of
c                    small primes.
c
c       alpha        the constant alpha in the equation.
c
c       beta         the constant beta in the equation.
c
c       bda          array of dimension at least n.
c                    the dimension must be at least m if a is zero.
c                    if a is non-zero, then
c                    bda(j) contains the derivative with respect
c                    to r of the function at the inner boundary.
c                    (at r = a and t = j*2*pi/n  ,j=1,2,..n.)
c                    if the parameter a is zero then this array
c                    of length at least m is used as an extra
c                    workspace.
c
c       bdb          array of dimension at least n.
c                    bdb(j) contains the derivative with respect
c                    to r of the function at the outer boundary.
c                    (at r = b and t = j*2*pi/n  ,j=1,2,..n.)
c
c       f            array of dimension at least ( m+2,n ) .
c                    f(i,j) must contain the values of the right
c                    hand side function f(r,t) evaluated at the points
c                    (r ,t ) , where
c                      i  j
c                    r = a+(i-1)*(b-a)/(m+1)  , i=2,3,.. m+1.
c                     i
c                    t = j*2*pi/n             , j=1,2,.. n.
c                     j
c                    f(1,j)  contains the boundary values at the
c                            inner boundary. (r =a , t = t  )
c                                                         j
c                    if a is zero then f(1,j) should contain the value
c                    of the function f(0,t). (the element f(1,1) will
c                    be read by the code).
c
c                    f(m+2,j)  contains the boundary values at the
c                              outer boundary. (r = b , t = t )
c                                                            j
c
c       idf          rowdimension of f as declared in the
c                    calling program.
c
c       iflag   = 1  if the discrete system is positive definite.
c                    this is always the case if alpha is less or equal
c                    zero and beta is greater or equal zero. this option
c                    is somewhat faster than using iflag = 2.
c               = 2  this option should handle all nonsingular discrete
c                    versions of the equations.
c
c       description of output variables.
c
c       f            contains the solution  u of the discrete
c                    system in all gridpoints (r ,t ) , i=1,2,.. m+2.
c                                               i  j  , j=1,2,.. n.
c
c       iflag        error indicator.
c               = 0  normal return.
c               =-1  n and/or m was less than 3.
c               =-2  a is greater or equal b or a is negative.
c               =-3  idf is less than m+2 or lw is too small.
c               =-4  linpack failure in cholesky factorization.
c                    this can happen if the routine is called with
c                    iflag=1 and the linear system is indefinite.
c                    try calling with iflag=2.
c                    this error should never occur if alpha is less or
c                    equal zero and beta is greater or equal zero.
c               =-5  linpack detected a computationally singular
c                    system using gaussian elimination with pivoting.
c
c       description of workspace.
c
c       lw           the length of the user supplied workspace w.
c                    lw depends on iflag in the following way:
c
c                    iflag          lw must be at least
c
c                      1            2*m+n+max(8*m+4,2*n+15)
c                      2            2*m+n+max(13*m,2*n+15)
c
c       w            a one dimensional array containing at least lw
c                    elements. w can be used for other purposes
c                    between calls to sbiplr.
c
c       subroutines and functions
c       needed by this program.
c
c            1. biharmonic:   sbipl
c            2. fourier
c               transform :   srffti, srfti1, srfftf, srftf1,
c                             sradf2, sradf3, sradf4, sradf5,
c                             sradfg, srfftb, drftb1, sradb2,
c                             sradb3, sradb4, sradb5, sradbg.
c
c            3. linpack   :   spbfa, spbsl, sgbfa, sgbsl.
c
c            4. blas      :   sscal, scopy, sdot, saxpy, isamax
c
c            5. fortran   :   min, max, mod, cos, atan, sqrt
c                             sin
c
c
c      local.
c
c      biharmonic:   sbipl
c      fortran   :   max
c
      integer i1,i2,i3,i4,i5,i6,i7,i8
c
c     check for input errors.
c
      i1 = 2*m + n + max(8*m+4, 2*n+15)
      i2 = 2*m + n + max(13*m , 2*n+15)
c
      if(n     .lt. 3 .or. m .lt. 3)    iflag = -1
      if(a     .ge. b)                  iflag = -2
      if(a     .lt. 0.0e0)              iflag = -2
      if(idf   .lt. m+2)                iflag = -3
      if(iflag .eq. 1 .and. lw .lt. i1) iflag = -3
      if(iflag .eq. 2 .and. lw .lt. i2) iflag = -3
      if(iflag .lt. 0)                  go to 100
c
c     assign workspace.
c
      i1 = 1
      i2 = m+i1
      i3 = m+i2
      i4 = n+i3
      i5 = m+i4
      i6 = m+i5
      i7 = m+i6
      i8 = m+i7
c
      call sbipl(a,b,m,bda,bdb,n,f,idf,alpha,beta,iflag,
     +           w(i1),w(i2),w(i3),w(i4),w(i5),w(i6),
     +           w(i7),w(i8))
c
      if(iflag.ge.0) return
100   write(6,1) iflag
    1 format(1x,'error return from sbiplr , iflag= ',i4)
      end
