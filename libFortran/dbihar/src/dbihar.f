      subroutine dbihar(a,b,m,bda,bdb,bdc,bdd,c,d,n,f,idf,
     +                  alpha,beta,iflag,tol,itcg,w,lw)
c
      integer m,n,idf,iflag,itcg,lw
      double precision a,b,c,d,alpha,beta,tol
      double precision bda(*),bdb(*),bdc(*),bdd(*),f(idf,*),w(*)
c
c
c     this subroutine solves the equation
c
c     u    + 2*u    + u    + alpha*(u  +u  ) + beta*u = f(x,y)
c      xxxx     xxyy   yyyy         xx  yy
c
c
c     in the rectangle   a < x < b and c < y < d  ,
c
c     subject to boundary conditions of the first
c     kind, u and the exterior normal derivative
c     of u must be specified at the boundary.
c
c     the solution time is essentially proportional to the
c     number of gridpoints when the conjugate gradient
c     option is used.
c
c     this is version 2.0, april 1984.
c
c     this code is written by petter bjorstad, the author
c     acknowledges paul n swarztrauber for the use of his
c     fast fourier transform package and the linpack project
c     for the use of linear equation solvers.
c
c     this code and the mathematical theory behind it, is
c     described in detail in:
c
c          numerical solution of the biharmonic equation.
c
c          ph.d. dissertation, stanford university  1980.
c          copyright 1980, petter e. bjorstad.
c
c
c     any comments, questions or suggestions are most welcome
c     and should be directed to :
c
c                      petter e. bjorstad
c
c                      veritas research
c                      p.o. box 300
c                      n-1322 hovik
c                      norway
c
c     description of input variables.
c
c     a,b,c and d.    defines the rectangle.
c                     see above.
c
c     m               the number of interior grid points
c                     in the x-direction. the interval
c                     a.le.x.le.b is divided in (m+1) panels
c                     each of width (b-a)/(m+1).
c                     m must be odd and at least 3.
c                     the method is more efficient
c                     if m+1 is a product of small primes.
c
c     n               the number of interior grid points
c                     in the y-direction. the interval
c                     c.le.y.le.d is divided in (n+1) panels
c                     each of width (d-c)/(n+1).
c                     n must be odd and at least 3. the method is
c                     somewhat faster if n+1 is a product of small
c                     primes. it is also advisable to take n less
c                     than or equal to m.
c                     (for max speed and min storage.)
c
c     alpha           the constant alpha in the equation.
c
c     beta            the constant beta in the equation.
c
c     bda             array of dimension n containing the values
c                     of the exterior normal derivative on the
c                     side x=a ,y=c+j*(d-c)/(n+1)  ,j=1,..n .
c
c     bdb             array of dimension n containing the values
c                     of the exterior normal derivative on the
c                     side x=b ,y=c+j*(d-c)/(n+1)  ,j=1,..n .
c
c     bdc             array of dimension m containing the values
c                     of the exterior normal derivative on the
c                     side y=c ,x=a+i*(b-a)/(m+1)  ,i=1,..m .
c
c     bdd             array of dimension m containing the values
c                     of the exterior normal derivative on the
c                     side y=d ,x=a+i*(b-a)/(m+1)  ,i=1,..m .
c
c     f               array of dimension at least (m+2,n+2).
c                     f(i,j) must contain the values of the right
c                     hand side function f(x,y) evaluated at the
c                     points (x ,y ), where
c                              i  j
c                     x =a+(i-1)*(b-a)/(m+1) ,i=2,..m+1 .
c                      i
c                     y =c+(j-1)*(d-c)/(n+1) ,j=2,..n+1 .
c                      j
c
c                     f(1,j)   must contain the boundary values along
c                     the side x=a, y=y  ,j=2,..n+1 .
c                                      j
c                     f(m+2,j) must contain the boundary values along
c                     the side x=b, y=y  ,j=2,..n+1 .
c                                      j
c                     f(i,1)   must contain the boundary values along
c                     the side y=c, x=x  ,i=2,..m+1 .
c                                      i
c                     f(i,n+2) must contain the boundary values along
c                     the side y=d, x=x  ,i=2,..m+1 .
c                                      i
c
c     idf             rowdimension of array f as declared in
c                     the calling program.
c
c     tol             the conjugate gradient iteration is
c                     terminated using this parameter.
c                     the error in the solution of the discrete
c                     approximation (not to be confused with the
c                     truncation error) can be expected to be of the
c                     same order of magnitude provided that the
c                     function f has norm of magnitude unity.
c                     a good choice for tol is taking it a few
c                     orders of magnitude less than the expected
c                     truncation error.
c                     if the right hand side f in the given problem
c                     is many orders of magnitude smaller or bigger
c                     than one, then it is recommended that the user
c                     scales his problem.
c                     tol is a dummy variable if iflag is 3 or 4.
c
c     iflag
c              =1     this option is not available in version 2.0 .
c
c              =2     if this is the first solution on this grid
c                     using conjugate gradients with fourier trans-
c                     forms in both coordinate directions.
c                     there are no parameter restrictions, but
c                     the routine is only guaranteed to work when
c                     the discrete approximation is positive
c                     definite. this is always the case when alpha
c                     is nonpositive and beta is nonnegative.
c                     in other cases the user should monitor the
c                     output parameter itcg. if it stays larger
c                     than 15, the alternative iflag=4 is
c                     recommended.
c                     a call with iflag=2 restarts the low rank
c                     approximations used to speed up convergence.
c
c              =3     if this is the first solution on this grid
c                     using cholesky factorization.
c                     the parameters must satisfy alpha.le.0 ,
c                     beta.ge.0. (if this is violated the code
c                     will change iflag to 4.)
c
c              =4     if this is the first solution on this grid
c                     using an indefinite symmetric factorization.
c                     there are no parameter restrictions, but
c                     an error return will occur if the discrete
c                     system is computationally singular.
c
c     description of output variables.
c
c     f               contains the solution u of the discrete
c                     system in all gridpoints (x ,y ),i=1,..m+2 ,
c                                                i  j  j=1,..n+2 .
c
c     iflag
c
c     normal returns.  iflag will return with a value appropriate
c                      for repeated calls. it need normally not
c                      be changed by the user between calls to
c                      dbihar.
c                      however, if a sequence of problems is being
c                      solved where the parameters alpha or beta
c                      changes then itcg should be monitored and
c                      iflag should be reset to 2 when itcg
c                      increases (say above 15).
c                      also a change between the direct and
c                      iterative version must be set explicitly.
c
c                      if the code changes the value of iflag due
c                      to a change in the user input, a warning is
c                      printed indicating the reason for the change.
c
c                      iflag=2 new initial solution (see above).
c                      iflag=3 new initial solution (see above).
c                      iflag=4 new initial solution (see above).
c                      iflag=6 repeated solution after iflag=2.
c                      iflag=7 repeated solution after iflag=3.
c                      iflag=8 repeated solution after iflag=4.
c
c     error returns.   if iflag returns a negative value then
c                      an error was detected. the computed f(i,j)
c                      should be considered wrong. an error message
c                      giving the value of iflag is printed.
c
c              =-1     n and/or m was even or less than 3.
c              =-2     a.ge.b and/or c.ge.d .
c              =-3     idf.lt.m+2 or lw is too small.
c              =-4     linpack failure in cholesky-factorization.
c                      this should not occur,check input carefully.
c              =-5     linpack detected a computationally singular
c                      system using the symmetric indefinite
c                      factorization.
c              =-6     the conjugate gradient iteration failed to
c                      converge in 30 iterations. the probable
c                      cause is an indefinite or near singular
c                      system. try using iflag=4. note that tol
c                      returns an estimate of the residual in
c                      the current conjugate gradient iteration.
c
c     tol              an upper bound on the residuals obtained in
c                      the conjugate gradient iterations.
c                      tol will therefore normally be unchanged.
c
c     itcg             the number of conjugate gradient iterations.
c                      if this is large (say 20) then a direct
c                      solution using iflag =4 may be considered.
c
c     description of workspace.
c
c     lw              integer indicating the number of elements
c                     supplied in the workspace w.
c                     lw depends on iflag in the following way:
c                     (this can be improved slightly,version 3.0 ??)
c
c          iflag=  2: lw must be at least max(7*n,3*m)+2*(n+m)+19.
c                     if only one problem is solved on the grid then
c                     lw should be given its minimum value. (for
c                     maximum speed.)
c                     if several problems are to be solved then
c                     any larger lw will reduce the execution
c                     time for subsequent problems.(a low rank
c                     correction will be computed and used to
c                     improve the preconditioning that is used.)
c                     the code will not make use of lw larger than
c                     max(7*n,3*m)+2*(n+m)+19+20*(n+3) under any
c                     circumstance.
c
c         iflag=3,4:  lw must be at least max(3*m,4*n)+4*n+2*m
c                                        +0.5*(n+1)**2+19  .
c
c
c     w               a one dimensional array with at least
c                     lw elements.
c
c     subroutines and functions
c     included in this package.
c
c            1. biharmonic:    dbihar, dstart, dftrnx, dftrny,
c                              dbislf, dbisld, dpentf, dmatge,
c                              dtrigi, dhzeri, dhzero, dconju,
c                              dcmult, dpreco, dupdat.
c
c            2. fourier        dsinti, drffti, drfti1,
c               transform:     dsint,  drfftf, drftf1, dradf2,
c                              dradf3, dradf4, dradf5, dradfg.
c                              (by p.n. swarztrauber.)
c
c            3. linpack:       dspfa, dspsl, dppfa, dppsl.
c
c            4. blas:          idamax, daxpy, dcopy, ddot,
c                              dscal, dswap.
c
c            5. fortran:       mod,min,max,abs,
c                              atan,cos,sin,sqrt.
c
c     local.
c
c     biharmonic:         dstart, dftrnx, dftrny, dbislf, dbisld
c     fortran:            mod,min,max
c
      integer i1,i2,i3,i4,i5,i6,i7,i8
      integer nold,mold,iwf,iwl,il1,il2
      integer maxi
      double precision dx,dy,del,alf,bet
      double precision dxo,dyo,alfo,beto,wfo,wlo
c
c     check input for mistakes.
c
      il1  = max(7*n,3*m)+2*(n+m)
      il2  = max(4*n,3*m)+4*n+2*m+(n+1)**2/2+19
      if(n.lt.3.or.m.lt.3)               iflag = -1
      if(mod(m,2).eq.0.or.mod(n,2).eq.0) iflag = -1
      if(a.ge.b.or.c.ge.d)               iflag = -2
      if(idf.lt.m+2.or.lw.lt.il1)        iflag = -3
      if(iflag.lt.0)                     go to 1000
c
      dx   = (b-a)/(m+1.0d0)
      dy   = (d-c)/(n+1.0d0)
      del  = (dy/dx)**2
      alf  = alpha*dy*dy
      bet  = beta*dy**4
      maxi = (lw-il1)/(2*n+6)
      maxi = min(maxi,10)
c
      call dstart(m,n,alf,bet,f,idf,bda,bdb,bdc,bdd,dx,dy,del)
      call dftrnx(m,n,f(2,2),idf,w)
c
      if(iflag.eq.3.and.lw.lt.il2) go to 10
      if(iflag.le.4) go to 40
      if(n .ne.nold) go to 20
      if(m .ne.mold) go to 20
      if(dx.ne.dxo)  go to 20
      if(dy.ne.dyo)  go to 20
      if(w(iwf).ne.wfo.or.w(iwl).ne.wlo)  go to 30
      if(iflag.le.6) go to 50
      if(alfo.eq.alpha.and.beto.eq.beta) go to 50
      iflag= iflag-4
      write(6,2) iflag,alpha,beta
      go to 40
   10 iflag= 2
      write(6,5) iflag
      go to 40
   20 iflag= iflag-4
      write(6,1) iflag,n,m,a,b,c,d
      go to 40
   30 iflag= iflag-4
      write(6,3) iflag,iwf,iwl,wfo,wlo,w(iwf),w(iwl)
   40 nold = n
      mold = m
      dxo  = dx
      dyo  = dy
      alfo = alpha
      beto = beta
   50 go to (60,70,80,90,100,110,120,130),iflag
      go to 1000
   60 iwf  = max(4*n,3*m)+1
      iwl  = il1+(n+3)*maxi
      iflag= 2
      write(6,2) iflag,alpha,beta
   70 iwf  = max(4*n,3*m)+1
      iwl  = il1+(n+3)*maxi
      go to 200
   80 iwf  = max(3*m,4*n)+1
      iwl  = il2
      if(alpha.le.0.0d0.and.beta.ge.0.0d0) go to 200
      iflag= 4
      write(6,2) iflag,alpha,beta
   90 iwf  = max(3*m,4*n)+1
      iwl  = il2
      go to 200
  100 iflag= 2
      write(6,2) iflag,alpha,beta
  110 go to 200
  120 go to 200
  130 go to 200
  200 call dftrny(m,n,f(2,2),idf,w)
      if(iflag.eq.6) go to 250
      if(iflag.ne.2) go to 300
      i1   = 1
      i2   = i1+(n+1)/2
      i3   = i2+(n+1)/2
      i4   = i3+(n+1)/2
      i5   = i4+(n+1)/2
      i6   = max((5*m)/2,(7*n)/2)+17
      i7   = i6+2*(m+n)
      i8   = i7+2*maxi*(n+3)
  250 call dbislf(m,n,maxi,iflag,del,tol,alf,bet,itcg,idf,f(2,2),
     +     w(i1),w(i2),w(i3),w(i4),w(i5),w(i6),w(i7),w(i8))
      if(iflag.lt.0) go to 1000
      if(iflag.eq.2) iflag = 6
      go to 400
  300 if(iflag.eq.7.or.iflag.eq.8) go to 350
      i1   = 1
      i2   = i1+(n+1)/2
      i3   = i2+(n+1)/2
      i4   = max((5*m)/2,(7*n)/2)+17
      i5   = i4+2*(m+n)
  350 call dbisld(m,n,iflag,del,alf,bet,idf,f(2,2),
     +            w(i1),w(i2),w(i3),w(i4),w(i5))
      if(iflag.lt.0) go to 1000
      if(iflag.eq.3) iflag = 7
      if(iflag.eq.4) iflag = 8
  400 call dftrny(m,n,f(2,2),idf,w)
      call dftrnx(m,n,f(2,2),idf,w)
      wfo  = w(iwf)
      wlo  = w(iwl)
      return
 1000 write(6,4) iflag
    1 format(1x,'*warning*,iflag changed to ',i3,'n,m,maxi,a,b,c,d=',
     +       2i6,2x,4e12.2)
    2 format(1x,'*warning*,iflag changed to ',i3,'alpha,beta=',2e16.6)
    3 format(1x,'*warning*,iflag changed to ',i3,/1x,'element no ',
     +       i6,' and ',i6,' of w changed from ',2e16.6,' to ',
     +       2e16.6,' by user. ')
    4 format(/5x,'***error in dbihar, iflag= ',i6/)
    5 format(1x,'*warning*,iflag changed to ',i3,
     +       ' workspace needed, given : ',2i8)
      return
      end
