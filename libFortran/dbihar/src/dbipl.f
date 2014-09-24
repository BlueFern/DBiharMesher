      subroutine dbipl(a,b,m,bda,bdb,n,f,idf,alpha,beta,iflag,
     +                 d1,amat,ev,a1,a2,a3,d2,w)
c
      integer           m,n,idf,iflag
      double precision  a,b,alpha,beta
      double precision  bda(*),bdb(*),f(idf,*)
      double precision  amat(*),ev(*),a1(*),a2(*),a3(*),w(*)
      double precision  d1(*),d2(*)
c
c     biharmonic on a disk.
c     see description in routine dbiplr.
c
c     petter bjorstad january 1984.
c
c     local.
c
c     fourier transform : drffti, drfftf, drfftb
c     linpack           : dpbfa, dpbsl, dgbfa, dgbsl
c     blas              : dcopy, dscal
c     fortran           : mod, atan, cos, sqrt
c
      integer           i,j,jp,jm,nh,mm,info
      double precision  dr,drh,dt,dt4,alf,bet,del2,del4
      double precision  tpi,zero(1)
      double precision  t1,t2,t3,t4,t5,t6,t7,t8,t9
c
      zero(1) = 0.0d0
      nh      = (n-1)/2
      mm      = m-1
      tpi     = 8.0d0*atan(1.0d0)
      dr      = (b-a)/(m+1)
      drh     = .50d0*dr
      dt      = tpi/n
      dt4     = dt**4
      del2    = (dt/dr)**2
      del4    = del2*del2
      alf     = dt*dt*alpha
      bet     = dt4*beta
c
c     set up the correct right hand side vector,
c     incorporating the boundary data.
c
      do 20 j=1,n
         call dscal(m,dt4,f(2,j),1)
   20    continue
c
      if(a.eq.0.0d0) go to 30
      t1 = del2*(a+drh)*(4.0d0*del2+2.0d0/
     +     (a+dr)**2+2.0d0/a**2-alf)/(a+dr)
      t2 = del2*(a+drh)*(1.0d0/(a+dr)**2+1.0d0/a**2)/(a+dr)
      t3 = 2.0d0*dr*del4*(a+drh)*(a-drh)/(a*(a+dr))
      t4 = del4*(a+1.50d0*dr)*(a+.50d0*dr)/((a+2.0d0*dr)*(a+dr))
      go to 40
   30 f(1,1) = dt4 * f(1,1)
      call dcopy(m-2,zero,0,bda(3),1)
   40 t5 = del2*(b-drh)*(4.0d0*del2+2.0d0/
     +     (b-dr)**2+2.0d0/b**2-alf)/(b-dr)
      t6 = del2*(b-drh)*(1.0d0/(b-dr)**2+1.0d0/b**2)/(b-dr)
      t7 = 2.0d0*dr*del4*(b-drh)*(b+drh)/(b*(b-dr))
      t8 = del4*(b-1.50d0*dr)*(b-.50d0*dr)/((b-2.0d0*dr)*(b-dr))
      do 50 j=1,n
         jp = mod(j,n)+1
         jm = mod(j-2+n,n)+1
         if(a.eq.0.0d0) go to 45
         f(2,j) = f(2,j)+t1*f(1,j)-t2*(f(1,jp)+f(1,jm))+t3*bda(j)
         f(3,j) = f(3,j)- t4*f(1,j)
   45    f(m+1,j) = f(m+1,j)+t5*f(m+2,j)-t6*(f(m+2,jp)+f(m+2,jm))-
     +              t7*bdb(j)
         f(m,j) = f(m,j)-t8*f(m+2,j)
   50    continue
c
c     set up elementary matrices and scale f.
c     this is the first multiply of f by (i tensorprod d**(-1/2)).
c
c     d1 is the diagonal matrix d from the theory.
c     amat(i) is the scalar ( a sub i ) from the theory.
c     amat is scaled by del2.
c
      t1 = a
      do 100 i=1,m
         t1   = t1+dr
         t2   = dr/t1
         t3   = sqrt(t1)
         d1(i)= 1.0d0/t1
         amat(i) = del2 * (1.0d0+.50d0*t2)/sqrt(1.0d0+t2)
         call dscal(n,t3,f(i+1,1),idf)
  100    continue
c
c     forward real transform.
c
      call drffti(n,a1)
      do 150 i=1,m
         call dcopy(n,f(i+1,1),idf,ev,1)
         call drfftf(n,ev,a1)
         call dcopy(n,ev,1,f(i+1,1),idf)
  150    continue
c
c     compute eigenvalues of periodic tridiag(-1,2,-1).
c     the matrix r in the theory.
c
      ev(1) = 0.0d0
      do 120 j=1,nh
         ev(2*j)  = 2.0d0*(1.0d0-cos(tpi*j/n))
  120    continue
      call dcopy(nh,ev(2),2,ev(3),2)
      if(mod(n,2).eq.0) ev(n) = 4.0d0
c
c     precompute parts of the pentadiagonal coefficients.
c
      if(a.gt.0.0d0) then
         a1(1) = del4*4.0d0+amat(1)**2+del4*2.0d0*(a+drh)/(a+dr)+bet
      else
         a1(1) = del4*4.0d0+amat(1)**2+bet
      end if
      a1(m) = del4*4.0d0+amat(m-1)**2+del4*2.0d0*(b-drh)/(b-dr)+bet
c
      do 200 i=2,mm
         a1(i) = del4*4.0d0+amat(i-1)**2+amat(i)**2+bet
  200    continue
       do 210 i=1,mm
         a2(i)    = - 4.0d0*del2*amat(i)
         a3(i)    =   amat(i)*amat(i+1)
  210    continue
      do 220 i=1,m
         d2(i) = d1(i)**2
  220    continue
c
      t5 = 2.0d0*del2*alf
      t9 = 4.0d0 * del2 - alf
c
      do  400 j=1,n
         t6 = ev(j)
         if(iflag.eq.2) go to 350
         t8 = d2(1)*t6
         do 300 i=1,m
            t7 = t8
            t8 = d2(i+1) * t6
            w(3*i)   = a1(i) + t7 * (t9 + t7)-t5
            w(3*i+2) = a2(i) - amat(i) * ((t7 + t8) - alf)
            w(3*i+4) = a3(i)
  300       continue
         if(j .eq. 1 .and. a .eq. 0.0d0) w(3) = w(3) + 2.0d0*del4
         call dpbfa(w,3,m,2,info)
         if(info.ne.0) go to 900
         call dpbsl(w,3,m,2,f(2,j))
c
c    if we have a disk, then solve special system as well.
c
         if(j.eq.1 .and. a.eq.0.0d0) then
            bda(1)= n*(-3.0d0*del4+.50d0*alf*del2)/sqrt(d1(1))
            bda(2)= n*(del4*3.0d0/8.0d0)/sqrt(d1(2))
            call dpbsl(w,3,m,2,bda)
         end if
c
c     end loop if positive definite system.
c
         go to 400
c
c     start indefinite solution path.
c
  350    t8 = d2(1) * t6
         do 360 i=1,m
            t7 = t8
            t8 = d2(i+1)*t6
            w(7*i-2) = a1(i)+t7*(4.0d0*del2+t7-alf)-t5
            w(7*i-1) = a2(i)-amat(i) * ((t7+t8)-alf)
            w(7*i)   = a3(i)
  360       continue
         if(j.eq.1 .and. a.eq.0.0d0) w(5) = w(5) + 2.0d0*del4
         call dcopy(m-1,w(6),7,w(11),7)
         call dcopy(m-2,w(7),7,w(17),7)
         call dgbfa(w,7,m,2,2,w(7*m+1),info)
         if(info.ne.0) go to 910
         call dgbsl(w,7,m,2,2,w(7*m+1),f(2,j),0)
c
c    check if we have a disk
c
         if(j.eq.1 .and. a.eq.0.0d0) then
            bda(1)= n*(-3.0d0*del4+.50d0*alf*del2)/sqrt(d1(1))
            bda(2)= n*(del4*3.0d0/8.0d0)/sqrt(d1(2))
            call dgbsl(w,7,m,2,2,w(7*m+1),bda,0)
         end if
c
c     end loop 400
c
  400    continue
c
c     inverse real transform.
c
      call drffti(n,a1)
      do 450 i=1,m
         call dcopy(n,f(i+1,1),idf,ev,1)
         call drfftb(n,ev,a1)
         call dcopy(n,ev,1,f(i+1,1),idf)
  450    continue
      do 550 i=1,m
         t1 = sqrt(d1(i))/n
         call dscal(n,t1,f(i+1,1),idf)
  550    continue
c
c     if a is greater than zero then exit.
c
      if(a.gt.0.0d0) return
c
c     if we have a=0, ie if the equation is solved in a
c     disk, then update the solution.
c
c    scale special solution vector by sqrt(d1)/n.
c
      do 600 i=1,m
         bda(i) = bda(i)*sqrt(d1(i))/n
  600    continue
c
c    the two first components of y are t1 and t2,
c    the two first components of x are t3 and t4.
c
           t1 = (4.0d0*del2*alf- 64.0d0*del4/3.0d0)/n
           t2 = (16.0d0*del4/3.0d0)/n
           t3 = -3.0d0*del4 + .50d0*alf*del2
           t4 = 3.0d0*del4/8.0d0
c
c    we can now compute the scalar y transpose * z
c    (z is stored in bda).
c
      t5 = n*(bda(1)*t1 + bda(2)*t2)
c
c    the inner product between y and f
c
      t6 = 0.0d0
      do 650 j=1,n
         t6=t6+t1*f(2,j)+t2*f(3,j)
  650    continue
c
c    we now compute the value of f at the center point.
c
      t7 = 16.0d0 * del4 - 4.0d0 * del2 * alf + bet
      f(1,1) = (f(1,1) - t6)/(t7 - t5)
c
c    copy the center point solution to all of f(1,j).
c
      call dcopy(n-1,f(1,1),0,f(1,2),idf)
c
c    update the rest of f.
c
      do 700 j = 1,n
         call daxpy(m,-f(1,1),bda,1,f(2,j),1)
  700    continue
      return
c
c     error returns
c
  900 continue
      iflag=-4
      return
  910 continue
      iflag=-5
      return
      end
