      subroutine strigi(n,del,trig,w)
c
      integer n
      real    del
      real    trig(*),w(*)
c
c     this routine computes trigonometric information needed to
c     solve the biharmonic equation.
c
c     input.
c            n   is assumed to be odd and at least 3.
c            del is (dy/dx)**2.
c     output.
c            trig(i), i=1,2.. 2*n.
c
c     workspace.
c            w  must have at least n/2+(n/2+1)/2 elements.
c
c            sin(i*pi/(n+1)) and 2*del*(1-cos(i*pi/(n+1)) are
c            computed for i=1,2,..n and stored in trig in the
c            following way.
c
c            odd  sin   i=1,3,5,...
c            odd  cos   i=1,3,5,...
c            even sin   i=2,4,6,...
c            even cos   i=2,4,6,...
c
c     local.
c
c     fortran:        atan,sin
c
      integer n2,n4,i
      real    pi,ar,del2,del4
c
      pi  = 4.0e0*atan(1.0e0)
      ar  = pi/(n+1.0e0)
      del2= 2.0e0*del
      del4= 2.0e0*del2
      n2  = n/2
      n4  = (n2+1)/2
      do 10 i=1,n2
         w(i)=sin(i*ar)
   10    continue
      ar  = .50e0*ar
      do 20 i=1,n4
         w(n2+i)=del4*sin((2*i-1)*ar)**2
   20    continue
      trig(n4+1)    = 1.0e0
      trig(n2+n4+2) = del2
      do 30 i=1,n4
         trig(i)      = w(2*i-1)
         trig(n2+2-i) = w(2*i-1)
         trig(n2+1+i) = w(n2+i)
         trig(n+2-i)  = del4-w(n2+i)
   30    continue
      trig(n+n4+1)   = 1.0e0
      trig(n+n2+n4+1)= del2
      n4  = n2/2
      if(n4.eq.0) return
      do 40 i=1,n4
         trig(n+1+i)    = w(2*i)
         trig(n+n2+2-i) = w(2*i)
         trig(n+n2+1+i) = del4*w(i)**2
         trig(2*n+1-i)  = del4-trig(n+n2+1+i)
   40    continue
      return
      end
