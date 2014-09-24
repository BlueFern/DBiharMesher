      subroutine sefftf (n,r,azero,a,b,wsave)
c
c                       version 3  june 1979
c
      integer n
      real r(*), azero, a(*), b(*), wsave(*)
      real cf, cfm
c
      if (n-2) 101,102,103
  101 azero = r(1)
      return
c
  102 azero = .5e0*(r(1)+r(2))
      a(1) = .5e0*(r(1)-r(2))
      return
c
c     to supress repeated initialization, remove the following statement
c     ( call seffti(n,wsave) ) from both defftf and defftb and insert it
c     at the beginning of your program following the definition of n.
c
  103 call seffti (n,wsave)
c
      do 104 i=1,n
         wsave(i) = r(i)
  104 continue
c
      call srfftf (n,wsave,wsave(n+1))
c
      cf = 2.e0/float(n)
      cfm = -cf
      azero = .5e0*cf*wsave(1)
      ns2 = (n+1)/2
      ns2m = ns2-1
      do 105 i=1,ns2m
         a(i) = cf*wsave(2*i)
         b(i) = cfm*wsave(2*i+1)
  105 continue
      if (mod(n,2) .eq. 0) a(ns2) = .5e0*cf*wsave(n)
c
      return
      end
