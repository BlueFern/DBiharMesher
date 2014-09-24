      subroutine sefftb (n,r,azero,a,b,wsave)
      integer n
      real r(*), azero, a(*), b(*), wsave(*)
c
      if (n-2) 101,102,103
  101 r(1) = azero
      return
c
  102 r(1) = azero+a(1)
      r(2) = azero-a(1)
      return
c
c     to supress repeated initialization, remove the following statement
c     ( call seffti(n,wsave) ) from both defftf and defftb and insert it
c     at the beginning of your program following the definition of n.
c
  103 call seffti (n,wsave)
c
      ns2 = (n-1)/2
      do 104 i=1,ns2
         r(2*i) = .5e0*a(i)
         r(2*i+1) = -.5e0*b(i)
  104 continue
      r(1) = azero
      if (mod(n,2) .eq. 0) r(n) = a(ns2+1)
c
      call srfftb (n,r,wsave(n+1))
c
      return
      end
