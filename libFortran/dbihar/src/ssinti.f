      subroutine ssinti (n,wsave)
      integer n
      real wsave(*)
      real dt, fk, pi
      data pi /  3.141592653 5897932384 6264338327 950 e0 /
c
      if (n .le. 1) return
      np1 = n+1
      ns2 = n/2
      dt = pi/float(np1)
      fk = 0.e0
      do 101 k=1,ns2
         fk = fk+1.e0
         wsave(k) = 2.e0*sin(fk*dt)
  101 continue
c
      call srffti (np1,wsave(ns2+1))
c
      return
      end
