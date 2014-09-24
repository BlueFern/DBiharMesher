      subroutine scosti (n,wsave)
      integer n
      real wsave(*)
      real dt, fk, pi
      data pi /  3.141592653 5897932384 6264338327 950 e0 /
c
      if (n .le. 3) return
c
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      dt = pi/float(nm1)
      fk = 0.e0
      do 101 k=2,ns2
         kc = np1-k
         fk = fk+1.e0
         wsave(k) = 2.e0*sin(fk*dt)
         wsave(kc) = 2.e0*cos(fk*dt)
  101 continue
c
      call srffti (nm1,wsave(n+1))
c
      return
      end
