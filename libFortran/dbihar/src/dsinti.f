      subroutine dsinti (n,wsave)
      integer n
      double precision wsave(*)
      double precision dt, fk, pi
      data pi /  3.141592653 5897932384 6264338327 950 d0 /
c
      if (n .le. 1) return
      np1 = n+1
      ns2 = n/2
      dt = pi/dfloat(np1)
      fk = 0.d0
      do 101 k=1,ns2
         fk = fk+1.d0
         wsave(k) = 2.d0*dsin(fk*dt)
  101 continue
c
      call drffti (np1,wsave(ns2+1))
c
      return
      end
