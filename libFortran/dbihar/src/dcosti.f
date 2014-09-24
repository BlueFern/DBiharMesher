      subroutine dcosti (n,wsave)
      integer n
      double precision wsave(*)
      double precision dt, fk, pi
      data pi /  3.141592653 5897932384 6264338327 950 d0 /
c
      if (n .le. 3) return
c
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      dt = pi/dfloat(nm1)
      fk = 0.d0
      do 101 k=2,ns2
         kc = np1-k
         fk = fk+1.d0
         wsave(k) = 2.d0*dsin(fk*dt)
         wsave(kc) = 2.d0*dcos(fk*dt)
  101 continue
c
      call drffti (nm1,wsave(n+1))
c
      return
      end
