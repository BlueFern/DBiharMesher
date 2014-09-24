      subroutine dcosqi (n,wsave)
      integer n
      double precision wsave(*)
      double precision dt, fk, pih
      data pih /  1.570796326 7948966192 3132169163 975 d0 /
c
      dt = pih/dfloat(n)
      fk = 0.d0
      do 101 k=1,n
         fk = fk+1.d0
         wsave(k) = dcos(fk*dt)
  101 continue
c
      call drffti (n,wsave(n+1))
c
      return
      end
