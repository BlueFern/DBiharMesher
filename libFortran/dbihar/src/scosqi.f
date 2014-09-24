      subroutine scosqi (n,wsave)
      integer n
      real wsave(*)
      real dt, fk, pih
      data pih /  1.570796326 7948966192 3132169163 975 e0 /
c
      dt = pih/float(n)
      fk = 0.e0
      do 101 k=1,n
         fk = fk+1.e0
         wsave(k) = cos(fk*dt)
  101 continue
c
      call srffti (n,wsave(n+1))
c
      return
      end
