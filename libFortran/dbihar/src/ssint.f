      subroutine ssint (n,x,wsave)
      integer n
      real x(*), wsave(*)
      real sqrt3, t1, t2, x1, xh, xim1
      data sqrt3 /  1.7320508075 6887729352 7446341505 87237e0/
c
      if (n-2) 101,102,103
  101 x(1) = x(1)+x(1)
      return
c
  102 xh = sqrt3*(x(1)+x(2))
      x(2) = sqrt3*(x(1)-x(2))
      x(1) = xh
      return
c
  103 np1 = n+1
      ns2 = n/2
      x1 = x(1)
      x(1) = 0.e0
      do 104 k=1,ns2
         kc = np1-k
         t1 = x1-x(kc)
         t2 = wsave(k)*(x1+x(kc))
         x1 = x(k+1)
         x(k+1) = t1+t2
         x(kc+1) = t2-t1
  104 continue
      modn = mod(n,2)
      if (modn .ne. 0) x(ns2+2) = 4.e0*x1
c
      call srfftf (np1,x,wsave(ns2+1))
c
      x(1) = .5e0*x(1)
      do 105 i=3,n,2
         xim1 = x(i-1)
         x(i-1) = -x(i)
         x(i) = x(i-2)+xim1
  105 continue
      if (modn.eq.0) x(n) = -x(n+1)
c
      return
      end
