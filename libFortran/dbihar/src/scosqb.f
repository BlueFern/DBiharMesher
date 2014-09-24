      subroutine scosqb (n,x,wsave)
      integer n
      real x(*), wsave(*)
      real tsqrt2, x1
      data tsqrt2 /  2.828427124 7461900976 0337744841 94 e0 /
c
      if (n-2) 101,102,103
  101 x(1) = 4.e0*x(1)
      return
c
  102 x1 = 4.e0*(x(1)+x(2))
      x(2) = tsqrt2*(x(1)-x(2))
      x(1) = x1
      return
c
  103 call scsqb1 (n,x,wsave,wsave(n+1))
c
      return
      end
