      subroutine scosqf (n,x,wsave)
      integer n
      real x(*), wsave(*)
      real sqrt2, tsqx
      data sqrt2 /  1.414213562 3730950488 0168872420 970 e0 /
c
      if (n-2) 102,101,103
  101 tsqx = sqrt2*x(2)
      x(2) = x(1)-tsqx
      x(1) = x(1)+tsqx
  102 return
c
  103 call scsqf1 (n,x,wsave,wsave(n+1))
c
      return
      end
