      subroutine srfftf (n,r,wsave)
      integer n
      real r(*), wsave(*)
c
      if (n .eq. 1) return
c
      call srftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
c
      return
      end
