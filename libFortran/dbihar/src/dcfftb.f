      subroutine dcfftb (n,c,wsave)
      integer n
      double precision c(*), wsave(*)
c
      if (n .eq. 1) return
c
      iw1 = n+n+1
      iw2 = iw1+n+n
      call dcftb1 (n,c,wsave,wsave(iw1),wsave(iw2))
c
      return
      end
