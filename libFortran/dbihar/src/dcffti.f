      subroutine dcffti (n,wsave)
      integer n
      double precision wsave(*)
c
      if (n .eq. 1) return
c
      iw1 = n+n+1
      iw2 = iw1+n+n
      call dcfti1 (n,wsave(iw1),wsave(iw2))
c
      return
      end
