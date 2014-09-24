      subroutine seffti (n,wsave)
      integer n
      real wsave(*)
c
      if (n .eq. 1) return
c
      call sefft1 (n,wsave(2*n+1),wsave(3*n+1))
c
      return
      end
