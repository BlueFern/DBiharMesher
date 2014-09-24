      subroutine dsinqi (n,wsave)
      integer n
      double precision wsave(*)
c
      call dcosqi (n,wsave)
c
      return
      end
