      subroutine ssinqf (n,x,wsave)
      integer n
      real x(*), wsave(*)
      real xhold
c
      if (n .eq. 1) return
c
      ns2 = n/2
      do 101 k=1,ns2
         kc = n-k
         xhold = x(k)
         x(k) = x(kc+1)
         x(kc+1) = xhold
  101 continue
c
      call scosqf (n,x,wsave)
c
      do 102 k=2,n,2
         x(k) = -x(k)
  102 continue
c
      return
      end
