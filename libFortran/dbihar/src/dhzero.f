      subroutine dhzero(nn,kj,x,y,d,ws)
c
      integer nn,kj
      double precision    x(*),y(*),d(*),ws(*)
c
c     this routine preconditions the cg - iteration
c     by multipying with an approximation to the
c     inverse of the capacitance matrix.
c
c     x is input vector. x is not changed.
c     y is output vector. ie y=h0*x .
c     d is a diagonal matrix of dimension nn needed to define h0.
c     kj and ws are dummy parameters in this version of dhzero.
c
c     local.
c
      integer i
c
      do 10 i=1,nn
         y(i)=d(i)*x(i)
   10    continue
      return
      end
