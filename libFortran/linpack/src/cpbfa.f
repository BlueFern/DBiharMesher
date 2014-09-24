      subroutine cpbfa(abd,lda,n,m,info)
      integer lda,n,m,info
      complex abd(lda,*)
c
c     cpbfa factors a complex hermitian positive definite matrix
c     stored in band form.
c
c     cpbfa is usually called by cpbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     complex(lda, n)
c                the matrix to be factored.  the columns of the upper
c                triangle are stored in the columns of abd and the
c                diagonals of the upper triangle are stored in the
c                rows of abd .  see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. m + 1 .
c
c        n       integer
c                the order of the matrix  a .
c
c        m       integer
c                the number of diagonals above the main diagonal.
c                0 .le. m .lt. n .
c
c     on return
c
c        abd     an upper triangular matrix  r , stored in band
c                form, so that  a = ctrans(r)*r .
c
c        info    integer
c                = 0  for normal return.
c                = k  if the leading minor of order  k  is not
c                     positive definite.
c
c     band storage
c
c           if  a  is a hermitian positive definite band matrix,
c           the following program segment will set up the input.
c
c                   m = (band width above diagonal)
c                   do 20 j = 1, n
c                      i1 = max0(1, j-m)
c                      do 10 i = i1, j
c                         k = i-j+m+1
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas cdotc
c     fortran aimag,cmplx,conjg,max0,real,sqrt
c
c     internal variables
c
      complex cdotc,t
      real s
      integer ik,j,jk,k,mu
c     begin block with ...exits to 40
c
c
         do 30 j = 1, n
            info = j
            s = 0.0e0
            ik = m + 1
            jk = max0(j-m,1)
            mu = max0(m+2-j,1)
            if (m .lt. mu) go to 20
            do 10 k = mu, m
               t = abd(k,j) - cdotc(k-mu,abd(ik,jk),1,abd(mu,j),1)
               t = t/abd(m+1,jk)
               abd(k,j) = t
               s = s + real(t*conjg(t))
               ik = ik - 1
               jk = jk + 1
   10       continue
   20       continue
            s = real(abd(m+1,j)) - s
c     ......exit
            if (s .le. 0.0e0 .or. aimag(abd(m+1,j)) .ne. 0.0e0)
     *         go to 40
            abd(m+1,j) = cmplx(sqrt(s),0.0e0)
   30    continue
         info = 0
   40 continue
      return
      end
