      subroutine sbislf(m,n,maxi,iflag,del,tol,alpha,beta,itcg,idf,f,
     +                  p,s,y,z,ws,trig,w3,diag)
c
      integer m,n,maxi,iflag,idf,itcg
      real    del,tol,alpha,beta
      real    f(idf,*)
      real    p(*),s(*),y(*),z(*),diag(*)
      real    trig(*),ws(*),w3(*)
c
c     conjugate gradient based solver of generalized
c     biharmonic equation.
c
c     p,s,y,z      must have at least (n/2+1) elements each.
c     ws           must have at least n+1 elements.
c     diag         must have at least 2*n elements.
c     w3           must have at least 2*max*(n+3) elements.
c
c     local.
c
c     biharmonic:         strigi, shzeri, spentf, sconju
c     blas:               scopy,  sscal,  saxpy
c
      integer i,ki,kj,i1,i2,itc
      integer m2,n2,ip
      real    scal1,scal2,x1,zero(1)
c
      zero(1)= 0.0e0
      itcg   = 0
      if(iflag.eq.6) go to 50
      x1     = 2.0e0/(n+1.0e0)
      scal1  = x1*(del/(m+1.0e0))**2
      scal2  = x1*.1250e0/(m+1.0e0)
      call strigi(m,del,trig,p)
      if(m.ne.n.or.del.ne.1.0e0) go to 30
      call scopy(2*n,trig,1,trig(2*m+1),1)
      go to 40
   30 call strigi(n,1.0e0,trig(2*m+1),p)
   40 call shzeri(m,n,1,del,alpha,beta,diag,trig,p)
c
   50 ip=1
      do 800 kj=1,2
        n2 = n/2+2-kj
        i2 = 2*m+1+(n+1)*(kj-1)
        do 700 ki=1,2
          m2 = m/2+2-ki
          i1 = (m+1)*(ki-1)
          call scopy(n2,zero,0,z,1)
          call scopy(n2,zero,0,y,1)
          do 400 i = 1,m2
            call scopy(n2,f(2*i+ki-2,kj),2*idf,s,1)
            x1 = scal1*trig(i1+i)
            call spentf(n2,kj,trig(i1+m2+i),alpha,beta,trig(i2),s,s,ws)
            call saxpy(n2,x1,s,1,z,1)
            call sscal(n2,scal2,s,1)
            call scopy(n2,s,1,f(2*i+ki-2,kj),2*idf)
  400       continue
c
c     the capacitance matrix equation is solved
c     here using preconditioned conjugate gradients.
c
          call sconju(m2,n2,ki,kj,maxi,iflag,itc,del,tol,alpha,beta,
     +            z,s,p,y,trig,ws,diag(ip),w3)
c
c
          itcg = itcg+itc
          do 600 i = 1,m2
            call spentf(n2,kj,trig(i1+m2+i),alpha,beta,trig(i2),y,s,ws)
            call saxpy(n2,-trig(i1+i),s,1,f(2*i+ki-2,kj),2*idf)
  600       continue
            ip=ip+n2
  700     continue
  800   continue
      itcg = itcg/4
      return
      end
