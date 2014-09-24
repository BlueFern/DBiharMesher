      subroutine sconju(mm,nn,ki,kj,maxi,iflag,itcg,del,tol,alpha,
     +                  beta,z,s,p,y,trig,ws,diag,w3)
c
      integer mm,nn,ki,kj,maxi,iflag,itcg
      real    del,tol,alpha,beta
      real    z(*),s(*),p(*),y(*),trig(*),diag(*)
      real    ws(*),w3(*)
c
c     conjugate gradient routine with preconditioning.
c     z,s,p,y and diag must have n elements.
c     (z,y and diag is input,y assumed zero.)
c
c     local.
c
c     biharmonic:         spreco, scmult, supdat
c     blas:               isamax, sdot,   saxpy, sscal
c     fortran             abs,sqrt
c
      integer i1,i2,kr,it
      integer km(2,2)
      integer isamax
      real    alf,aa,bet,bb,xx
      real    sdot
c
      itcg = 0
      if(kj.eq.0) go to 5
      kr = kj
      i1 = maxi*(ki+2*kj-3)+1
      i2 = maxi*(ki+2*kj-3)*(nn+kj-1)+4*maxi+1
      go to 10
    5 kr = 1
      i1 = (ki-1)*maxi+1
      i2 = (ki-1)*maxi*nn+2*maxi+1
   10 if(iflag.le.2) km(ki,kr) = 0
      it = isamax(nn,z,1)
      if(abs(z(it)).lt.tol*tol) return
      itcg = 1
c
      call spreco(nn,kj,iflag,maxi,km(ki,kr),
     +            p,z,diag,w3(i1),w3(i2),ws)
      alf  = sdot(nn,z,1,p,1)
      call scmult(mm,nn,ki,kj,del,alpha,beta,p,s,trig,ws)
      aa   = alf/sdot(nn,p,1,s,1)
      call saxpy(nn,aa,p,1,y,1)
      call supdat(nn,kj,maxi,iflag,km(ki,kr),tol,
     +            s,p,diag,w3(i1),w3(i2),ws)
c
      do 600 it=1,30
         call saxpy(nn,-aa,s,1,z,1)
         xx   = sqrt(sdot(nn,z,1,z,1))
         if(xx.lt.tol) return
         itcg = it+1
         bet  = alf
         call spreco(nn,kj,iflag,maxi,km(ki,kr)-1,
     +               s,z,diag,w3(i1),w3(i2),ws)
         alf  = sdot(nn,z,1,s,1)
         bb   = alf/bet
         call sscal(nn,bb,p,1)
         call saxpy(nn,1.0e0,s,1,p,1)
         call scmult(mm,nn,ki,kj,del,alpha,beta,p,s,trig,ws)
         aa   = alf/sdot(nn,p,1,s,1)
         call saxpy(nn,aa,p,1,y,1)
         call supdat(nn,kj,maxi,iflag,km(ki,kr),tol,
     +              s,p,diag,w3(i1),w3(i2),ws)
  600    continue
      iflag=-6
      tol  = xx
      return
      end
