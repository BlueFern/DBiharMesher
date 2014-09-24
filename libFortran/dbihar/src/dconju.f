      subroutine dconju(mm,nn,ki,kj,maxi,iflag,itcg,del,tol,alpha,
     +                  beta,z,s,p,y,trig,ws,diag,w3)
c
      integer mm,nn,ki,kj,maxi,iflag,itcg
      double precision    del,tol,alpha,beta
      double precision    z(*),s(*),p(*),y(*),trig(*),diag(*)
      double precision    ws(*),w3(*)
c
c     conjugate gradient routine with preconditioning.
c     z,s,p,y and diag must have n elements.
c     (z,y and diag is input,y assumed zero.)
c
c     local.
c
c     biharmonic:         dpreco, dcmult, dupdat
c     blas:               idamax, ddot,   daxpy, dscal
c     fortran             abs,sqrt
c
      integer i1,i2,kr,it
      integer km(2,2)
      integer idamax
      double precision    alf,aa,bet,bb,xx
      double precision    ddot
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
      it = idamax(nn,z,1)
      if(abs(z(it)).lt.tol*tol) return
      itcg = 1
c
      call dpreco(nn,kj,iflag,maxi,km(ki,kr),
     +            p,z,diag,w3(i1),w3(i2),ws)
      alf  = ddot(nn,z,1,p,1)
      call dcmult(mm,nn,ki,kj,del,alpha,beta,p,s,trig,ws)
      aa   = alf/ddot(nn,p,1,s,1)
      call daxpy(nn,aa,p,1,y,1)
      call dupdat(nn,kj,maxi,iflag,km(ki,kr),tol,
     +            s,p,diag,w3(i1),w3(i2),ws)
c
      do 600 it=1,30
         call daxpy(nn,-aa,s,1,z,1)
         xx   = sqrt(ddot(nn,z,1,z,1))
         if(xx.lt.tol) return
         itcg = it+1
         bet  = alf
         call dpreco(nn,kj,iflag,maxi,km(ki,kr)-1,
     +               s,z,diag,w3(i1),w3(i2),ws)
         alf  = ddot(nn,z,1,s,1)
         bb   = alf/bet
         call dscal(nn,bb,p,1)
         call daxpy(nn,1.0d0,s,1,p,1)
         call dcmult(mm,nn,ki,kj,del,alpha,beta,p,s,trig,ws)
         aa   = alf/ddot(nn,p,1,s,1)
         call daxpy(nn,aa,p,1,y,1)
         call dupdat(nn,kj,maxi,iflag,km(ki,kr),tol,
     +              s,p,diag,w3(i1),w3(i2),ws)
  600    continue
      iflag=-6
      tol  = xx
      return
      end
