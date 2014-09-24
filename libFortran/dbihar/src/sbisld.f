      subroutine sbisld(m,n,iflag,del,alpha,beta,idf,f,
     +                  s,z,ws,trig,cmat)
c
      integer m,n,iflag,idf
      real del,alpha,beta
      real f(idf,*),cmat(*)
      real s(*),z(*),ws(*),trig(*)
c
c     direct solver for generalized biharmonic problem.
c     note that a real array is passed to linpack to store
c     pivot information if iflag is 4 or 8.
c     s,z and ws must have n2/2+1  elements each.
c     trig       must have 2*(m+n) elements.
c     cmat       must have (n+1)*(n+1)/2   if iflag is 3 or 7.
c     cmat       must have (n+1)*(n+1)/2+2*n if iflag is 4 or 8.
c
c     local.
c
c     biharmonic:         strigi, spentf, smatge
c     linpack:            sspfa,  sspsl,  sppfa, sppsl
c     blas:               scopy,  saxpy,  sscal
c
      integer i,ki,kj,i1,i2
      integer m2,n2,ip,iq,info
      real    scal1,scal2,x1,zero(1)
c
      zero(1) = 0.0e0
      if(iflag.eq.7.or.iflag.eq.8) go to 50
      x1     = 2.0e0/(n+1.0e0)
      scal1  = x1*(del/(m+1.0e0))**2
      scal2  = x1*.1250e0/(m+1.0e0)
      call strigi(m,del,trig,s)
      if(m.ne.n.or.del.ne.1.0e0) go to 40
      call scopy(2*n,trig,1,trig(2*m+1),1)
      go to 50
   40 call strigi(n,1.0e0,trig(2*m+1),s)
   50 ip     = 1
      iq     = 0
      do 800 kj=1,2
        n2 = n/2+2-kj
        i2 = 2*m+1+(n+1)*(kj-1)
        if(iflag.eq.4.or.iflag.eq.8) iq=n2
        do 700 ki=1,2
          i1 = (m+1)*(ki-1)
          m2 = m/2+2-ki
          call scopy(n2,zero,0,z,1)
          do 400 i = 1,m2
            call scopy(n2,f(2*i+ki-2,kj),2*idf,s,1)
            x1 = scal1*trig(i1+i)
            call spentf(n2,kj,trig(i1+m2+i),alpha,beta,trig(i2),s,s,ws)
            call saxpy(n2,x1,s,1,z,1)
            call sscal(n2,scal2,s,1)
            call scopy(n2,s,1,f(2*i+ki-2,kj),2*idf)
  400       continue
c
c     the capacitance matrix equation is solved using linpack.
c
          if(iflag.eq.7) go to 450
          if(iflag.eq.8) go to 440
          call smatge(m2,n2,ki,kj,del,alpha,beta,trig,cmat(ip+iq),ws)
          if(iflag.eq.3) go to 445
          call sspfa(cmat(ip+iq),n2,cmat(ip),info)
          if(info.ne.0) go to 1020
  440     call sspsl(cmat(ip+iq),n2,cmat(ip),z)
          go to 460
  445     call sppfa(cmat(ip),n2,info)
          if(info.ne.0) go to 1010
  450     call sppsl(cmat(ip),n2,z)
c
  460     do 600 i = 1,m2
            call spentf(n2,kj,trig(i1+m2+i),
     +                 alpha,beta,trig(i2),z,s,ws)
            call saxpy(n2,-trig(i1+i),s,1,f(2*i+ki-2,kj),2*idf)
  600       continue
          ip = ip+iq+(n2*(n2+1))/2
  700     continue
  800   continue
      return
 1010 iflag  = -4
      return
 1020 iflag  = -5
      return
      end
