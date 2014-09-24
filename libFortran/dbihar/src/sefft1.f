      subroutine sefft1 (n,wa,ifac)
      integer n
      real wa(*)
      integer ifac(*)
      real arg1, argh, ch1, ch1h, dch1, dsh1, sh1, tpi
      integer ntryh(4)
      data ntryh(1), ntryh(2), ntryh(3), ntryh(4) /4, 2, 3, 5/
      data tpi   /  6.2831853071 7958647692 5286766559 00577e0/
c
      nl = n
      nf = 0
      j = 0
c
  101 j = j+1
      if (j.le.4) ntry = ntryh(j)
      if (j.gt.4) ntry = ntry + 2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr.ne.0) go to 101
c
  105 nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
  106 continue
      ifac(3) = 2
c
  107 if (nl .ne. 1) go to 104
c
      ifac(1) = n
      ifac(2) = nf
      argh = tpi/float(n)
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
c
      do 111 k1=1,nfm1
         ip = ifac(k1+2)
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         arg1 = float(l1)*argh
         ch1 = 1.e0
         sh1 = 0.e0
         dch1 = cos(arg1)
         dsh1 = sin(arg1)
c
         do 110 j=1,ipm
            ch1h = dch1*ch1-dsh1*sh1
            sh1 = dch1*sh1+dsh1*ch1
            ch1 = ch1h
            i = is+2
            wa(i-1) = ch1
            wa(i) = sh1
            if (ido .lt. 5) go to 109
            do 108 ii=5,ido,2
               i = i+2
               wa(i-1) = ch1*wa(i-3)-sh1*wa(i-2)
               wa(i) = ch1*wa(i-2)+sh1*wa(i-3)
  108       continue
  109       is = is+ido
  110    continue
c
         l1 = l2
  111 continue
c
      return
      end
