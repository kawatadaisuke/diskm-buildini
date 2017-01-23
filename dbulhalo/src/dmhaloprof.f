c /***************************************************
c   dmhalogen.f  for dbulhalo
c  18 Jan. 2011  written by D.Kawata
c ****************************************************/
c70*******************************************************************


c /*****   generating NFW DM halo or mext.dat *****/
      subroutine dmhaloprof(omg0,omz,hub,cnfw,rvir,zt
     &  ,xcut,rt,fb,flag)
      include 'define.f'

      integer flag
      double precision omg0,omz,hub,cnfw,rvir,zt
      double precision fb
      double precision rhoc,delc
c for work
      integer i
      character filen*60
c *** for radial profile ***
      double precision ri,ro,xcut,lri,lro,dlr,dr,rt,rs
      double precision mr,pmr,facnfw
      double precision rhoi,rhoo,x
c *** external function 
      integer idum
      real ran1
      external ran1
      double precision rhonfwsw,mrnfw
      external rhonfwsw,mrnfw

      idum = -111

      rhoc = (3.0d0*HUB0*HUB0*hub*hub/(8.0d0*M_PI*G))
     &  *(omg0/omz)*((1.0d0+zt)**3)
      delc = (200.0d0/3.0d0)*(cnfw**3/(dlog(1.0d0+cnfw)
     &  -cnfw/(1.0d0+cnfw)))  
      facnfw = delc*rhoc
      if(myrank.eq.0) then
        write(6,*) ' rhoc,delc,facnfw=',rhoc,delc,facnfw
      endif

c /*** integral of NFW profile ****/
      ri = frini0*rt
      ro = xcut*rt
      lri = dlog10(ri)
      lro = dlog10(ro)
      dlr = (lro-lri)/dble(ndr)
      ro = 0.0d0
      mr = 0.0d0
      do i=0,ndr
        if(i.ge.0) then
          rhoi = rhoo
          pmr = mr
        endif
        ri = ro
        ro = 10.0d0**(lri+dlr*dble(i))
        dr = ro-ri
        rhoo = rhonfwsw(ro,cnfw,rvir,rt,flag)
        if(i.lt.0) then
          rhoi = rhoo
        endif
        mr = mr+4.0d0*M_PI*0.5d0*(rhoo*(ro**2)+rhoi*(ri**2))
     &    *dr
        if(i.lt.0) then
          pmr = mr
        endif
        mrdmh(i) = (1.0d0-fb)*mr*facnfw
        rdmh(i) = ro
        rhodmh(i)=(1.0d0-fb)*rhoo*facnfw
      enddo

c *** set common values
      lrodmh=lro
      lridmh=lri
      dlrdmh=dlr
c NFW rs
      rs=rvir/cnfw

      if(myrank.eq.0) then
        open(62,file='mext-dmh.dat',status='unknown'
     &   ,form='unformatted')
        open(63,file='mexta-dmh.dat',status='unknown')
        write(62) ndr+1
        write(62) dlr,lri,lro
        do i=0,ndr
          x=rdmh(i)
          write(63,'(6(1pE13.5))') rdmh(i)*LUKPC,mrdmh(i)
     &     ,rhodmh(i)
     &     ,(1.0d0-fb)*facnfw*mrnfw(x,cnfw,rvir)
c escape velocity BT2nd p.71
     &     ,VUKMS*dsqrt(2.0d0*(1.0d0-fb)*facnfw*4.0d0*M_PI*(rs**2)
     &     *dlog(1.0d0+x/rs)/(x/rs)),x
          write(62) rdmh(i),mrdmh(i)
        enddo
        close(62)
        close(63)
      endif

      return
      end
