c /***************************************************
c   hgext.f  for dbulhalo
c  21 Mar. 2015  written by D.Kawata
c ****************************************************/
c70*******************************************************************


c /*****   add hot gas density profile to mext.dat   *****/
c also set mrhg, rmhg, rhorhg
      subroutine hgext(mtothg,cnfw,rvir,xcut,rt)
      include 'define.f'

      double precision mtothg
      double precision omg0,omz,hub,cnfw,rvir,zt
      double precision fb
      double precision rhoc,delc
      double precision gcv,rho0hg
      double precision rs
c for work
      integer i,flag
      character filen*60
c *** for radial profile ***
      double precision ri,ro,xcut,lri,lro,dlr,dr,rt
      double precision rhoi,rhoo,xo,mro
c *** external function 
      double precision rhombhg,mrmbhg
      external rhombhg,mrmbhg

c hot gas mass profile parameter from thermal core NFW from Maller & Bullock (2004)

      rs=rvir/cnfw
      gcv=9.0d0*dlog(1.0d0+(4.0d0/3.0d0)*cnfw)
     &   -8.0d0*dlog(1.0d0+cnfw)-(4.0d0*cnfw/(1.0d0+cnfw))
      rho0hg=mtothg/(4.0d0*M_PI*(rs**3)*gcv)

      if(myrank.eq.0) then
        open(60,file='mbhgprof.dat',status='unknown')
      endif

      ri = frini0*rt
      ro = xcut*rt
      lri = dlog10(ri)
      lro = dlog10(ro)
      dlr = (lro-lri)/dble(ndr)
      ro = 0.0d0
      do i=0,ndr
        ri=ro
        ro=10.0d0**(lri+dlr*dble(i))
        dr=ro-ri
c Maller & Bullock profile
        xo=ro/rs
        rhoo=rho0hg*rhombhg(xo)
        mro=mtothg*mrmbhg(xo)/gcv
        if(myrank.eq.0) then
          write(60,'(3(1pE13.5))') ro*LUKPC,rhoo,mro
        endif
        mrdmh(i)=mrdmh(i)+mro
        rdmh(i) = ro
        rhodmh(i)=rhodmh(i)+rhoo
c set mass profile
        mrhg(i)=mro
        rmhg(i)=ro
        rhorhg(i)=rhoo
      enddo
      if(myrank.eq.0) then
        close(60)
      endif

c *** set common values
      lrodmh=lro
      lridmh=lri
      dlrdmh=dlr

      if(myrank.eq.0) then
        open(63,file='mexta_hg.dat',status='unknown')
        do i=0,ndr
          write(63,'(3(1pE13.5))') rdmh(i)*LUKPC,mrdmh(i)
     &     ,rhodmh(i)
        enddo
        close(63)
      endif

      return
      end
