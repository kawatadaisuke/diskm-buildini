c /***************************************************
c   hgtemprof.f
c  20 Apr. 2015  written by D.Kawata
c ****************************************************/
c70*******************************************************************


c ***  generating hot gas with a constant entropy profile ***
      subroutine hgtemprof(mtothg,cnfw,rvir,rcuthg,mtothgrc)
      include 'define.f'

      double precision mtothg,cnfw,rvir,rcuthg
      double precision mtothgrc
c for work
      integer i,ndrh,ir
      character filen*60
c *** for radial profile ***
      double precision gcv,rho0hg,rs
      double precision tmr,pdtmr,dtmr,dr
      double precision ltmr,ltmri
      double precision mrtot
c *** for integration 
      double precision lpri,rhoro,rhor,ndenr,pr,ur,gkmm,temr,fr,lrhor
      double precision ri,rr,ro,xr,rp
      double precision dlr,lri,lro
      double precision dltmr,dltmr1,dltmr2,dltmr3,dltmr4
      double precision frho0
      double precision fb,facnfw
c for extention of gas density profile
      double precision dlrhodlr,lr0,lrho0,lrr
c *** external function 
      double precision rhombhg,mrmbhg
      external rhombhg,mrmbhg


      gkmm=(GAM-1.0d0)/((KCGS/(MP*MYU))/K_MU)
      rs=rvir/cnfw
      gcv=9.0d0*dlog(1.0d0+(4.0d0/3.0d0)*cnfw)
     &   -8.0d0*dlog(1.0d0+cnfw)-(4.0d0*cnfw/(1.0d0+cnfw))
      rho0hg=mtothg/(4.0d0*M_PI*(rs**3)*gcv)
      ndrh=ndr

c *** integral hot gas profile ****
c *** reset the radial grid for hot gas profile ***
      ri=frini0*rcuthg
      ro=rcuthg
      lri=dlog(ri)
      lro=dlog(ro)
      dlr=(lro-lri)/dble(ndr)
      do i=0,ndr
        rmhg(i)=dexp(lri+dble(i)*dlr)
        rr=rmhg(i)
        xr=rr/rs
        rhorhg(i)=rho0hg*rhombhg(xr)
        mrhg(i)=mtothg*mrmbhg(xr)/gcv
      enddo
! total mass with rcuthg
      mtothgrc=mrhg(ndr)
      if(myrank.eq.0) then
        write(6,*)
        write(6,*) 'Total hot gas mass <rcuthg =',mtothgrc
        write(6,*) 'rcuthg =',rcuthg
        write(6,*)
      endif

c     *** getting temperature profile
c first step from the outer boundary
      temr=0.0d0
      rr=rmhg(ndr)
c total mass, mrtot set
      rp=rr
      ir=int((dlog10(rp)-lridmh)/dlrdmh)
      if(ir.lt.0) then
        ir=0
      else if(ir.ge.ndr) then
        ir=ndr-1
      endif
      mrtot=mrdmh(ir)+(rp-rdmh(ir))
     &  *(mrdmh(ir+1)-mrdmh(ir))/(rdmh(ir+1)-rdmh(ir))
c 
      pdtmr=rhorhg(ndr)*G*mrtot/(rr**2)

      i=ndr-1
      rr=rmhg(i)
c total mass, mrtot set
      rp=rr
      ir=int((dlog10(rp)-lridmh)/dlrdmh)
      if(ir.lt.0) then
        ir=0
      else if(ir.ge.ndr) then
        ir=ndr-1
      endif
      mrtot=mrdmh(ir)+(rp-rdmh(ir))
     &  *(mrdmh(ir+1)-mrdmh(ir))/(rdmh(ir+1)-rdmh(ir))
c 
      dtmr=rhorhg(i)*G*mrtot/(rr**2)
      dr=rmhg(i+1)-rmhg(i)
c first step with trapezium law
      temr=temr+dr*0.5d0*(dtmr+pdtmr)
      temhg(i)=(1.0d0/((KCGS/(MP*MYU))/K_MU))
     &  *rhorhg(i)*temr
      temhg(ndr)=temhg(i)
c log integration  
      ltmr=dlog(temr)
c RK integration
      do i=ndr-2,0,-1
        ltmri=ltmr
c *** step 1
        rr=rmhg(i+1)
c total mass, mrtot set
        rp=rr
        ir=int((dlog10(rp)-lridmh)/dlrdmh)
        if(ir.lt.0) then
          ir=0
        else if(ir.ge.ndr) then
          ir=ndr-1
        endif
        mrtot=mrdmh(ir)+(rp-rdmh(ir))
     &    *(mrdmh(ir+1)-mrdmh(ir))/(rdmh(ir+1)-rdmh(ir))
c 
        tmr=mrtot
        rhor=rhorhg(i+1)
        dltmr1=dlr*(rhor*G*tmr/rr)/temr
c *** at r i+0.5
        rr=dexp(lri+(dble(i)+0.5)*dlr)
c total mass, mrtot set
        rp=rr
        ir=int((dlog10(rp)-lridmh)/dlrdmh)
        if(ir.lt.0) then
          ir=0
        else if(ir.ge.ndr) then
          ir=ndr-1
        endif
        mrtot=mrdmh(ir)+(rp-rdmh(ir))
     &    *(mrdmh(ir+1)-mrdmh(ir))/(rdmh(ir+1)-rdmh(ir))
c 
        tmr=mrtot
        xr=rr/rs
        rhor=rho0hg*rhombhg(xr)
c *** step 2
        ltmr=ltmri+0.5d0*dltmr1
        temr=dexp(ltmr)
        dltmr2=dlr*(rhor*G*tmr/rr)/temr
c *** step 3
        ltmr=ltmri+0.5d0*dltmr2
        temr=dexp(ltmr)
        dltmr3=dlr*(rhor*G*tmr/rr)/temr
c *** at ri
        rr=dexp(lri+dble(i)*dlr)
c total mass, mrtot set
        rp=rr
        ir=int((dlog10(rp)-lridmh)/dlrdmh)
        if(ir.lt.0) then
          ir=0
        else if(ir.ge.ndr) then
          ir=ndr-1
        endif
        mrtot=mrdmh(ir)+(rp-rdmh(ir))
     &    *(mrdmh(ir+1)-mrdmh(ir))/(rdmh(ir+1)-rdmh(ir))
c 
        tmr=mrtot
        rhor=rhorhg(i)
c *** step 4
        ltmr=ltmri+dltmr3
        temr=dexp(ltmr)
        dltmr4=dlr*(rhor*G*tmr/rr)/temr
c *** final
        dltmr=dltmr1/6.0d0+dltmr2/3.0d0+dltmr3/3.0d0+dltmr4/6.0d0
        ltmr=ltmri+dltmr
        temr=dexp(ltmr)
        temhg(i)=(1.0d0/((KCGS/(MP*MYU))/K_MU))
     &    *rhorhg(i)*temr
      enddo

      open(62,file='hgtemprof.dat',status='unknown')
      do i=0,ndrh 
        rr=rmhg(i)
        ndenr=rhorhg(i)*(DU/(MYU*MP))
c temhg is 1e4 K unit
        uhg(i)=temhg(i)/gkmm
c total mass, mrtot set
        rp=rr
        ir=int((dlog10(rp)-lridmh)/dlrdmh)
        if(ir.lt.0) then
          ir=0
        else if(ir.ge.ndr) then
          ir=ndr-1
        endif
        mrtot=mrdmh(ir)+(rp-rdmh(ir))
     &    *(mrdmh(ir+1)-mrdmh(ir))/(rdmh(ir+1)-rdmh(ir))
        write(62,'(7(1pE13.5))') rmhg(i)*LUKPC,mrtot
     &     ,rhorhg(i),uhg(i),temhg(i)*TUK,ndenr
     &     ,rhorhg(i)*mrtot/(rmhg(i)**2)
      enddo
      close(62)

      return
      end
