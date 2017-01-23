c /***************************************************
c   expdgen.f  for dbulhalo
c  11 June. 2015  written by D.Kawata
c ****************************************************/
c70*******************************************************************

c /*****   generating exponential disk *****/
      subroutine expdgen(id,nd,mtdisk,mdisk,hd,zd0,rdlim,zdlim
     &   ,ngr,ngz,flagb
     &   ,omg0,omz,hub,cnfw,rvir,zt,fb,flagmexth
     &   ,flagflzd,Rflzd)
      include 'define.f'

      integer nd,id,flagb,flagflzd
      integer ngr,ngz
      character filen*60
      double precision mtdisk,mdisk,hd,zd0,rdlim,zdlim,Rflzd
      integer i,j
      integer nrp,ndrm
      double precision ri,ro,lri,lro,dlr,rp,mdp,rr,mrp,zdr
      double precision pfr,ph
c *** Exp disk profile ***
      double precision mrexpd(MNR),rexpd(MNR)
      double precision mzexpd(MNR),zexpd(MNR)
c *** for particles ***
      double precision xbp,ybp,zbp,vxbp,vybp,vzbp,mbp,rhobp
c *** for velocity ***
      integer ir,iz
      integer pnrizi,pnrizo,pnrozi,pnrozo
      double precision wrizi,wrizo,wrozi,wrozo
      double precision zp,dr,dz
      double precision vsigz,vsigr,vsigph,vphm
      double precision vrp,vphp
c for NFW halo parameter
      integer flagmexth       
      double precision omg0,omz,hub,cnfw,rvir,zt,fb
      double precision rhoc,delc,rs
      double precision fvesc,vesc2r,v2p,r3dp

c *** external function ***
      double precision mrexpdf
      external mrexpdf
      double precision rhonfwsw,mrnfw
      external rhonfwsw,mrnfw
c *** external function 
      integer idum
      real ran1,gasdev
      external ran1,gasdev

      idum = -111+id

c NFW profile constant
      rhoc = (3.0d0*HUB0*HUB0*hub*hub/(8.0d0*M_PI*G))
     &  *(omg0/omz)*((1.0d0+zt)**3)
      delc = (200.0d0/3.0d0)*(cnfw**3/(dlog(1.0d0+cnfw)
     &  -cnfw/(1.0d0+cnfw)))  
      rs=rvir/cnfw
c escape velocity constant from BT2nd p.71
      fvesc=2.0d0*4.0d0*M_PI*G*(1.0d0-fb)*delc*rhoc*(rs**2)
      r3dp=0.1d0
      if(myrank.eq.0) then
        write(6,*) ' rhoc,delc,fb=',rhoc,delc,fb
        write(6,*) ' esc velocity at 10 kpc=',VUKMS*dsqrt(fvesc
     &    *dlog(1.0d0+r3dp/rs)/(r3dp/rs))
      endif

      ri = 1.0e-6*rdlim
      ro = rdlim
      lri = dlog10(ri)
      lro = dlog10(ro)
      dlr = (lro-lri)/dble(ndr-1)
c *** set mr and r ***
      rexpd(1) = 0.0d0
      mrexpd(1) = 0.0d0
      do i=1,ndr
        rexpd(i) = 10.0d0**(lri+dlr*dble(i-1))
        rr = rexpd(i)
c        mrexpd(i) = 1.0d0-(1.0d0+rr/hd)*dexp(-rr/hd)
        mrexpd(i)=mrexpdf(rr,hd,mtdisk)/mdisk
      enddo
c *** set mz and z ***
      ri = 1.0e-7*zdlim
      ro = zdlim
      lri = dlog10(ri)
      lro = dlog10(ro)
      dlr = (lro-lri)/dble(ndr-1)
c *** set mr and r ***
      do i=1,ndr
        zexpd(i) = 10.0d0**(lri+dlr*dble(i-1))
        rr = zexpd(i)
c Integral sech(z/zd)^2 is zd tanh(z/zd)
        mzexpd(i) = dtanh(rr/zd0)
      enddo
      mdp = mdisk/dble(nd)

      write(6,*) id,' disk particle mass=',mdp

c *** generating particles ***
      write(filen,'(a12,i1,a4)') 'output/diskc',id,'.dat'
      if(flagb.eq.0) then
        open(60,file=filen,status='unknown',form='unformatted')
        write(60) nd,0
      else
        open(60,file=filen,status='unknown')
        write(60,'(a8,2I8)') '#n,flag=',nd,0
      endif

      do i=0,nd-1
 67     pfr = dble(ran1(idum))
c /*** search radius ***/
        do j=1,ndr
          if(mrexpd(j).gt.pfr) then
            goto 72
          endif
        enddo
c        write(6,*) ' Warning: in finding radius for disk'
c        write(6,*) '  pfr =',pfr
        goto 67
 72     if(j.eq.1) then
c          write(6,*) ' Warning: in finding radius for disk: j,mr,prf='
c     &      ,j,mrexpd(j),pfr
          j = 2
        endif
        rp=10.0d0**(dlog10(rexpd(j-1))
     &    +(dlog10(rexpd(j))-dlog10(rexpd(j-1)))
     &    *(dlog10(pfr)-dlog10(mrexpd(j-1)))
     &    /(dlog10(mrexpd(j))-dlog10(mrexpd(j-1))))
        ph =  2.0d0*M_PI*dble(ran1(idum))
        mbp = mdp
        xbp = rp*dcos(ph)
        ybp = rp*dsin(ph)
c *** z-axis ****
 66     pfr = dble(ran1(idum))
c /*** search z ***/
        do j=1,ndr
          if(mzexpd(j).gt.pfr) then
            goto 73
          endif
        enddo
c        write(6,*) ' Warning: in finding z for disk'
c        write(6,*) '  pfr =',pfr
        goto 66
 73     if(j.eq.1) then
c          write(6,*) ' Warning: in finding z for disk: j=1,mz,prf'
c     &    ,mzexpd(j),pfr
          j = 2
        endif
        zbp = 10.0d0**(dlog10(zexpd(j-1))
     &    +(dlog10(zexpd(j))-dlog10(zexpd(j-1)))
     &    *(dlog10(pfr)-dlog10(mzexpd(j-1)))
     &    /(dlog10(mzexpd(j))-dlog10(mzexpd(j-1))))
        if(flagflzd.ne.0) then
c z for scale-height zd0, with zdr z2=zdr z1/zd1
          zbp=zbp*exp(rp/Rflzd)
          zdr=zd0*exp(rp/Rflzd)
        else
          zdr=zd0
        endif
c *** absolute value of z
        zp=zbp
        if(dble(ran1(idum)).lt.0.5d0) then
          zbp=-zbp
        endif
c *** density ***
        rhobp=(mtdisk/(4.0d0*M_PI*zdr*(hd**2)))
     &    *(1.0d0/(dcosh(zbp/zdr)**2))
     &    *dexp(-rp/hd)
c *** velocity ***
c *** find the point in r_rzg ***
        if(rp.gt.0.0d0) then
          ir=int((dlog10(rp)-lri_rzg)/dlr_rzg)+1
        else
          ir=0
        endif
        if(ir.lt.0) then
          ir=0
        else if(ir.gt.ngr-1) then
          ir=ngr-1
        endif
c *** find the point in z_rzg ***
        if(zp.gt.0.0d0) then
          iz=int((dlog10(zp)-lzi_rzg)/dlz_rzg)+1
        else
          iz=0
        endif
        if(iz.lt.0) then
          iz=0
        else if(iz.gt.ngz-1) then
          iz=ngz-1
        endif
c *** grid ***
        pnrizi=id_rzg(ir,iz)
        pnrizo=id_rzg(ir,iz+1)
        pnrozi=id_rzg(ir+1,iz)
        pnrozo=id_rzg(ir+1,iz+1)
        dr=x_rzg(pnrozi)-x_rzg(pnrizi)
        dz=z_rzg(pnrizo)-z_rzg(pnrizi)
c weight 
        wrizi=(x_rzg(pnrozi)-rp)*(z_rzg(pnrizo)-zp)/(dr*dz)
        wrizo=(x_rzg(pnrozo)-rp)*(zp-z_rzg(pnrizi))/(dr*dz)
        wrozi=(rp-x_rzg(pnrizi))*(z_rzg(pnrozo)-zp)/(dr*dz)
        wrozo=(rp-x_rzg(pnrizo))*(zp-z_rzg(pnrozi))/(dr*dz)
c *** vphm
        vphm=wrizi*vph_rzg(pnrizi)+wrizo*vph_rzg(pnrizo)
     &   +wrozi*vph_rzg(pnrozi)+wrozo*vph_rzg(pnrozo)
        if(vphm.lt.0.0d0) then
          vphm=0.0d0
        endif
c *** vsigph
        vsigph=wrizi*vsigph_rzg(pnrizi)+wrizo*vsigph_rzg(pnrizo)
     &   +wrozi*vsigph_rzg(pnrozi)+wrozo*vsigph_rzg(pnrozo)
        if(vsigph.lt.0.0d0) then
          vsigph=0.0d0
        endif
c *** vsigr
        vsigr=wrizi*vsigr_rzg(pnrizi)+wrizo*vsigr_rzg(pnrizo)
     &   +wrozi*vsigr_rzg(pnrozi)+wrozo*vsigr_rzg(pnrozo)
        if(vsigr.lt.0.0d0) then
          vsigr=0.0d0
        endif
c *** vsigz
        vsigz=wrizi*vsigz_rzg(pnrizi)+wrizo*vsigz_rzg(pnrizo)
     &   +wrozi*vsigz_rzg(pnrozi)+wrozo*vsigz_rzg(pnrozo)
        if(vsigz.lt.0.0d0) then
          vsigz=0.0d0
        endif
c *** get vph,vr,vz
 74     vphp=vphm+vsigph*dble(gasdev(idum))
        vrp=vsigr*dble(gasdev(idum))
        vzbp=vsigz*dble(gasdev(idum))
c *** convert vr,vph -> vx,vy ***
        if(rp.gt.0.0d0) then
          vxbp=vrp*xbp/rp-vphp*ybp/rp
          vybp=vrp*ybp/rp+vphp*xbp/rp
        else
          vxbp=0.0d0
          vybp=0.0d0
        endif
c check escape velocity (from Halo)
        if(flagmexth.le.1) then
          r3dp=dsqrt(xbp**2+ybp**2+zbp**2)
          vesc2r=fvesc*dlog(1.0d0+r3dp/rs)/(r3dp/rs)
          v2p=vxbp**2+vybp**2+vzbp**2
          if(v2p.gt.vesc2r) then
            goto 74
          endif
c escape velocity
        endif
      
c *** output 
        if(flagb.eq.0) then
          write(60) xbp,ybp,zbp,vxbp,vybp,vzbp,mbp,rhobp
        else
          write(60,160) xbp,ybp,zbp,vxbp,vybp,vzbp,mbp,rhobp
     &     ,dsqrt(v2p),dsqrt(vesc2r),dsqrt(xbp**2+ybp**2)
 160      format(11(1pE13.5))
        endif
      enddo
      close(60)

      return
      end
