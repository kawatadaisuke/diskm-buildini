c /***************************************************
c   dmhpgen.f  for dbulhalo
c   16 Jun. 2011  written by D.Kawata
c ****************************************************/
c70*******************************************************************

c /*****   generating exponential disk *****/
      subroutine dmhpgen(nh,mdmhp,ngr,ngz,omg0,omz,hub,cnfw
     &      ,rvir,zt,rcutdmh,rt,fb,flagb,flagmexth)
      include 'define.f'

      integer nh,ngr,ngz,flagb,flagmexth
      double precision omg0,omz,hub,cnfw,rvir,zt,rcutdmh,rt,fb
      character filen*60
      integer i,j
      integer nrp,ndrm
      double precision rhoc,delc,facnfw
      double precision ri,ro,lri,lro,dlr,rp,mdmhp,rr,mrp,rp2d
      double precision pfr,ph,th
c *** for particles ***
      double precision xbp,ybp,zbp,vxbp,vybp,vzbp,mbp,rhobp
c *** for velocity ***
      integer ir,iz
      integer pnrizi,pnrizo,pnrozi,pnrozo
      double precision wrizi,wrizo,wrozi,wrozo
      double precision dr,dz
      double precision vsigz,vsigr,vsigph,vphm
      double precision vrp,vphp

c to measure profiles 
      integer ndpr,ipr
      double precision lripr,lropr,mdrpr(0:MNR),dlrpr,dvdr
c for escape velocity 
      double precision fvesc,vesc2r,v2p,r3dp,rs

c *** external function 
      double precision rhonfwsw
      external rhonfwsw
      integer idum
      real ran1,gasdev
      external ran1,gasdev

c *** NFW profile constant
      rhoc = (3.0d0*HUB0*HUB0*hub*hub/(8.0d0*M_PI*G))
     &  *(omg0/omz)*((1.0d0+zt)**3)
      delc = (200.0d0/3.0d0)*(cnfw**3/(dlog(1.0d0+cnfw)
     &  -cnfw/(1.0d0+cnfw)))
      facnfw = delc*rhoc*(1.0d0-fb)
      rs=rvir/cnfw
c escape velocity constant from BT2nd p.71
      fvesc=2.0d0*4.0d0*M_PI*G*(1.0d0-fb)*delc*rhoc*(rs**2)
      r3dp=0.1d0
      if(myrank.eq.0) then
        write(6,*) ' in dmhpgen rhoc,delc,fb=',rhoc,delc,fb
        write(6,*) ' esc velocity at 10 kpc=',VUKMS*dsqrt(fvesc
     &    *dlog(1.0d0+r3dp/rs)/(r3dp/rs))
      endif

      idum = -111

c *** generating particles ***
      if(myrank.eq.0) then
        write(filen,'(a15)') 'output/halo.dat'
        if(flagb.eq.0) then
          open(60,file=filen,status='unknown',form='unformatted')
          write(60) nh
        else
          open(60,file=filen,status='unknown')
          write(60,'(a3,I10)') '#n=',nh
        endif
      endif

c for profile
      lripr=dlog10(0.01d0)
      lropr=dlog10(0.6d0)
      ndpr=10
      dlrpr=(lropr-lripr)/dble(ndpr)
      do i=0,ndpr
        mdrpr(i)=0.0d0
      enddo

      do i=0,nh-1
 67     pfr = dble(ran1(idum))
c /*** search radius ***/
        do j=1,ndr
          if(mrdmh(j)/mrdmh(ndr).gt.pfr) then
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
        pfr=pfr*mrdmh(ndr)
        rp=10.0d0**(dlog10(rdmh(j-1))
     &    +(dlog10(rdmh(j))-dlog10(rdmh(j-1)))
     &    *(dlog10(pfr)-dlog10(mrdmh(j-1)))
     &    /(dlog10(mrdmh(j))-dlog10(mrdmh(j-1))))
        if(rp.gt.rcutdmh) then
          goto 67
        endif
c to measure profile
        ir=int((dlog10(rp)-lripr)/dlrpr)
        if(ir.ge.0.and.ir.le.ndpr) then
          mdrpr(ir)=mdrpr(ir)+mdmhp
        endif
c this does not work. not random in z
c        th =-0.5d0*M_PI+M_PI*dble(ran1(idum))
c        zbp = rp*dsin(th)
        zbp=-rp+2.0d0*rp*dble(ran1(idum))
        ph =  2.0d0*M_PI*dble(ran1(idum))
        rp2d=dsqrt(rp**2-zbp**2)
        mbp = mdmhp
        xbp = rp2d*dcos(ph)
        ybp = rp2d*dsin(ph)

c *** density ***
        rhobp=rhonfwsw(rp,cnfw,rvir,rt,flagmexth)*facnfw
c *** velocity ***
c *** find the point in r_rzg ***
        if(rp2d.gt.0.0d0) then
          ir=int((dlog10(rp2d)-lri_rzg)/dlr_rzg)+1
        else
          ir=0
        endif
        if(ir.lt.0) then
          ir=0
        else if(ir.gt.ngr-1) then
          ir=ngr-1
        endif
c *** find the point in z_rzg ***
        if(dabs(zbp).gt.0.0d0) then
          iz=int((dlog10(dabs(zbp))-lzi_rzg)/dlz_rzg)+1
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
        wrizi=(x_rzg(pnrozi)-rp2d)*(z_rzg(pnrizo)-dabs(zbp))/(dr*dz)
        wrizo=(x_rzg(pnrozo)-rp2d)*(dabs(zbp)-z_rzg(pnrizi))/(dr*dz)
        wrozi=(rp2d-x_rzg(pnrizi))*(z_rzg(pnrozo)-dabs(zbp))/(dr*dz)
        wrozo=(rp2d-x_rzg(pnrizo))*(dabs(zbp)-z_rzg(pnrozi))/(dr*dz)


c *** vphm
        vphm=0.0d0
c *** vsigph
        vsigph=wrizi*vsigph_rzg(pnrizi)+wrizo*vsigph_rzg(pnrizo)
     &   +wrozi*vsigph_rzg(pnrozi)+wrozo*vsigph_rzg(pnrozo)
        if(vsigph.lt.0.0d0) then
          vsigph=0.0d0
        endif

c        if(xbp.gt.0.01d0.and.xbp.lt.0.1d0) then
c          if(zbp.gt.0.01d0.and.zbp.lt.0.1d0) then
c          if(dabs(vsigph-vsigph_rzg(pnrizi))/vsigph_rzg(pnrizi)
c     &      .lt.0.6) then
c        write(6,*) vsigph,ir,iz,rp2d,zbp,dr,dz
c     &     ,x_rzg(pnrozi),x_rzg(pnrizi),z_rzg(pnrizo),z_rzg(pnrizi)
c     &     ,vsigph_rzg(pnrizi),vsigph_rzg(pnrizo),vsigph_rzg(pnrozi)
c     &     ,vsigph_rzg(pnrozo),wrizi,wrizo,wrozi,wrozo
c          endif
c          endif
c        endif
c        if(vsigph.lt.1.0e-5) then
c          write(6,*) vsigph,ir,iz,rp2d,zbp,dr,dz
c     &     ,x_rzg(pnrozi),x_rzg(pnrizi),z_rzg(pnrizo),z_rzg(pnrizi)
c        endif

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
        if(rp2d.gt.0.0d0) then
          vxbp=vrp*xbp/rp2d-vphp*ybp/rp2d
          vybp=vrp*ybp/rp2d+vphp*xbp/rp2d
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
        if(myrank.eq.0) then
          if(flagb.eq.0) then
            write(60) xbp,ybp,zbp,vxbp,vybp,vzbp,mbp,rhobp
          else
            write(60,160) xbp,ybp,zbp,vxbp,vybp,vzbp
     &       ,mbp,rhobp,dsqrt(xbp**2+ybp**2+zbp**2)*LUKPC
     &       ,vsigph,vsigr,vsigz,rp2d
     &       ,dsqrt(v2p),dsqrt(vesc2r)
 160        format(15(1pE13.5))
          endif
        endif
      enddo
      close(60)

      if(myrank.eq.0) then
        open(60,file='dmhpprof.dat',status='unknown')
        do i=0,ndpr 
          rp=10.0d0**(lripr+dlrpr*(dble(i)+0.5d0))
          dvdr=(4.0d0*M_PI/3.0d0)
     &      *((10.0d0**(lripr+dlrpr*dble(i+1)))**3
     &       -(10.0d0**(lripr+dlrpr*dble(i)))**3)
          write(60,'(3(1pE13.5))') rp*LUKPC,mdrpr(i)/dvdr
     &     ,rhonfwsw(rp,cnfw,rvir,rt,flagmexth)*facnfw
        enddo
        close(60)
      endif

      return
      end
