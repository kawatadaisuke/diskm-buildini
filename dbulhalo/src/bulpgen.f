c /***************************************************
c   bulpgen.f  for dbulhalo
c   9 Feb. 2014  written by D.Kawata
c ****************************************************/
c70*******************************************************************

c /*****   generating exponential disk *****/
      subroutine bulpgen(nb,mbp,ngr,ngz,mbul,ah,rcut,flagb)
      include 'define.f'

      integer nb,ngr,ngz,flagb
      double precision mbul,ah,mbp
      character filen*60
      integer i,j
      integer nrp,ndrm
      double precision rhoc,delc,facnfw,rcut
      double precision ri,ro,lri,lro,dlr,rp,rr,mrp,rp2d
      double precision pfr,ph,th
c *** for particles ***
      double precision xbp,ybp,zbp,vxbp,vybp,vzbp,rhobp
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

c *** external function 
      double precision rmrhp,rhohp
      external rmrhp,rhohp
      integer idum
      real ran1,gasdev
      external ran1,gasdev

      idum = -111

c *** generating particles ***
      if(myrank.eq.0) then
        write(filen,'(a16)') 'output/bulge.dat'
        if(flagb.eq.0) then
          open(60,file=filen,status='unknown',form='unformatted')
          write(60) nb
        else
          open(60,file=filen,status='unknown')
          write(60,'(a3,I10)') '#n=',nb
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

      do i=0,nb-1
 67     pfr = dble(ran1(idum))
c *** search radius ***
        rp=rmrhp(pfr,mbul,ah)
        if(rp.gt.rcut) then
c          write(6,*) ' try again: myrank,pfr,rp>rcut=',myrank,pfr,rp
          goto 67
        endif
c to measure profile
        ir=int((dlog10(rp)-lripr)/dlrpr)
        if(ir.ge.0.and.ir.le.ndpr) then
          mdrpr(ir)=mdrpr(ir)+mbp
        endif
c this does not work. not random in z
c        th =-0.5d0*M_PI+M_PI*dble(ran1(idum))
c        zbp = rp*dsin(th)
        zbp=-rp+2.0d0*rp*dble(ran1(idum))
        ph =  2.0d0*M_PI*dble(ran1(idum))
        rp2d=dsqrt(rp**2-zbp**2)
        mbp = mbp
        xbp = rp2d*dcos(ph)
        ybp = rp2d*dsin(ph)

c *** density ***
        rhobp=rhohp(rp,mbul,ah)
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
        vphp=vphm+vsigph*dble(gasdev(idum))
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
c *** output 
        if(myrank.eq.0) then
          if(flagb.eq.0) then
            write(60) xbp,ybp,zbp,vxbp,vybp,vzbp,mbp,rhobp
          else
            write(60,160) xbp,ybp,zbp,vxbp,vybp,vzbp
     &       ,mbp,rhobp,dsqrt(xbp**2+ybp**2+zbp**2)*LUKPC
     &       ,vsigph,vsigr,vsigz,rp2d
 160        format(13(1pE13.5))
          endif
        endif
      enddo
      close(60)

      if(myrank.eq.0) then
        open(60,file='bulpprof.dat',status='unknown')
        do i=0,ndpr 
          rp=10.0d0**(lripr+dlrpr*(dble(i)+0.5d0))
          dvdr=(4.0d0*M_PI/3.0d0)
     &      *((10.0d0**(lripr+dlrpr*dble(i+1)))**3
     &       -(10.0d0**(lripr+dlrpr*dble(i)))**3)
          write(60,'(3(1pE13.5))') rp*LUKPC,mdrpr(i)/dvdr
     &     ,rhohp(rp,mbul,ah)
        enddo
        close(60)
      endif

      return
      end
