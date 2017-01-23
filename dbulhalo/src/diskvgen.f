c /***************************************************
c   diskvgen.f  for dbulhalo
c  23 Jan. 2011  written by D.Kawata
c ****************************************************/
c70*******************************************************************

c /*****   generating exponential disk *****/
      subroutine diskvgen(ids,ide,ngr,ngz,id)
      include 'define.f'

      integer ids,ide,ngr,ngz,id
      integer i,ir,iz
      integer pnrizi,pnrizo,pnrozi,pnrozo
      double precision wrizi,wrizo,wrozi,wrozo
      double precision rp,zp,dr,dz
      double precision vsigz,vsigr,vsigph,vphm
      double precision vrp,vphp
      character filen*60
c *** external function 
      integer idum
      real ran1,gasdev
      external ran1,gasdev

      idum = -111

c      write(filen,'(a6,i1.1,a4)') 'diskvc',id,'.dat'
c      open(60,file=filen,status='unknown') 

      do i=ids,ide
c *** radius and abs(z) ***
        rp=dsqrt(x_bp(i)**2+y_bp(i)**2)
        zp=dabs(z_bp(i))
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
        vphp=vphm+vsigph*dble(gasdev(idum))
        vrp=vsigr*dble(gasdev(idum))
        vz_bp(i)=vsigz*dble(gasdev(idum))
c *** convert vr,vph -> vx,vy ***
        if(rp.gt.0.0d0) then
          vx_bp(i)=vrp*x_bp(i)/rp-vphp*y_bp(i)/rp
          vy_bp(i)=vrp*y_bp(i)/rp+vphp*x_bp(i)/rp
        else
          vx_bp(i)=0.0d0
          vy_bp(i)=0.0d0
        endif
c *** output ***
c        write(60,160) x_bp(i),y_bp(i),z_bp(i)
c     &   ,vx_bp(i),vy_bp(i),vz_bp(i)
c     &   ,rp,vphm,vsigph,vsigr,vsigz
c     &   ,dabs(z_bp(i)),vphp
c     &   ,vph_rzg(pnrizi),vph_rzg(pnrizo)
c     &   ,vph_rzg(pnrozi),vph_rzg(pnrozo)
c     &   ,x_rzg(pnrizi),x_rzg(pnrizo)
c     &   ,x_rzg(pnrozi),x_rzg(pnrozo)
c     &   ,pnrizi,pnrizo,pnrozi,pnrozo,ir,iz
c 160    format(21(1pE13.5),6I6)
c 160    format(13(1pE13.5))
      enddo
c      close(60)

      return
      end
