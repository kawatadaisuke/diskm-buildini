c /***************************************************
c   disksgen.f  for dbulhalo
c  1 Feb. 2011  written by D.Kawata
c ****************************************************/
c70*******************************************************************

c /*****   generating exponential disk *****/
      subroutine diskspgen(np,rdlim,zdlim)
      include 'define.f'

      integer np
      character filen*60
      double precision rdlim,zdlim
      integer i,j
      integer nr,nth,nz,ir,ith,iz,izp
      double precision rlim,zlim,dr,dz,dth
      double precision lzi,lzo,dlz,zip,zop
      double precision thp,zp,dvr,rhop,rp
      double precision mtot
      integer nrp,ndrm
      double precision ri,ro,lri,lro,dlr,rr,mrp

c      nr=1024
c      nth=64
c      nz=64
      nr=256
      nth=32
      nz=64

      rlim=rdlim
      zlim=zdlim
      dr=rlim/dble(nr)
      dth=2.0d0*M_PI/dble(nth)
      mtot=0.0d0
c *** set dz by constant d log z
      lzi= dlog10(zdlim*frini0dz)
      lzo = dlog10(zdlim)
      dlz = (lzo-lzi)/(dble(nz/2))

      write(6,*) ' sample disk z range=',10.0d0**lzi,10.0d0**lzo

c *** generate particles
      np=0
      open(60,file='disks.dat',status='unknown')
      do ith=0,nth-1
        thp=0.5d0*dth+dth*dble(ith)
        do iz=0,nz
          if(iz.lt.nz/2) then
            izp=iz
          else
            izp=iz-nz/2
          endif
          zp=10.0d0**(lzi+dble(izp)*dlz)
          if(izp.gt.0) then
            zip=10.0d0**(lzi+dble(izp-1)*dlz)
          else
            zip=0.0d0
          endif
          zop=10.0d0**(lzi+dble(izp+1)*dlz)
          dz=0.5d0*(zop-zip)
          if(iz.lt.nz/2) then
            zp=-zp
          endif         
          do ir=0,nr-1
            rp=dble(ir+1)*dr
            x_p(np)=rp*dcos(thp)
            y_p(np)=rp*dsin(thp)
            z_p(np)=zp
            dvr=dz*(dth*rp)*dr
            dvr_p(np)=dvr
            m_p(np)=0.0d0
            h_p(np)=ETAH*(dvr**(1.0d0/3.0d0))
c
            write(60,'(3(1pE13.5))') x_p(np),y_p(np),z_p(np) 
c
            np=np+1
          enddo
        enddo
      enddo
      close(60)
      write(6,*) 'np sample disk =',np
      write(6,*) 'r,z lim=',rlim,zlim
      write(6,*) 'nr,nth,nz=',nr,nth,nz

      return
      end

c /*****   generating exponential disk *****/
      subroutine disksdenm(np,mdisk,hd,zd)
      include 'define.f'

      integer np,np0,nd,idisk
      character filen*60
      double precision mdisk,hd,zd
      integer i,j
      integer nr,nth,nz,ir,ith,iz
      double precision rlim,zlim,dr,dz,dth
      double precision thp,zp,dvr,rhop,rp
      double precision mtot
      integer nrp,ndrm
      double precision ri,ro,lri,lro,dlr,rr,mrp
      mtot=0.0d0
c *** assign density
      do i=0,np-1
        rp=dsqrt(x_p(i)**2+y_p(i)**2)
c *** density ***
        rhop=(mdisk/(4.0d0*M_PI*zd*(hd**2)))
     &   *(1.0d0/(dcosh(z_p(i)/zd)**2))
     &   *dexp(-rp/hd)
        m_p(i)=m_p(i)+rhop*dvr_p(i)
        mtot=mtot+rhop*dvr_p(i)
      enddo
      write(6,*) ' Mtot sdisk=',mtot

      return
      end
