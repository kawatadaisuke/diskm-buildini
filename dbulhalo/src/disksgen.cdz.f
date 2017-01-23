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
      integer nr,nth,nz,ir,ith,iz
      double precision rlim,zlim,dr,dz,dth
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
      dz=2.0d0*zlim/dble(nz)
      dth=2.0d0*M_PI/dble(nth)
      mtot=0.0d0

c *** generate particles
      np=0
c      open(60,file='disks.dat',status='unknown')
      do ith=0,nth-1
        thp=0.5d0*dth+dth*dble(ith)
        do iz=-nz/2,nz/2
          zp=-dble(iz)*dz
          do ir=0,nr-1
            rp=dble(ir+1)*dr
            x_p(np)=rp*dcos(thp)
            y_p(np)=rp*dsin(thp)
            z_p(np)=zp
            dvr=dz*(dth*rp)*dr
            dvr_p(np)=dvr
            m_p(np)=0.0d0
            h_p(np)=ETAH*(dvr**(1.0d0/3.0d0))
            np=np+1
          enddo
        enddo
      enddo
      write(6,*) 'generated sample disk particle np=',np
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
      write(6,*) ' total mass in sample particles=',mtot

      return
      end
