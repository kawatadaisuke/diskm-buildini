c /***************************************************
c   exp2dmr.f for dbulhalo
c  5 Jun. 2014  written by D.Kawata
c ****************************************************/
c70*******************************************************************

c calculate and tabulate exp(-R/h-h/R) profile
      subroutine setexp2dmr(mtdisk,mdisk,hd,rdlim)
      include 'define.f'
      double precision mtdisk,mdisk,hd,rdlim
c for work
      integer i
      character filen*60
c *** for radial profile ***
      double precision ri,ro,dr
      double precision rhoro,rhori
      double precision mr,pmr

c integration
      ri=0.0d0
      ro=rdlim
      ndrexp2d=10000
      dr=(ro-ri)/dble(ndrexp2d)
   
      if(myrank.eq.0) then
        write(6,*) ' adjusting mtdisk for exp2d profile'
        write(6,*) ' ri,ro,dr=',ri,ro,dr
      endif

      ro = 0.0d0
      mr = 0.0d0
      rhoro=0.0d0
      rexp2d(0)=0.0d0
      mrexp2d(0)=0.0d0
      rhorexp2d(0)=0.0d0
      do i=1,ndrexp2d
        rhori=rhoro
        pmr = mr
        ri = ro
        ro = ri+dr
        rhoro=dexp(-ro/hd-hd/ro)
        mr = mr+2.0d0*M_PI*(rhoro*ro+rhori*ri)*dr*0.5d0
        mrexp2d(i)=mr
        rexp2d(i)=ro
        rhorexp2d(i)=rhoro
      enddo
      ndrexp2d=ndrexp2d+1

c *** set common values
      roexp2d=ro/LUKPC
      riexp2d=ri/LUKPC
      drexp2d=dr/LUKPC

      if(myrank.eq.0) then
        write(6,*) ' exp2d profile total vs. <Rlim='
     &   ,3.18884*(hd**2),mrexp2d(ndrexp2d-1)
     &   ,3.18884*(hd**2)/mrexp2d(ndrexp2d-1)
      endif
c increase total mass (R<Inf) to M(<Rlim)=mdisk for constant particle mass
      mtdisk=mdisk*3.18884*(hd**2)/mrexp2d(ndrexp2d-1)
c normalisation for M(<R)
      open(60,file='exp2d.dat',status='unknown')
      do i=0,ndrexp2d-1
        mrexp2d(i)=mrexp2d(i)*mtdisk/(3.18884*(hd**2))
        rhorexp2d(i)=mrexp2d(i)*mtdisk/(3.18884*(hd**2))
c change unit
        rexp2d(i)=rexp2d(i)/LUKPC
        mrexp2d(i)=mrexp2d(i)/MUSM
        write(60,'(3(1pE13.5))') rexp2d(i),rhorexp2d(i),mrexp2d(i)
      enddo
      close(60)

      return
      end

      function mrexp2dr(rp)
      include 'define.f'
      double precision mrexp2dr,rp
      integer idr

      idr=int(rp/drexp2d)
      if(idr.ge.0.and.idr.lt.ndrexp2d-1) then
        mrexp2dr=mrexp2d(idr)+(mrexp2d(idr+1)-mrexp2d(idr))
     &   *(rp-rexp2d(idr))/drexp2d
      else if(idr.lt.0) then
        mrexp2dr=0.0d0
      else
        mrexp2dr=mrexp2d(ndrexp2d-1)
      endif

      return
      end
