c /***************************************************
c   disksgen.f  for dbulhalo
c  1 Feb. 2011  written by D.Kawata
c ****************************************************/
c70*******************************************************************

c /*****   generating exponential disk *****/
      subroutine diskspgen(np,rdlim,zdlim,flaglog,id,nr,nth,nz)
      include 'define.f'

      integer np,flaglog,id
      character filen*60
      double precision rdlim,zdlim
      integer i,j
      integer nr,nth,nz,ir,ith,iz,isz,iez,izp
      double precision rlim,zlim,dr,dz,dth,dmax
      double precision thp,zp,dvr,rhop,rp
      double precision mtot
      integer nrp,ndrm
      double precision ri,ro,lri,lro,dlr,rr,mrp
      double precision lzi,lzo,dlz,zip,zop

c      nr=1024
c      nth=64
c      nz=64
c      nr=256
c      nth=32
c      nz=64

      if(nprocs.gt.nth) then
        if(myrank.eq.0) then 
          write(6,*) ' Error: nprocs>nth=',nth
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE()
        stop
      endif

      rlim=rdlim
      zlim=zdlim
      dr=rlim/dble(nr)
      dz=2.0d0*zlim/dble(nz)
      dth=2.0d0*M_PI/dble(nth)
      mtot=0.0d0
      if(flaglog.ne.0) then
        lzi=dlog10(zdlim*frini0dz)
c 1 pc?
c        lzi=-1.0e-5
        lzo=dlog10(zdlim)
        dlz=(lzo-lzi)/dble(nz/2-1)
        if(myrank.eq.0) then
          write(6,*) id,' log sample disk z range='
     &    ,10.0d0**lzi,10.0d0**lzo
        endif
      endif

c *** generate particles
c      open(60,file='disks.dat',status='unknown')
c get range of ith for each proc
      call para_range(0,nth-1,nprocs,myrank,jsta,jend)

      if(flaglog.ne.0) then
        isz=0
        iez=nz-1
      else
        isz=-nz/2
        iez=nz/2
      endif
      do ith=jsta,jend
        thp=0.5d0*dth+dth*dble(ith)
        do iz=isz,iez
          if(flaglog.ne.0) then
            if(iz.lt.nz/2) then
              izp=iz
            else
              izp=iz-nz/2
            endif
            zp=10.0d0**(lzi+dble(izp)*dlz)
            if(izp.gt.0) then
              zip=10.0d0**(lzi+dble(izp-1)*dlz)
            else
              zip=-10.0d0**(lzi)
            endif
            zop=10.0d0**(lzi+dble(izp+1)*dlz)
            dz=0.5d0*(zop-zip)
            if(iz.lt.nz/2) then
              zp=-zp
            endif
          else
            zp=-dble(iz)*dz
          endif
          do ir=0,nr-1
            rp=dble(ir+1)*dr
            x_p(np)=rp*dcos(thp)
            y_p(np)=rp*dsin(thp)
            z_p(np)=zp
            dvr=dz*(dth*rp)*dr
            dvr_p(np)=dvr
            m_p(np)=0.0d0
            dmax=rp*dth
            if(dz.gt.dmax) then
              dmax=dz
            endif
            if(dr.gt.dmax) then
              dmax=dr
            endif
c            h_p(np)=ETAH*(dvr**(1.0d0/3.0d0))
            h_p(np)=ETAH*dmax
            np=np+1
          enddo
        enddo
      enddo
 
      if(np.gt.MNS-1) then
        write(6,*) ' np > MNS np,myrank,MNS=',np,myrank,MNS
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE()
        stop
      endif

      return
      end

c /*****   generating exponential disk *****/
      subroutine disksdenm(ids,ide,mtdisk,hd,zd0,id,flaggd
     &  ,flagflzd,Rflzd)
      include 'define.f'

      integer ids,ide,nd,idisk,id,flaggd,flagflzd
      character filen*60
      double precision mtdisk,hd,zd0,Rflzd
      integer i,j
      integer nr,nth,nz,ir,ith,iz
      double precision rlim,zlim,dr,dz,dth,zdr
      double precision thp,zp,dvr,rhop,rp
      double precision mtot,mtotall
      integer nrp,ndrm,npp,npt
      double precision ri,ro,lri,lro,dlr,rr,mrp
      double precision rhoc

      mtot=0.0d0

      npp=ide-ids+1
      npt=0
      call MPI_ALLREDUCE(npp,npt,1,MPI_INTEGER
     &  ,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myrank.eq.0) then
        write(6,*) id,' disk Ntot,sample=',npt
      endif
c *** for constant density if flaggd>1
      rhoc=mtdisk/(M_PI*(hd**2)*(2.0d0*zd0))

c *** assign density
      do i=ids,ide
        rp=dsqrt(x_p(i)**2+y_p(i)**2)
c *** density ***
        if(flaggd.le.1) then
          if(flagflzd.eq.0) then
            rhop=(mtdisk/(4.0d0*M_PI*zd0*(hd**2)))
     &       *(1.0d0/(dcosh(z_p(i)/zd0)**2))
     &       *dexp(-rp/hd)
          else
c flaring
            zdr=zd0*dexp(rp/Rflzd)
            rhop=(mtdisk/(4.0d0*M_PI*zdr*(hd**2)))
     &       *(1.0d0/(dcosh(z_p(i)/zdr)**2))
     &       *dexp(-rp/hd)
          endif
        else if(flaggd.eq.2) then
c *** constant density case ***
          if(rp.lt.hd.and.dabs(z_p(i)).lt.zd0) then
            rhop=rhoc
          else
            rhop=RHO_LIM
          endif
        else
          rhop=(mtdisk/(3.18884*2.0d0*zd0*(hd**2)))
     &     *(1.0d0/(dcosh(z_p(i)/zd0)**2))
     &     *dexp(-rp/hd-hd/rp)
        endif
        m_p(i)=m_p(i)+rhop*dvr_p(i)
        mtot=mtot+rhop*dvr_p(i)
      enddo
      mtotall=0.0d0
      call MPI_ALLREDUCE(mtot,mtotall,1,MPI_DOUBLE_PRECISION
     &  ,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myrank.eq.0) then
        write(6,*) id,' disk Mtot=',mtotall
      endif

      return
      end
