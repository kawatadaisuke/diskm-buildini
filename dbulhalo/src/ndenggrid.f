c *****************************************************
c    ndenggrid.f
c 27 Jan., 2011   written by D. Kawata
c ***************************************************** 
c70*******************************************************************

      subroutine ndenggrid(ngr,ngz,mdisk,hd,ids,ide,fpjm,mdp
     &  ,zred,dsig,flagexit,flaggd)
      include './define.f'
      integer ngr,ngz,npg,ids,ide,flagexit,flaggd
      double precision fpjm,mdp
      double precision mdisk,hd
      integer i,j,k,pn,pnpi,pnpo,ir,iz,pnp2i,pnrp
      double precision sdenr
      double precision sumd,fn,prhoz,rhoz
c *** for finding the density ***
      integer pnrizi,pnrizo,pnrozi,pnrozo
      double precision wrizi,wrizo,wrozi,wrozo
      double precision rp,zp,dr,dz,mtot,mtotall,rhop
      double precision zmax
      double precision drho,drhomax
      integer npp,npt
c *** smoothing 
      integer ip,jp,ndi
      integer ipi,jpj,pnj
      double precision trho_rzg(0:MN-1),dsig,dx,dy,rp2
      double precision totm1,totm2
c *** getting p ***
      double precision lmetp,lnhp,temp,myup,pp,gamp,phlim
      double precision zred,sdenc,hp

c *** dendity for constant density case
      sdenc=mdisk/(M_PI*(hd**2))

      zmax=0.0d0
      drhomax=0.0d0
      do i=0,ngr
        j=0
        pn=id_rzg(i,j)
        rp=x_rzg(pn)
c *** surface density target ***
        if(flaggd.le.1) then
          sdenr=(mdisk/(2.0d0*M_PI*(hd**2)))
     &      *dexp(-rp/hd)
        else if(flaggd.eq.2) then
          if(rp.lt.hd) then
            sdenr=sdenc
          else
            sdenr=0.0d0
          endif
        else
c exp(-R/h-h/R)
          sdenr=(mdisk/(3.18884*(hd**2)))
     &      *dexp(-rp/hd-hd/rp)
        endif
        sumd=0.0d0
        rhoz=rho_rzg(pn)
        if(rhoz.lt.RHO_LIM*1.0001d0) then
          rhoz=0.0d0
        endif
        do j=1,ngz
          prhoz=rhoz
          pn=id_rzg(i,j)
          pnpi=id_rzg(i,j-1)
          dz=z_rzg(pn)-z_rzg(pnpi)
          rhoz=rho_rzg(pn)
          if(rhoz.lt.RHO_LIM*1.0001d0) then
c           if(j.ge.2) then
c              pnp2i=id_rzg(i,j-2)      
c              pnpi=id_rzg(i,j-1)
c              drho=dlog10(rho_rzg(pnpi))-dlog10(rho_rzg(pnp2i))
c              rho_rzg(pn)=10.0d0**(dlog10(rho_rzg(pnpi))+drho)
c              if(rho_rzg(pn).lt.RHO_LIM*1.0001d0) then
c                pnrp=id_rzg(i-1,j)
c                if(rho_rzg(pnrp).gt.RHO_LIM*1.0001d0) then
c                  rho_rzg(pn)=rho_rzg(pnrp)
c                  rhoz=rho_rzg(pn)
c                else
                  rhoz=0.0d0
c                endif
c              else
c                rhoz=rho_rzg(pn)
c              endif
c            endif
          endif
c *** integration ***
          sumd=sumd+dz*0.5d0*(rhoz+prhoz)
        enddo
c *** normalisation factor
        sumd=sumd*2.0d0
        if(sumd.gt.0.0d0) then
          fn=sdenr/sumd
        else
          fn=1.0d0
        endif
        if(flagexit.eq.0) then
          do j=0,ngz
            pn=id_rzg(i,j)
            if(rho_rzg(pn)*fn.gt.RHO_LIM*1.0001d0) then
              rho_rzg(pn)=rho_rzg(pn)*fn
              if(zmax.lt.z_rzg(pn)) then
                zmax=z_rzg(pn)
              endif
            else
              rho_rzg(pn)=RHO_LIM
            endif
          enddo
        else
          do j=0,ngz
            pn=id_rzg(i,j)
            if(rho_rzg(pn).gt.RHO_LIM*1.0001d0) then
              rho_rzg(pn)=rho_rzg(pn)*fn
              if(zmax.lt.z_rzg(pn)) then
                zmax=z_rzg(pn)
              endif
            else
              rho_rzg(pn)=RHO_LIM
            endif
          enddo
        endif
      enddo
c      open(60,file='ndengas.dat',status='unknown')

c gaussian smoothing
      do i=0,ngr
        do j=0,ngz
          pn=id_rzg(i,j)
          trho_rzg(pn)=rho_rzg(pn)
        enddo
      enddo
c      dsig=3.0d0
      ndi=5*int(dsig)
      do i=0,ngr
        do j=0,ngz
          pn=id_rzg(i,j)
          rho_rzg(pn)=0.0d0
          do ip=i-ndi,i+ndi
            dx=dabs(dble(ip-i))
            ipi=ip
            if(ip.lt.0) then
              ipi=-ip
            else if(ip.gt.ngr) then
              ipi=ngr
            endif
            do jp=j-ndi,j+ndi            
              dy=dabs(dble(jp-j))    
              jpj=jp
              if(jp.lt.0) then
                jpj=-jp
              else if(jp.gt.ngr) then
                jpj=ngz
              endif
              rp2=dx**2+dy**2
              if(dsqrt(rp2).lt.dble(ndi)*1.1d0) then
                pnj=id_rzg(ipi,jpj)
c *** gaussian smoothing of density ***
                rho_rzg(pn)=rho_rzg(pn)+(1.0d0/(2.0d0*M_PI*dsig**2))
     &            *dexp(-rp2/(2.0d0*dsig**2))*trho_rzg(pnj)
              endif
            enddo
          enddo
c *** get p for new rho ***
          rhop=rho_rzg(pn)
          if(fpjm.ge.0.0d0) then          
            lmetp=lmet_rzg(pn)
            lnhp=dlog10(rhop*fh_rzg(pn)*(DU/MP))
            call pgameos(rhop,lmetp,pp,gamp,temp,lnhp,myup)
c from eq. (1) of Hopkins et al. (2011)
            hp=ETAH*(mdp/rhop)**THIRD
            phlim=1.2d0*(fpjm**(2.0d0/3.0d0)*G*(hp**2)
     &        *(rhop**2))/GAM
            if(pp.lt.phlim) then
              pp=phlim
              gamp=4.0d0/3.0d0
              call cool(rhop,pp,lnhp,lmetp,zred,myup)
              call cool(rhop,pp,lnhp,lmetp,zred,myup)
              temp=pp*(myup*TUK)/(rhop*TPRHO*MYU)
            endif
          else
c *** isothermal case ***
            gamp=1.0d0
            pp=rhop*TPRHO*(-fpjm/TUK)
            temp=-fpjm
            myup=0.6d0
          endif
          tm_rzg(pn)=temp
          p_rzg(pn)=pp
          myu_rzg(pn)=myup
        enddo
      enddo

c *** change the mass of the sampled disk particles ***
      mtot=0.0d0
      do i=ids,ide
        rp=dsqrt(x_p(i)**2+y_p(i)**2)
        zp=dabs(z_p(i))
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
c density
        rhop=wrizi*rho_rzg(pnrizi)+wrizo*rho_rzg(pnrizo)
     &   +wrozi*rho_rzg(pnrozi)+wrozo*rho_rzg(pnrozo)
c        if(rho_p(i).lt.RHO_LIM*1.0001d0) then
c          rho_p(i)=0.0d0
c          m_p(i)=0.0d0
c          h_p(i)=hmin_p
c        else
c mass
         m_p(i)=rhop*dvr_p(i)
         mtot=mtot+rhop*dvr_p(i)
      enddo
      close(60)

c adjust h
c      npp=ide-ids+1
c      npt=0
c      call MPI_ALLREDUCE(npp,npt,1,MPI_INTEGER
c     &  ,MPI_SUM,MPI_COMM_WORLD,ierr)
c      mdp=mdisk/dble(npt)
c      if(myrank.eq.0) then
c        write(6,*) ' in ndenggrid mdp=',mdp
c      endif

c      do i=ids,ide
c        h_p(i)=h_p(i)*((m_p(i)/mdp)**(1.0d0/3.0d0))
c        if(h_p(i).lt.eps) then
c          h_p(i)=eps
c        endif
c      enddo
c
c 
      mtotall=0.0d0
      call MPI_ALLREDUCE(mtot,mtotall,1,MPI_DOUBLE_PRECISION
     &  ,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myrank.eq.0) then
        write(6,*) ' M, zmax for gas disk=',mtotall,zmax
      endif

      return
      end
