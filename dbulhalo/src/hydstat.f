c *****************************************************
c    hydstat.f
c 15 Mar., 2013   written by D. Kawata
c ***************************************************** 
c70*******************************************************************
      subroutine hydstat(ngr,ngz,errmax,errm,nit,fpjm,mdp,zred)
      include './define.f'
      integer ngr,ngz,npg
      integer i,j,k,pn,pnpi,pnpo,nit,nerr
      double precision fpjm,mdp,zred
      integer pnrp
      character filen*60
      double precision gamp,pp,dz,errmax,temp,myup,rho0,errm
      double precision hp
      double precision phlim
      double precision pgamp,ppp,ptemp
      double precision rpmax,zpmax,rhomax,rho0max
      double precision pdrho,drho,lnhp,lmetp,errp
      double precision rhopi,rhopp,rhocc
      double precision lnhpi,tmpi
c /*** Solar Abundance (meteorites from WW95) in terms of mass ***/
c /* 1  2   3   4   5   6    7    8    9    10   11  12   13   14   15
c /* 1H,4He,12C,14N,16O,20Ne,23Na,24Mg,27Al,28Si,32S,36Ar.40Ca,56Fe,58Ni */
      double precision xsol(15),XSZ
      data xsol/0.706d0,0.275d0,3.03d-3,1.11d-3,9.59d-3,1.62d-3,
     & 3.34d-5,5.15d-4,5.81d-5,6.53d-4,3.96d-4,7.74d-5,5.99d-5,
     & 1.17d-3,4.94d-5/
      parameter (XSZ=0.019d0)

      errmax=0.0d0
      errm=0.0d0
      nerr=0
c
c      if(myrank.eq.0) then
c        write(filen,'(a5,i3.3)') 'hydst',nit
c        open(60,file=filen,status='unknown')
c      endif
c
      do i=ngr,0,-1
        do j=1,ngz+1
          errp=0.0d0
          pnpi=id_rzg(i,j-1)
          if(j.le.ngz) then
            pn=id_rzg(i,j)
            dz=z_rzg(pn)-z_rzg(pnpi)
            rho0=rho_rzg(pn)
          endif
c *** store initial metallicity and density
          lmetp=lmet_rzg(pnpi)
c *** use properties at j-1 and integrate
          rhopi=rho_rzg(pnpi)
          if(fpjm.ge.0.0d0) then
c *** N_H for j-1
            lnhp=dlog10(rhopi*fh_rzg(pnpi)*(DU/MP))
            lnhpi=lnhp
c *** (drho/dz)_j-1
            call pgameos(rhopi,lmet_rzg(pnpi),pp,gamp,temp,lnhp,myup)
c constant gamma
c            gamp=4.0d0/3.0d0
c from eq. (1) of Hopkins et al. (2011)
            hp=ETAH*(mdp/rhopi)**THIRD
            phlim=1.2d0*(fpjm**(2.0d0/3.0d0)*G*(hp**2)
     &        *(rhopi**2))/GAM
            if(pp.lt.phlim) then
              pp=phlim
              gamp=4.0d0/3.0d0
              call cool(rhopi,pp,lnhp,lmet_rzg(pnpi),zred,myup)
              call cool(rhopi,pp,lnhp,lmet_rzg(pnpi),zred,myup)
              temp=pp*(myup*TUK)/(rhopi*TPRHO*MYU)
            endif
          else
c *** isothermal case ***
            gamp=1.0d0
c fpjm is negative temperature of the gas disk
            pp=rhopi*TPRHO*(-fpjm/TUK)
            temp=fpjm
            myup=0.6d0
          endif
          pgamp=gamp
          ppp=pp
          ptemp=temp
          pdrho=-((rhopi**2)/(gamp*pp))*(-fz_rzg(pnpi))
c *** set temperature and gamma for j-1
          tm_rzg(pnpi)=temp
          p_rzg(pnpi)=pp
          myu_rzg(pnpi)=myup
          if(j.le.ngz) then
c *** get predicted density for j ***
            rho_rzg(pn)=rhopi+dz*pdrho
            if(rho_rzg(pn).lt.RHO_LIM) then
              rho_rzg(pn)=RHO_LIM
            endif
c *** (drho/dz)_j
            rhopp=rho_rzg(pn)
            lmetp=lmet_rzg(pn)
            if(fpjm.ge.0) then
              call pgameos(rhopp,lmet_rzg(pn),pp,gamp,temp,lnhp,myup)
c constant gamma
c            gamp=4.0d0/3.0d0
c from eq. (1) of Hopkins et al. (2011)
              hp=ETAH*(mdp/rhopp)**THIRD
              phlim=1.2d0*(fpjm**(2.0d0/3.0d0)*G*(hp**2)
     &          *(rhopp**2))/GAM
              if(pp.lt.phlim) then
                pp=phlim
                gamp=4.0d0/3.0d0
c since this temp is not used later, no need to call cool again
              endif
            else
c *** isothermal case ***
              gamp=1.0d0
              pp=rhopp*TPRHO*(-fpjm/TUK)
              temp=-fpjm
              myup=0.6d0
            endif
            drho=-((rhopp**2)/(gamp*pp))*(-fz_rzg(pn))
c *** corrector ***
c *** get predicted density for j ***
            rho_rzg(pn)=rhopi+dz*0.5d0*(drho+pdrho)
            if(rho_rzg(pn).lt.RHO_LIM) then
              if(i.lt.ngr.and.j.lt.ngz) then
                pnrp=id_rzg(i+1,j)
                if(rho_rzg(pnrp).gt.RHO_LIM*1.0001d0) then
                  rho_rzg(pn)=10.0d0**(0.5d0*(dlog10(rho_rzg(pnrp))
     &              +dlog10(RHO_LIM)))
                else
                  rho_rzg(pn)=RHO_LIM
                endif
              else
               rho_rzg(pn)=RHO_LIM
              endif
            endif
            lnhp=dlog10(rho_rzg(pn)*fh_rzg(pn)*(DU/MP))
            errp=dabs(rho0-rho_rzg(pn))/rho0
            if(rho_rzg(pn).gt.RHO_LIM*100.0d0
     &        .and.rho0.gt.RHO_LIM*100.0d0) then
c *** ignore central region
            if(x_rzg(pn).gt.0.0001d0.and.x_rzg(pn).lt.0.3d0) then
            if(z_rzg(pn).gt.0.0001d0.and.z_rzg(pn).lt.0.05d0) then
              if(errp.gt.errmax) then
                errmax=errp
                rpmax=x_rzg(pn)
                zpmax=z_rzg(pn)
                rhomax=rho_rzg(pn)
                rho0max=rho0
              endif
              errm=errm+errp
              nerr=nerr+1
            endif
            endif
           endif
c check
c            if(myrank.eq.0) then
c              temp=tm_rzg(pnpi)
c              write(60,'(18(1pE13.5),I8)') x_rzg(pn),z_rzg(pn)
c     &         ,rho_rzg(pn),rho0,lnhp,errp,dlog10(temp)
c     &         ,gamp,pp,fz_rzg(pn)
c     &         ,pdrho,pgamp,ppp,drho,rhopi,ptemp
c     &         ,x_rzg(pnpi),z_rzg(pnpi),pnpi
c            endif
          endif
        enddo
      enddo
      if(nerr.gt.0) then
        errm=errm/dble(nerr)
      endif

      if(myrank.eq.0) then
        write(6,*) ' err max at r,z,rho,rho0=',rpmax,zpmax
     &    ,rhomax,rho0max
      endif

c *** set rho(r=0)=rho(1)
c      i=1
c      do j=0,ngz
c        pnpi=id_rzg(i,j)
c        pn=id_rzg(i+1,j)
c        rho_rzg(pnpi)=rho_rzg(pn)
c        lmet_rzg(pnpi)=lmet_rzg(pn)
c        fh_rzg(pnpi)=fh_rzg(pn)
c        tm_rzg(pnpi)=tm_rzg(pn)
c        p_rzg(pnpi)=p_rzg(pn)
c      enddo

      return
      end
