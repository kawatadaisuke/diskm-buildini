c /***************************************************
c   gasdgen.f  for dbulhalo
c   13 Sep. 2012  written by D.Kawata
c ****************************************************/
c70*******************************************************************

c /*****   generating exponential disk *****/
      subroutine gasdgen(np,mtdisk,mdisk,hd,rdlim,zdlim,ngr,ngz,id
     &  ,flagb,dfeh,feh0,dalfe,alfe0,flaggd)
      include 'define.f'
      integer np,nd,ngr,ngz,id,flagb,flaggd
      character filen*60
      double precision mtdisk,mdisk,hd,rdlim,zdlim
      double precision dfeh,feh0,dalfe,alfe0
      integer i,j,pn,pno,pni
      integer nrp,ndrm
      double precision pfr,ph,dr
c *** Exp disk profile ***
      integer irlim
      double precision mrexpd(MNR),rexpd(MNR),rlim
      double precision ri,ro,lri,lro,dlr,rp,mdp,rr,mrp
c *** for particle data ***
      double precision xbp,ybp,zbp,vxbp,vybp,vzbp,mbp,rhobp,pbp
     & ,tmbp,myubp,nhp
c *** to get total mass for z profile ***
      integer ir,iz
      double precision fsum,rhoz,prhoz,pzp,zp,dz,rmz
      integer pnrizi,pnrizo,pnrozi,pnrozo
      double precision wrizi,wrizo,wrozi,wrozo
      double precision wr1,wr2,mzzo,mzzi,dmz
c *** for velocity ***
      double precision vsigz,vsigr,vsigph,vphm
      double precision vrp,vphp
c *** for checking density profile ***
      integer ndrpr,ndzpr,irpr,izpr,nppr(0:MN2D,0:MN2D)
      double precision rlimpr(0:1),zlimpr(0:1)
      double precision logripr,logzipr,dlogrpr,dlogzpr
      double precision rhopr(0:MN2D,0:MN2D),zi,zo,dvpr
c *** external function ***
      double precision mrexpdf,mrexp2dr
      external mrexpdf,mrexp2dr
c *** external function 
      integer idum
      real ran1,gasdev
      external ran1,gasdev

      idum = -111-id
c use mdisk (M(<rdlim)) for particle mass
      mdp=mdisk/dble(np)
      if(myrank.eq.0) then  
        write(6,*) id,' disk particle mass=',mdp
        write(6,*) 'flaggd=',flaggd
      endif

c *** find rlim
      rlim=0.0d0
      irlim=0
      if(flaggd.le.2) then
        do i=0,ngr
          pn=id_rzg(i,0)
          rhoz=rho_rzg(pn)
          if(rhoz.lt.RHO_LIM*1.001d0) then
            irlim=i
            rlim=x_rzg(pn)
            goto 68
          endif
        enddo
      endif
      irlim=ngr
      pn=id_rzg(ngr,0)
      rlim=x_rzg(pn)
c *** for radial profile use analytic formula
  68  ri = 1.0e-6*rdlim
      ro = rdlim
      lri = dlog10(ri)
      lro = dlog10(ro)
      dlr = (lro-lri)/dble(ndr-2)
c *** set mr and r ***
      rexpd(1) = 0.0d0
      mrexpd(1) = 0.0d0
      if(flaggd.le.1) then
        do i=2,ndr
          rexpd(i) = 10.0d0**(lri+dlr*dble(i-2))
          rr = rexpd(i)
          if(rr.lt.rlim) then
c        mrexpd(i) = 1.0d0-(1.0d0+rr/hd)*dexp(-rr/hd)
c mass fraction
            mrexpd(i)=mrexpdf(rr,hd,mtdisk)/mdisk
            if(mrexpd(i).gt.1.0d0) then
              if(myrank.eq.0) then
                write(6,*) ' warning mrexpd>1, i,r,mrexpd=',i,rr
     &            ,mrexpd(i)
              endif
            endif
            irlim=i
          else
            mrexpd(i)=1.0d0
          endif
        enddo
      else if(flaggd.eq.2) then
        do i=2,ndr
          rexpd(i) = 10.0d0**(lri+dlr*dble(i-2))
          rr = rexpd(i)
          if(rr.lt.hd) then
c        mrexpd(i) = 1.0d0-(1.0d0+rr/hd)*dexp(-rr/hd)
            mrexpd(i)=(rr/hd)**2
            irlim=i
          else
            mrexpd(i)=1.0d0
          endif
        enddo
      else
        do i=2,ndr
          rexpd(i) = 10.0d0**(lri+dlr*dble(i-2))
          rr = rexpd(i)
          if(rr.lt.rlim) then
c        mrexpd(i) = 1.0d0-(1.0d0+rr/hd)*dexp(-rr/hd)
c mass fraction
            mrexpd(i)=mrexp2dr(rr)/mdisk
            if(mrexpd(i).gt.1.0d0) then
              if(myrank.eq.0) then
                write(6,*) ' warning mrexpd>1, i,r,mrexpd=',i,rr
     &            ,mrexpd(i)
              endif
            endif
            irlim=i
          else
            mrexpd(i)=1.0d0
          endif
        enddo
      endif

      if(myrank.eq.0) then
        write(6,*) ' rlim,irlim,mr=',rlim,irlim,mrexpd(irlim)
      endif
c *** set M(<z) use _rzg
      do i=0,ngr
        j=0
        pn=id_rzg(i,j)
        mz_rzg(pn)=0.0d0
        rhoz=rho_rzg(pn)
        if(rhoz.lt.RHO_LIM*1.001d0) then
          rhoz=0.0d0
        endif
        zp=z_rzg(pn)
        fsum=0.0d0
        do j=1,ngz
          prhoz=rhoz
          pn=id_rzg(i,j)
          rhoz=rho_rzg(pn)
          if(rhoz.lt.RHO_LIM*1.001d0) then
            rhoz=0.0d0
          endif
          pzp=zp
          zp=z_rzg(pn)
          dz=zp-pzp
          fsum=fsum+dz*0.5d0*(rhoz+prhoz)
          mz_rzg(pn)=fsum
        enddo
c *** normalise ***
        if(fsum.gt.0.0d0) then
          do j=1,ngz
            pn=id_rzg(i,j)
            mz_rzg(pn)=mz_rzg(pn)/fsum        
          enddo
        endif  
      enddo      

      open(60,file='mzgas.dat',status='unknown')
      do i=1,ngr     
        do j=1,ngr     
          pn=id_rzg(i,j)
          write(60,'(4(1pE13.5))') x_rzg(pn),z_rzg(pn),mz_rzg(pn)
     &     ,rho_rzg(pn)
        enddo
      enddo
      close(60)

c*** start generating particles ***
      if(myrank.eq.0) then

      write(filen,'(a12,i1,a4)') 'output/diskc',id,'.dat'
      if(flagb.eq.0) then
        open(60,file=filen,status='unknown',form='unformatted')
        write(60) np,1
        write(60) dfeh,feh0,dalfe,alfe0
      else      
        open(60,file=filen,status='unknown')
        write(60,'(a8,2I8)') '#n,flag=',np,1
        write(60,'(a1,4(1pE13.5))') '#'
     &   ,dfeh,feh0,dalfe,alfe0
      endif

c check profile
      ndrpr=50
      ndzpr=50
      rlimpr(1)=rdlim
      zlimpr(1)=zdlim
      rlimpr(0)=rdlim*1.0e-3
      zlimpr(0)=zdlim*1.0e-3
      dlogrpr=(dlog10(rlimpr(1))-dlog10(rlimpr(0)))/dble(ndrpr)
      dlogzpr=(dlog10(zlimpr(1))-dlog10(zlimpr(0)))/dble(ndzpr)
      logripr=dlog10(rlimpr(0))
      logzipr=dlog10(zlimpr(0))
      write(6,*) ' for gas density profile=',rlimpr(0),rlimpr(1)
     &  ,zlimpr(0),zlimpr(1),dlogrpr,dlogzpr,logripr,logzipr
      do j=0,ndzpr
        do i=0,ndrpr
          rhopr(i,j)=0.0d0
          nppr(i,j)=0
        enddo
      enddo

      do i=0,np-1
 67     pfr = dble(ran1(idum))
        if(pfr.gt.mrexpd(irlim)) then
          goto  67
        endif      
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
          write(6,*) ' Warning: in finding radius for disk: j,mr,prf='
     &      ,j,mrexpd(j),pfr
          pn=id_rzg(irlim,0)
        write(6,*) ' rlim,irlim,mr,rho=',rlim,irlim,mrexpd(irlim)
     &    ,rho_rzg(pn),RHO_LIM*1.001d0
          j = 2
        endif
        if(mrexpd(j)-mrexpd(j-1).gt.0.0d0) then
          rp=10.0d0**(dlog10(rexpd(j-1))
     &    +(dlog10(rexpd(j))-dlog10(rexpd(j-1)))
     &    *(dlog10(pfr)-dlog10(mrexpd(j-1)))
     &    /(dlog10(mrexpd(j))-dlog10(mrexpd(j-1))))
        else
          if(myrank.eq.0) then
            write(6,*) ' cannot find r for pfr=',pfr
          endif
          goto 67
        endif
c        if(myrank.eq.0) then
c          write(6,*) j,rexpd(j-1),rexpd(j),pfr,mrexpd(j-1),mrexpd(j)
c        endif
        rr=rp
        ph =  2.0d0*M_PI*dble(ran1(idum))
        mbp=mdp
        xbp=rp*dcos(ph)
        ybp=rp*dsin(ph)
c *** z-axis ****
 66     pfr = dble(ran1(idum))
c *** set r range ***
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
c *** search z for ir ***
        pni=id_rzg(ir,1)
        pno=id_rzg(ir+1,1)
c *** calculate weight for r direction ***
        dr=x_rzg(pno)-x_rzg(pni)
        wr1=(x_rzg(pno)-rr)/dr
        wr2=(rr-x_rzg(pni))/dr
        mzzo=0.0d0
        do j=0,ngz
          pni=id_rzg(ir,j)
          pno=id_rzg(ir+1,j)
          mzzi=mzzo
          mzzo=wr1*mz_rzg(pni)+wr2*mz_rzg(pno)
          if(mzzo.gt.pfr) then
            goto 73
          endif
        enddo
        if(myrank.eq.0) then
          write(6,*) ' Warning: in finding z for disk'
          write(6,*) '  pfr,r,mzmax=',pfr,rp,mzzo,wr1,wr2
        endif
        goto 66
 73     iz=j-1
        dmz=mzzo-mzzi
        pni=id_rzg(ir,iz)
        pno=id_rzg(ir,iz+1)
        dz=z_rzg(pno)-z_rzg(pni)
c *** z coordinate ***
        zbp=z_rzg(pni)+dz*(pfr-mzzi)/dmz
        zp=zbp
        if(dble(ran1(idum)).lt.0.5d0) then
          zbp=-zbp
        endif
c *** density ***
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
        pnrizi=id_rzg(ir,iz)
        pnrizo=id_rzg(ir,iz+1)
        pnrozi=id_rzg(ir+1,iz)
        pnrozo=id_rzg(ir+1,iz+1)
        dr=x_rzg(pnrozi)-x_rzg(pnrizi)
        dz=z_rzg(pnrizo)-z_rzg(pnrizi)
c weight 
        wrizi=(x_rzg(pnrozi)-rr)*(z_rzg(pnrizo)-zp)/(dr*dz)
        wrizo=(x_rzg(pnrozo)-rr)*(zp-z_rzg(pnrizi))/(dr*dz)
        wroz i=(rr-x_rzg(pnrizi))*(z_rzg(pnrozo)-zp)/(dr*dz)
        wrozo=(rr-x_rzg(pnrizo))*(zp-z_rzg(pnrozi))/(dr*dz)
c *** density ***
        rhobp=wrizi*rho_rzg(pnrizi)+wrizo*rho_rzg(pnrizo)
     &   +wrozi*rho_rzg(pnrozi)+wrozo*rho_rzg(pnrozo)
c *** pressure ***
        pbp=wrizi*p_rzg(pnrizi)+wrizo*p_rzg(pnrizo)
     &   +wrozi*p_rzg(pnrozi)+wrozo*p_rzg(pnrozo)
c *** myu ***
        myubp=wrizi*myu_rzg(pnrizi)+wrizo*myu_rzg(pnrizo)
     &   +wrozi*myu_rzg(pnrozi)+wrozo*myu_rzg(pnrozo)
c *** T ***
        tmbp=wrizi*tm_rzg(pnrizi)+wrizo*tm_rzg(pnrizo)
     &   +wrozi*tm_rzg(pnrozi)+wrozo*tm_rzg(pnrozo)
c *** velocity ***
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
        vzbp=vsigz*dble(gasdev(idum))
c *** convert vr,vph -> vx,vy ***
        if(rp.gt.0.0d0) then
          vxbp=vrp*xbp/rp-vphp*ybp/rp
          vybp=vrp*ybp/rp+vphp*xbp/rp
        else
          vxbp=0.0d0
          vybp=0.0d0
        endif
c *** output 
        if(flagb.eq.0) then
          write(60) xbp,ybp,zbp,vxbp,vybp,vzbp,mbp,rhobp,pbp,myubp
        else
          write(60,160) xbp,ybp,zbp,vxbp,vybp,vzbp,mbp,rhobp,pbp
     &     ,tmbp,myubp,pbp/((GAM-1.0d0)*rhobp),vsigz,vrp
     &     ,((GAM-1.0d0)*TUK/TPRHO)*pbp/((GAM-1.0d0)*rhobp)
 160      format(15(1pE13.5))
        endif
c *** profile
        if(rp.le.rlimpr(0)) then
          irpr=0
        else if(rp.gt.rlimpr(1)) then
          irpr=ndrpr+1
        else
          irpr=int((dlog10(rp)-logripr)/dlogrpr)+1
        endif
        if(dabs(zbp).le.zlimpr(0)) then
          izpr=0
        else if(dabs(zbp).gt.rlimpr(1)) then
          izpr=ndzpr+1
        else
          izpr=int((dlog10(dabs(zbp))-logzipr)/dlogrpr)+1
        endif
        rhopr(irpr,izpr)=rhopr(irpr,izpr)+mbp
        nppr(irpr,izpr)=nppr(irpr,izpr)+1
      enddo

      close(60)

c *** output profile
      open(60,file='gasprof.dat',status='unknown')
      zo=0.0d0
      do j=0,ndzpr
        zp=10.0d0**(logzipr+(dble(j)+0.5d0)*dlogzpr)
        zi=zo
        zo=10.0d0**(logzipr+(dble(j)+1.0d0)*dlogzpr)
        ro=0.0d0        
        do i=0,ndrpr
          rp=10.0d0**(logripr+(dble(i)+0.5)*dlogrpr)
          ri=ro
          ro=10.0d0**(logripr+(dble(i)+1.0d0)*dlogrpr)
          dvpr=(zo-zi)*M_PI*(ro**2-ri**2)
          rhopr(i,j)=0.5d0*rhopr(i,j)/dvpr
          write(60,'(3(1pE13.5),I10)') rp,zp,rhopr(i,j),nppr(i,j)
        enddo
      enddo

c *** only for myrank 0
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end
