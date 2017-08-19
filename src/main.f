c /***************************************
c      main.f  ver.14
c   to make an initial condition for the merger
c  23 Mar. 2015    written by D.KAWATA
c ***************************************/ 
c70*******************************************************************

      include 'define.f'
      integer MNGAL,MND
      parameter (MNGAL=2,MND=10)
      integer i,j,ig,ngal,i1,i2,id,im,npaskip
      integer flagascii,flagbi,flagrot(MNGAL),flagini
      integer nmgr,nejmgr(MNZMGR),nsniimgr,nzmgr
      integer ndisk(MNGAL),nbul(MNGAL),nhalo(MNGAL),nhgas(MNGAL)
c *** disk is the gas ***
      integer flagd(MNGAL,MND),flagbul(MNGAL),flagid(MNGAL,MND)
c *** for orbit ***
      double precision ecc,rperi,rini,rapo
      double precision vxini(MNGAL),vyini(MNGAL)
      double precision xini(MNGAL),zini(MNGAL)
      double precision etot,ltot,vx12,vy12,r12
      double precision igal(MNGAL),omgal(MNGAL),tdgal(2,MNGAL,MND)
     &  ,mgalt(MNGAL),myum,tbgal(2,MNGAL)
     &  ,metbgal(MNGAL),tagegal(2,MNGAL,MND),alfeb(MNGAL)
      double precision metgal(4,MNGAL,MND),metp,fehp,fehpmax
      double precision metpc,metpo,metpfe,vrotp
      double precision alfe(3,MNGAL,MND),alfep
      double precision mhalo,mbul,mdisk,mhgas
      integer mwnh(MNGAL),ntmp,mwnd(MNGAL,MND)
     &  ,mwnb(MNGAL),mwnhg(MNGAL)
      integer ngd,ndmd,nsd
c halo hot gas abundances
      double precision fehhg(MNGAL),alfehg(MNGAL)
c for orbit: flagini=1 smashing model, flagini=2 setting position and velocity for Gal2
      double precision bimpz,v0x
      double precision x0,y0,z0,th,vx0,vy0,vz0
c *** gas data ***
      integer np,ng,flagfdg
      double precision xp,yp,zp,vxp,vyp,vzp
     &  ,mpart,rhop,ug,hg,mzHeg,mzCg,mzNg
     &  ,mzOg,mzNeg,mzMgg,mzSig,mzFeg,mzZg
     &  ,mzHg,fehg,pg,nhp,mcg,myup,uthp
      double precision tx,ty,tz,rg,r2dg
c *** star data ***
      integer ns,flagfds
      double precision mzHes,mzCs
     &  ,mzNs,mzOs,mzNes,mzMgs,mzSis
     &  ,mzFes,mzZs,ages,rs,mzHs,fehs
      integer iz1
      double precision tepmgr(MNMGR),wz2,wz1
c *** for velocity correction ***
      integer flagvc
      double precision vxgr,vygr,vzgr
c *** for metallicity the end time of mass group of stars
      double precision zmet_smgr(0:MNZMGR)
     & ,te_smgr(MNMGR,0:MNZMGR-1)
      double precision llz_smgr,dlz_smgr,luz_smgr,zmetp
c *** DM data ***
      integer ndm,itmp
c *** number of particles within each galaxies ***
      integer npg(MNGAL),ndmg(MNGAL),nsg(MNGAL)
c /*** gamma: specific heat ****/
      double precision gam
c momentum and centre of the mass calculation  
      integer flagvch
      double precision tmassd(MND),cmxd(MND),cmyd(MND)
     &  ,cmzd(MND),vmxd(MND),vmyd(MND),vmzd(MND)
      double precision tmassh(MNGAL),cmxh(MNGAL)
     &  ,cmyh(MNGAL),cmzh(MNGAL),vmxh(MNGAL),vmyh(MNGAL),vmzh(MNGAL)
c /*** Initial parameter ***/
      double precision CEPS
c /* eps = pow(CEPS*rdm,1/3) */      
      parameter (CEPS=200.0d0)
c /*** for nomalization ***/
      double precision GCLKPC,GCVKMS
      parameter (GCLKPC=4.5d0)
      parameter (GCVKMS=220.0d0)
c /*** for Galactics D model 220 / 1.2 ***/
      double precision GCDVKMS
      parameter (GCDVKMS=220.0d0)
c /*** Solar Abundance (meteorites from WW95) in terms of mass ***/
c /* 1  2   3   4   5   6    7    8    9    10   11  12   13   14   15
c /* 1H,4He,12C,14N,16O,20Ne,23Na,24Mg,27Al,28Si,32S,36Ar.40Ca,56Fe,58Ni */
      double precision xsol(15),XSZ
      data xsol/0.706d0,0.275d0,3.03d-3,1.11d-3,9.59d-3,1.62d-3,
     & 3.34d-5,5.15d-4,5.81d-5,6.53d-4,3.96d-4,7.74d-5,5.99d-5,
     & 1.17d-3,4.94d-5/
      parameter (XSZ=0.019d0)
      character filen*20,ctmp*10,modname*20,cread*200
c *** for random number ***
      integer idum
      real ran1,gasdev
      external ran1,gasdev

      gam=5.0d0/3.0d0
      idum=-111
      ngal = 2
      npaskip=1
c *** 0: output ascii file ***
      flagascii=0
c/*****   Open Initial Data File ./ini/input.dat   *****/
      open(50,file='ini/input.dat',status='old')      
      read(50,*) flagbi,flagascii
      if(flagbi.ne.0) then
        write(6,*) ' flagbi has to be 0 and input data must be binary'
        stop
      endif
c *** number of galaxies ***
      read(50,*) ngal
      if(ngal.gt.2.or.ngal.lt.1) then
        write(6,*) ' number of galxies should be one or two'
        stop
      endif
c *** orbital parameter ***
      read(50,*) flagini
      if(flagini.eq.0) then
        read(50,*) ecc,rperi,rini
      else if(flagini.eq.1) then
        read(50,*) v0x,bimpz,rini
      else if(flagini.eq.2) then
        read(50,*) x0,z0,vx0,vy0
      else
        write(6,*) ' flagini has to be 0-2'
        stop
      endif
      do ig = 1,ngal
c *** since ver.2 only allow to use output of dbulhalo
        read(50,*) ndisk(ig),nbul(ig),nhalo(ig),nhgas(ig)
c /*** inclination angles with degree ***/
        read(50,*) flagrot(ig),flagvc 
        if(flagvc.ne.0) then
          write(6,*) ' Error: flagvc has to be 0'
          write(6,*) ' no V correction available after ver.9'
          stop
        endif
        read(50,*) igal(ig),omgal(ig)
        read(50,*) flagbul(ig)
        read(50,*) tbgal(1,ig),tbgal(2,ig)
        read(50,*) metbgal(ig),alfeb(ig)
        read(50,*) fehhg(ig),alfehg(ig)
        do id=1,ndisk(ig)
          read(50,*) flagd(ig,id)
          read(50,*) metgal(1,ig,id),metgal(2,ig,id)
     &     ,metgal(3,ig,id),metgal(4,ig,id)
          read(50,*) alfe(1,ig,id),alfe(2,ig,id),alfe(3,ig,id)
          read(50,*) tdgal(1,ig,id),tdgal(2,ig,id)
        enddo
      enddo
      read(50,*) flagvch

c *** output the input parameters ***
c orbit parameter
      if(flagini.eq.0) then
        write(6,*) ' merger model from ecc and rperi'
        write(6,*) ' ecc, rperi (100 kpc)= ', ecc,rperi
        if(ecc.ne.1.0d0) then
          rapo=rperi*(1.0d0+ecc)/(1.0d0-ecc)
          write(6,*) ' rapo (100 kpc)= ',rapo
          if(rini.gt.rapo) then
            write(6,*) ' rini must be smaller than rapo'
            stop
          endif
        endif
      else if(flagini.eq.1) then
        write(6,*) ' gal 2 smashing model'
        write(6,*) ' v0x (km/s), b (kpc)=',v0x,bimpz
c converting to GCD+ uinit
        v0x=v0x/VUKMS
        bimpz=bimpz/LUKPC
      else if(flagini.eq.2) then
        write(6,*) ' gal 2 initial position and velocity'
        write(6,*) ' x,z (kpc) =',x0,z0
        write(6,*) ' vx, vy (km/s) =',vx0,vy0
        x0=x0/LUKPC
        z0=z0/LUKPC
        vx0=vx0/VUKMS
        vy0=vy0/VUKMS
      endif
      write(6,*) ' rini (100 kpc) = ',rini
      do ig=1,ngal
        write(6,*) ' for galaxy',ig
        write(6,*) ' use disk+bulge+NFW halo model'
        write(6,*) '  read disk.',ig,'.dat bulge',ig
     &   ,'.dat halo',ig,'.dat '
        write(6,*) ' ndisk,nbul,nhalo,nhalogas ='
     &   ,ndisk(ig),nbul(ig),nhalo(ig),nhgas(ig)
        if(flagrot(ig).eq.0) then
          write(6,*) ' angles i,omega (degree) =',igal(ig),omgal(ig)
          igal(ig) = M_PI*igal(ig)/180.0d0
          omgal(ig) = M_PI*omgal(ig)/180.0d0
        else if(flagrot(ig).lt.0) then
          write(6,*) ' rotation axis is negative z axis'
        else
          write(6,*) ' rotation axis is positive z axis'
        endif
        if(flagbul(ig).eq.0) then
          write(6,*) ' bulge in dm.dat'
        else
          write(6,*) ' bulge in star.dat'
          write(6,*) '  age range=',tbgal(1,ig),tbgal(2,ig)
          write(6,*) ' [Fe/H],[al/Fe] bulge=',metbgal(ig),alfeb(ig)
        endif
        if(nhgas(ig).ge.1) then
          write(6,*) ' halo hot gas [Fe/H], [al/Fe]='
     &      ,fehhg(ig),alfehg(ig)
        endif
        do id=1,ndisk(ig)
          if(flagd(ig,id).gt.0) then
            write(6,*) id,' gas disk '
c            write(6,*) ' [Fe/H] slope,0,scatter =',metgal(1,ig,id)
c     &       ,metgal(2,ig,id),metgal(3,ig,id)
c            write(6,*) ' [al/Fe] slope,0,scatter =',alfe(1,ig,id)
c     &       ,alfe(2,ig,id),alfe(3,ig,id)
          else if(flagd(ig,id).lt.0) then
            write(6,*) id,' star disk age range=',tdgal(1,ig,id)
     &       ,tdgal(2,ig,id)
            tagegal(1,ig,id)=tdgal(1,ig,id)
            tagegal(2,ig,id)=tdgal(2,ig,id)
            write(6,*) ' [Fe/H] slope,0,scatter =',metgal(1,ig,id)
     &       ,metgal(2,ig,id),metgal(3,ig,id)
            write(6,*) ' slope of age metallicity =',metgal(4,ig,id)
          endif
        enddo
        if(flagvch.ne.0.and.nhalo(ig).gt.0) then
          write(6,*) ' velocity will be corrected in halo.'
        endif
      enddo

c *** read mass group info ***
      llz_smgr=-6.0d0
      luz_smgr=0.0d0
      open(50,file='ini/mgr.dat',status='old') 
      read(50,'(a6,I10)') ctmp,nmgr
      read(50,'(a6,I10)') ctmp,nzmgr
      write(6,*) ' nmgr,nzmgr=',nmgr,nzmgr
      if(nmgr.gt.MNMGR) then
        write(6,*) ' Error: increase MNMGR=',MNMGR
        stop
      else if(nzmgr.gt.MNZMGR) then
        write(6,*) ' Error: increase MNZMGR=',MNZMGR
        stop 
      endif
      dlz_smgr=(luz_smgr-llz_smgr)/dble(nzmgr-2)
      do j=0,nzmgr-1
c *** set metallicity ***
        if(j.eq.0) then
          zmet_smgr(j)=0.0d0
        else       
          zmet_smgr(j)=10.0d0**(llz_smgr+dlz_smgr*dble(j-1))
        endif
        do i=1,nmgr
          read(50,'(a200)') cread
          if(cread(1:1).eq.'#') then
            nejmgr(j)=i-1
            goto 90
c 155      forma((I8,(1pE13.5))
          else
            read(cread(9:21),'(1pE13.5)') te_smgr(i,j)
          endif
        enddo
 90   enddo
      read(50,'(a9,I10)') ctmp,nsniimgr
      write(6,*) ' nsniimgr=',nsniimgr
      do j=0,nzmgr-1
c *** set end time for SNe II stars
        do i=1,nsniimgr-1
          te_smgr(i,j)=te_smgr(nsniimgr,j)          
        enddo
c *** set end time for the rest of the mass group
        do i=nejmgr(j)+1,nmgr
          te_smgr(i,j)=2.0d0*te_smgr(nejmgr(j),j)
        enddo
      enddo
      open(60,file='testmgr.dat',status='unknown')
      do j=0,nzmgr-1
        do i=1,nmgr
          write(60,'(I8,2(1pE13.5))') i,te_smgr(i,j),zmet_smgr(j)
        enddo
      enddo
      close(60)

c /*** initialization ***/
      np=0
      ndm=0
      ns=0

c *** counting number of particles ***

      do ig=1,ngal
        ndmg(ig) = 0
        npg(ig) = 0
        mhalo = 0.0d0
        mbul = 0.0d0
        mdisk = 0.0d0
        mhgas=0.0d0
        mgalt(ig) = 0.0d0
        mwnh(ig) = 0
        mwnb(ig) = 0
        do id=1,ndisk(ig)
          tmassd(id)=0.0d0
          cmxd(id)=0.0d0
          cmyd(id)=0.0d0
          cmzd(id)=0.0d0
          vmxd(id)=0.0d0
          vmyd(id)=0.0d0
          vmzd(id)=0.0d0
          mwnd(ig,id) = 0
        enddo
        mwnhg(ig) = 0
c /*** Open Initial Data File ***/
c /*** disk ***/
        if(ndisk(ig).gt.0) then
          do id=1,ndisk(ig) 
            write(filen,'(a8,i1,a1,i1,a4)') 
     &        'ini/disk',ig,'-',id,'.dat'
            open(51,file=filen,status='old',form='unformatted')
            read(51) mwnd(ig,id),flagid(ig,id)
            if(flagid(ig,id).ne.0) then
              read(51) metgal(1,ig,id),metgal(2,ig,id)
     &          ,alfe(1,ig,id),alfe(2,ig,id)
              metgal(3,ig,id)=0.0d0
              metgal(4,ig,id)=0.0d0
              alfe(3,ig,id)=0.0d0
            endif
            if(flagid(ig,id).ne.0) then
              write(6,*) ' gas disk nd=',mwnd(ig,id)
              write(6,*) ' [Fe/H] slope,0,scatter =',metgal(1,ig,id)
     &          ,metgal(2,ig,id),metgal(3,ig,id)
              write(6,*) ' [al/Fe] slope,0,scatter =',alfe(1,ig,id)
     &         ,alfe(2,ig,id),alfe(3,ig,id)
              if(flagd(ig,id).le.0) then
                write(6,*) ' ***** Warning: gas disk date is read as'
                write(6,*) ' ***** flag=',flagd(ig,id)
              endif
            endif
c /*** Read and  Set Position ***/
            if(flagd(ig,id).eq.0) then
              ndm=ndm+mwnd(ig,id)
              ndmg(ig)=ndmg(ig)+mwnd(ig,id)
            else if(flagd(ig,id).gt.0) then
              if(flagid(ig,id).eq.0) then
                write(6,*) ' Error: disk,',ig,id,' flag setted'
     &           ,' to read gas data'
                write(6,*) ' but it is not gas data flag='
     &           ,flagid(ig,id)
                stop
              endif
              np=np+mwnd(ig,id)
              npg(ig)=npg(ig)+mwnd(ig,id)
            else
              ns=ns+mwnd(ig,id)
              nsg(ig)=nsg(ig)+mwnd(ig,id)
            endif

c /*** Read and  Set Position ***/
            if(flagd(ig,id).eq.1) then
c *** gas disk ***
              do i=1,mwnd(ig,id)
                read(51) xp,yp,zp,vxp,vyp,vzp,mpart,rhop,pg,myup
                mdisk = mdisk+mpart
                tmassd(id)=tmassd(id)+mpart
                cmxd(id)=cmxd(id)+mpart*xp
                cmyd(id)=cmyd(id)+mpart*yp
                cmzd(id)=cmzd(id)+mpart*zp
                vmxd(id)=vmxd(id)+mpart*vxp
                vmyd(id)=vmyd(id)+mpart*vyp
                vmzd(id)=vmzd(id)+mpart*vzp
              enddo
            else
c *** DM or stellar disk
              do i=1,mwnd(ig,id)
                read(51) xp,yp,zp,
     &           vxp,vyp,vzp,mpart,rhop
                mdisk = mdisk+mpart
                tmassd(id)=tmassd(id)+mpart
                cmxd(id)=cmxd(id)+mpart*xp
                cmyd(id)=cmyd(id)+mpart*yp
                cmzd(id)=cmzd(id)+mpart*zp
                vmxd(id)=vmxd(id)+mpart*vxp
                vmyd(id)=vmyd(id)+mpart*vyp
                vmzd(id)=vmzd(id)+mpart*vzp
              enddo
            endif
            close(51)
            if(tmassd(id).gt.0.0d0) then
              cmxd(id)=cmxd(id)/tmassd(id)
              cmyd(id)=cmyd(id)/tmassd(id)
              cmzd(id)=cmzd(id)/tmassd(id)
              vmxd(id)=vmxd(id)/tmassd(id)
              vmyd(id)=vmyd(id)/tmassd(id)
              vmzd(id)=vmzd(id)/tmassd(id)
            endif
          enddo
        endif
        write(6,*) ' after disk ng,ndm,ns=',np,ndm,ns
c /*** bulge ***/  
        if(nbul(ig).gt.0) then
          write(filen,'(a9,i1,a4)') 'ini/bulge',ig,'.dat'
          open(51,file=filen,status='old',form='unformatted')
          read(51) mwnb(ig)
          if(flagbul(ig).eq.0) then
            ndm = ndm+mwnb(ig)
            ndmg(ig) = ndmg(ig)+mwnb(ig)
          else
c *** stars ***
            ns=ns+mwnb(ig)
            nsg(ig)=nsg(ig)+mwnb(ig)
          endif
          do i=1,mwnb(ig)
c /*** Read and  Set Position ***/
            read(51) xp,yp,zp,vxp,vyp,vzp,mpart,rhop
            mbul = mbul+mpart
          enddo
          close(51)
        endif 
        write(6,*) ' after bulge ng,ndm,ns=',np,ndm,ns

c *** halo hot gas ***
        if(nhgas(ig).gt.0) then
          write(filen,'(a8,i1,a4)') 'ini/hgas,',ig,'.dat'
          open(51,file=filen,status='old',form='unformatted')
          read(51) mwnhg(ig)
          np = np+mwnhg(ig)
          do i=1,mwnhg(ig)
            read(51) xp,yp,zp,vxp,vyp,vzp,mpart,rhop,uthp
            mhgas=mhgas+mpart
          enddo
          npg(ig) = npg(ig)+mwnhg(ig)
          close(51)
        endif
        write(6,*) ' after halo gas ng,ndm,ns=',np,ndm,ns
c /*** Halo ***/  
        if(nhalo(ig).gt.0) then
          write(filen,'(a8,i1,a4)') 'ini/halo',ig,'.dat'
          open(51,file=filen,status='old',form='unformatted')
          read(51) mwnh(ig)
          ndm = ndm+mwnh(ig)
          ndmg(ig)=ndmg(ig)+mwnh(ig)
          tmassh(ig)=0.0d0
          cmxh(ig)=0.0d0
          cmyh(ig)=0.0d0
          cmzh(ig)=0.0d0
          vmxh(ig)=0.0d0
          vmyh(ig)=0.0d0
          vmzh(ig)=0.0d0
          do i=1,mwnh(ig)
            read(51) xp,yp,zp,vxp,vyp,vzp,mpart,rhop
            tmassh(ig)=tmassh(ig)+mpart
            cmxh(ig)=cmxh(ig)+xp*mpart
            cmyh(ig)=cmyh(ig)+yp*mpart
            cmzh(ig)=cmzh(ig)+zp*mpart
            vmxh(ig)=vmxh(ig)+vxp*mpart
            vmyh(ig)=vmyh(ig)+vyp*mpart
            vmzh(ig)=vmzh(ig)+vzp*mpart
            mhalo = mhalo+mpart
          enddo
          if(tmassh(ig).gt.0.0d0) then
            cmxh(ig)=cmxh(ig)/tmassh(ig)
            cmyh(ig)=cmyh(ig)/tmassh(ig)
            cmzh(ig)=cmzh(ig)/tmassh(ig)
            vmxh(ig)=vmxh(ig)/tmassh(ig)
            vmyh(ig)=vmyh(ig)/tmassh(ig)
            vmzh(ig)=vmzh(ig)/tmassh(ig)
          endif
        endif
        mgalt(ig) = mhalo+mbul+mdisk+mhgas

        write(6,*) ' after DM halo ng,ndm,ns=',np,ndm,ns
        write(6,*) '*** Galaxy',ig,' ***'
        write(6,*) ' number of particles and mass'
        write(6,*) '  halo =',mwnh(ig),mhalo
        write(6,*) ' bulge =',mwnb(ig),mbul
        write(6,*) ' halo gas=',mwnhg(ig),mhgas
        do id=1,ndisk(ig)
          write(6,*) ' disk',id,'=',mwnd(ig,id)
          write(6,*) ' cemtre of the mass=',cmxd(id),cmyd(id),cmzd(id)
          write(6,*) ' momentum=',vmxd(id),vmyd(id),vmzd(id)
        enddo
        write(6,*) ' disk mass total=',mdisk
        write(6,*) ' total galaxy mass=',mgalt(ig)
        if(nhalo(ig).gt.0.0d0) then
          write(6,*) ' Halo centre of the mass ='
     &       ,cmxh(ig),cmyh(ig),cmzh(ig)
          write(6,*) '  momentum =',vmxh(ig),vmyh(ig),vmzh(ig)
        endif
      enddo
      ng = np
      write(6,*) ' total number of gas, DM and star particles = '
     &  ,np,ndm,ns
      mhalo = 0.0d0
      mbul = 0.0d0
      mdisk = 0.0d0
      mhgas=0.0d0

c /*** calculate orbital parameter ***/

      do i=1,ngal
        xini(i)=0.0d0
        zini(i)=0.0d0
        vxini(i)=0.0d0
        vyini(i)=0.0d0
      enddo
      if(ngal.gt.1) then
        if(flagini.eq.0) then

c *** the reduced mass ***/
          myum = mgalt(1)*mgalt(2)/(mgalt(1)+mgalt(2))
c *** angular momentum ***
c      if(ecc.ne.1.0) then
c        ltot = dsqrt(rperi*(1.0-ecc**2)*myum*G*mgalt(1)*mgalt(2)
c     &    /(1.0d0-ecc))
          ltot = dsqrt(rperi*(1.0+ecc)*myum*G*mgalt(1)*mgalt(2))
c
c        ltot = 0.0d0
c      endif
c      ltot = dsqrt((ecc+1.0d0)*myum*G*mgalt(1)*mgalt(2))
c *** total energy ***
          if(ltot.ne.0.0d0.and.ecc.ne.1.0d0) then
            etot = (ecc**2-1.0d0)*(myum*(G*mgalt(1)*mgalt(2))**2)
     &      /(2.0d0*(ltot**2))
          else
            etot = 0.0d0
          endif
          write(6,*) ' ltot, etot = ',ltot,etot
c *** set initial position ***
          xini(1)=-rini*mgalt(2)/(mgalt(1)+mgalt(2))
          xini(2)=rini*mgalt(1)/(mgalt(1)+mgalt(2))
          write(6,*) ' rini = ',rini
          write(6,*) ' initial x-position for Gal 1 and 2=',xini(1)
     &     ,xini(2)
c *** y velocities ***
          vy12=ltot/(myum*rini)
          vyini(1)=mgalt(2)*vy12/(mgalt(1)+mgalt(2))
          vyini(2)=-mgalt(1)*vy12/(mgalt(1)+mgalt(2))
c ***  relative x-velocity for ***
          vx12=dsqrt((etot+G*mgalt(1)*mgalt(2)/rini
     &     -0.5d0*myum*(vy12**2))*2.0d0/myum)
c since 30/11/2010
c        vx12=dsqrt((etot+G*mgalt(1)*mgalt(2)/rini)*2.0d0/myum)

c *** set initial velocity ***
          vxini(1)=vx12*mgalt(2)/(mgalt(1)+mgalt(2))
          vxini(2)=-vx12*mgalt(1)/(mgalt(1)+mgalt(2))
          write(6,*) ' etot -Ekin,y=',etot+G*mgalt(1)*mgalt(2)/rini
     &     -0.5d0*myum*(vy12**2)
          write(6,*) ' vini = ',vx12,vy12
          write(6,*) ' initial velocity for Gal 1 =',vxini(1),vyini(1)
          write(6,*) ' initial velocity for Gal 2 =',vxini(2),vyini(2)
        else if(flagini.eq.1) then
          vxini(2)=v0x
          zini(2)=bimpz
          xini(2)=-rini
        else if(flagini.eq.2) then
          xini(2)=x0
          zini(2)=z0
          vxini(2)=vx0
          vyini(2)=vy0
        endif
      else
        xini(1)=0.0d0
        xini(2)=0.0d0
        vxini(1)=0.0d0
        vxini(2)=0.0d0
        vyini(1)=0.0d0
        vyini(2)=0.0d0
      endif

c *** open output files ***
c *** binary files
      itmp=1
      open(60,file='gas.dat',status='unknown'
     &  ,form='unformatted')
      write(60) itmp,itmp,ng,ng,1
      open(61,file='gas-metal.dat',status='unknown'
     &  ,form='unformatted')
      open(62,file='star.dat',status='unknown'
     &  ,form='unformatted')
      write(62) itmp,itmp,ns,ns
      open(63,file='star-metal.dat',status='unknown'
     &  ,form='unformatted')
      open(64,file='dm.dat',status='unknown'
     &  ,form='unformatted')
      write(64) itmp,itmp,ndm,ndm
c *** ascii data
      open(70,file='agas.dat',status='unknown')
      open(71,file='astar.dat',status='unknown')
      open(72,file='adm.dat',status='unknown')
      open(73,file='admhalo.dat',status='unknown')
      mcg=0.0d0
      do ig=1,ngal
c /*** Open Initial Data File ***/
c /*** disk ***/
        if(ndisk(ig).gt.0) then
          do id=1,ndisk(ig) 
            write(filen,'(a8,i1,a1,i1,a4)') 
     &        'ini/disk',ig,'-',id,'.dat'
            open(51,file=filen,status='old',form='unformatted')
            read(51) mwnd(ig,id),flagid(ig,id)
            if(flagid(ig,id).ne.0) then
              read(51) metgal(1,ig,id),metgal(2,ig,id)
     &          ,alfe(1,ig,id),alfe(2,ig,id)
              metgal(3,ig,id)=0.0d0
              metgal(4,ig,id)=0.0d0
              alfe(3,ig,id)=0.0d0
            endif
c /*** Read and  Set Position ***/
            if(flagd(ig,id).eq.0) then
c *** DM disk
              do i=1,mwnd(ig,id)
                read(51) xp,yp,zp,
     &           vxp,vyp,vzp,mpart,rhop
                mdisk = mdisk+mpart
c *** inclination setting ***
                call rot(xp,yp,zp,vxp,vyp,vzp,ig,igal,omgal,flagrot(ig))
c *** initial positon and velocity shift
                xp=xp+xini(ig)            
                zp=zp+zini(ig)
                vxp=vxp+vxini(ig)
                vyp=vyp+vyini(ig)
c *** ascii output
               if(mod(i,npaskip).eq.0) then
                 write(72,162) xp,yp,zp,vxp,vyp,vzp,mpart,rhop
 162  format(8(1pE13.5))
               endif    
c *** binary output        
               write(64) xp,yp,zp,vxp,vyp,vzp,mpart,rhop
             enddo
           else if(flagd(ig,id).gt.0) then
c *** gas disk ***
              do i=1,mwnd(ig,id)
                read(51) xp,yp,zp,vxp,vyp,vzp,mpart,rhop,pg,myup
                mdisk = mdisk+mpart
c set Galactic radius before changing the inclination and the centre
                rg=dsqrt(xp**2+yp**2)
c *** inclination setting ***
                call rot(xp,yp,zp,vxp,vyp,vzp,ig,igal,omgal,flagrot(ig)
     &           ,xini,vxini,vyini)
c *** initial positon and velocity shift
                xp=xp+xini(ig)            
                zp=zp+zini(ig)            
                vxp=vxp+vxini(ig)
                vyp=vyp+vyini(ig)
c *** gas properties ***
                ug=pg/((GAM-1.0d0)*rhop)
                hg=(mpart/rhop)**THIRD
                flagfdg=0
c *** gas disk metallicity, use 2D radius ***
                fehp=metgal(1,ig,id)*(rg*LUKPC)+metgal(2,ig,id)
     &            +metgal(3,ig,id)*dble(gasdev(idum))
                fehg=fehp
c *** [alpha/Fe] as a function of [Fe/H]
                alfep=alfe(1,ig,id)*fehp+alfe(2,ig,id)
     &            +alfe(3,ig,id)*dble(gasdev(idum))
                metp=10.0d0**fehp
c *** get H mass from scaling
                mzHg=MUSM*mpart
     &         *(0.76d0+(xsol(1)-0.76d0)*metp)
c *** total metal abundance 
                if(fehp.lt.fehp+alfep) then
                  metp=10.0d0**(fehp+alfep)
                endif
c for C and/or N enhanced
                if(fehp.lt.fehp-alfep) then
                  metp=10.0d0**(fehp-alfep)
                endif
c *** Z ***
c                mzZg(i)=musm*mg(i)*XSZ*metp
                mzZg=metp*(XSZ/xsol(1))*mzHg
                mzHeg=MUSM*mpart-mzHg-mzZg
                if(mzHeg.lt.0.0d0) then
                  write(6,*) ' He<0 for gas=',i,ig,id
     &             ,fehp,alfep
                endif
c /* 1  2   3   4   5   6    7    8    9    10   11  12   13   14   15
c /* 1H,4He,12C,14N,16O,20Ne,23Na,24Mg,27Al,28Si,32S,36Ar.40Ca,56Fe,58Ni */
c *** C/H and N/H ***
                metpc=10.0d0**(fehp-alfep)
                MzCg=metpc*mzHg*(xsol(3)/xsol(1))
                mzNg=metpc*mzHg*(xsol(4)/xsol(1))
c *** O/H and etc
                metpo=10.0d0**(fehp+alfep)
                mzOg=metpo*mzHg*(xsol(5)/xsol(1))
                mzNeg=metpo*mzHg*(xsol(6)/xsol(1))
                mzMgg=metpo*mzHg*(xsol(8)/xsol(1))
                mzSig=metpo*mzHg*(xsol(10)/xsol(1))
c *** Fe/H
                metpfe=10.0d0**(fehp)
                mzFeg=metpfe*mzHg*(xsol(14)/xsol(1))
c *** check ***
                if(mzZg.lt.mzCg+mzNg+mzOg+mzNeg
     &            +mzMgg+mzSig+mzFeg) then
                  write(6,*) ' Error: alfe is too big ig,id=',ig,id
                  write(6,*) ' metp,feh,alfe=',metp,fehp,alfep
                  write(6,*) ' MZ, MZsum=',mzZg
     &             ,mzCg+mzNg+mzOg+mzNeg
     &             +mzMgg+mzSig+mzFeg
                  stop
                endif
c *** Nh density
                nhp=(((mpart-(mzZg+mzHeg)/MUSM))/mpart)
     &          *rhop*(DU/MP)
                if(nhp.gt.0.1d0) then
                  mcg=mcg+mpart
                endif
c *** ascii output
                if(mod(i,npaskip).eq.0) then
                  vrotp=(vxp*yp-vyp*xp)/dsqrt(xp**2+yp**2)
                  write(70,160) xp,yp,zp,vxp,vyp,vzp
     &             ,mpart,rhop,ug,mzHeg,mzCg,mzNg
     &             ,mzOg,mzNeg,mzMgg,mzSig,mzFeg,mzZg
     &             ,ug*(GAM-1.0d0)*TUK*myup/(MYU*TPRHO)
     &             ,mzZg/(mpart*MUSM*XSZ)
     &             ,dlog10((mzFeg/mzHg)/(xsol(14)/xsol(1)))
     &             ,dlog10((mzOg/mzFeg)/(xsol(5)/xsol(14)))
     &             ,fehg,rg,vrotp,nhp
 160  format(26(1pE13.5))
                endif
c *** binary output ***
                write(60) xp,yp,zp,vxp,vyp,vzp,mpart
     &           ,rhop,ug,0.0d0,flagfdg
                write(61) mzHeg,mzCg,mzNg,mzOg,mzNeg,mzMgg,mzSig
     &           ,mzFeg,mzZg
              enddo
              write(6,*) ' M(nh>0.1)=',mcg
            else
c *** stars ***
              do i=1,mwnd(ig,id)
                read(51) xp,yp,zp,vxp,vyp,vzp,mpart,rhop
                mdisk = mdisk+mpart
c set radius before inclination and the centre changed. 
                rs=dsqrt(xp**2+yp**2)
c *** inclination setting ***
                call rot(xp,yp,zp,vxp,vyp,vzp,ig,igal,omgal,flagrot(ig)
     &           ,xini,vxini,vyini)
c *** initial positon and velocity shift
                xp=xp+xini(ig)            
                zp=zp+zini(ig)            
                vxp=vxp+vxini(ig)
                vyp=vyp+vyini(ig)
c *** set age ***
 91             ages=(tagegal(1,ig,id)+
     &            dble(ran1(idum))
     &            *(tagegal(2,ig,id)-tagegal(1,ig,id)))/TMUGYR
c *** stellar disk metallicity ***
                fehp=metgal(1,ig,id)*(rs*LUKPC)+metgal(2,ig,id)
     &            +ages*metgal(4,ig,id)*TMUGYR
     &            +metgal(3,ig,id)*dble(gasdev(idum))
                fehs=fehp
c *** [alpha/Fe] as a function of [Fe/H]
                alfep=alfe(1,ig,id)*fehp+alfe(2,ig,id)
     &            +alfe(3,ig,id)*dble(gasdev(idum))
c *** get H mass from scaling
                metp=10.0d0**fehp
                mzHs=MUSM*mpart*(0.76d0+(xsol(1)-0.76d0)*metp)
c *** total metal abundance 
                if(fehp.lt.fehp+alfep) then
                  metp=10.0d0**(fehp+alfep)
                endif
c for C and/or N enhanced
                if(fehp.lt.fehp-alfep) then
                  metp=10.0d0**(fehp-alfep)
                endif
                zmetp=metp
c *** check end time for group ***               
                if(zmetp.gt.0.0d0) then
                  iz1=int((dlog10(zmetp)-llz_smgr)/dlz_smgr)+1
                  if(iz1.lt.0) then
                    iz1=0
                  else if(iz1.ge.nzmgr) then
                    iz1=nzmgr-1
                  endif
                  wz2=(zmetp-zmet_smgr(iz1))
     &             /(zmet_smgr(iz1+1)-zmet_smgr(iz1))
                  wz1=1.0d0-wz2
                  do im=1,nmgr
                    tepmgr(im)=wz1*te_smgr(im,iz1)
     &                +wz2*te_smgr(im,iz1+1)
                  enddo
                else
                  do im=1,nmgr
                    tepmgr(im)=te_smgr(im,0)
                  enddo
                endif
c                flagfds=int(ran1(idum)*dble(nspgr-1))+2
                flagfds=int(ran1(idum)*dble(nmgr))+1
                if(flagfds.gt.nmgr) then
                   flagfds=nmgr
                endif
                if(ages.gt.tepmgr(flagfds)) then
                  goto 91
                endif
c *** Z ***
c                mzZs=musm*mpart*XSZ*metp
                mzZs=(XSZ/xsol(1))*mzHs*metp
                mzHes=MUSM*mpart-mzHs-mzZs
                if(mzHes.lt.0.0d0) then
                  write(6,*) ' He<0 for star=',i,ig,id
     &             ,fehp,alfep
                endif
c /* 1  2   3   4   5   6    7    8    9    10   11  12   13   14   15
c /* 1H,4He,12C,14N,16O,20Ne,23Na,24Mg,27Al,28Si,32S,36Ar.40Ca,56Fe,58Ni */
c *** C/H and N/H ***
                metpc=10.0d0**(fehp-alfep)
                MzCs=metpc*mzHs*(xsol(3)/xsol(1))
                mzNs=metpc*mzHs*(xsol(4)/xsol(1))
c *** O/H and etc
                metpo=10.0d0**(fehp+alfep)
                mzOs=metpo*mzHs*(xsol(5)/xsol(1))
                mzNes=metpo*mzHs*(xsol(6)/xsol(1))
                mzMgs=metpo*mzHs*(xsol(8)/xsol(1))
                mzSis=metpo*mzHs*(xsol(10)/xsol(1))
c *** Fe/H
                metpfe=10.0d0**(fehp)
                mzFes=metpfe*mzHs*(xsol(14)/xsol(1))
c *** check ***
                if(mzZs.lt.mzCs+mzNs+mzOs+mzNes
     &            +mzMgs+mzSis+mzFes) then
                  write(6,*) ' Error: alfe is too big star ig,id='
     &             ,ig,id,alfep,fehp
                  write(6,*) ' MZ, MZsum=',mzZs
     &             ,mzCs+mzNs+mzOs+mzNes
     &             +mzMgs+mzSis+mzFes
                  stop
                endif
c *** ascii output
                if(mod(i,npaskip).eq.0) then
                  vrotp=(vxp*yp-vyp*xp)/dsqrt(xp**2+yp**2)
                write(71,161) xp,yp,zp,vxp,vyp,vzp,mpart,mpart,rhop
     &            ,ages,mzHes,mzCs,mzNs
     &            ,mzOs,mzNes,mzMgs,mzSis,mzFes,mzZs
     &            ,mzZs/(mpart*MUSM*XSZ)
     &            ,dlog10((mzFes/mzHs)/(xsol(14)/xsol(1)))
     &            ,dlog10((mzOs/mzFes)/(xsol(5)/xsol(14)))
     &            ,fehs,rs,vrotp,flagfds
                endif
 161  format(25(1pE13.5),I5)
c *** binary output ***
                write(62) xp,yp,zp,vxp,vyp,vzp,mpart,mpart,rhop
     &           ,ages,flagfds
                write(63) mzHes,mzCs,mzNs,mzOs,mzNes,mzMgs
     &           ,mzSis,mzFes,mzZs
              enddo
            endif
          enddo
        endif
c *** end of the disk data

c *** bulge ***
        if(nbul(ig).gt.0) then
          write(filen,'(a9,i1,a4)') 'ini/bulge',ig,'.dat'
          open(51,file=filen,status='old',form='unformatted')
          read(51) mwnb(ig)

          if(flagbul(ig).eq.0) then
            do i=1,mwnb(ig)
c /*** Read and  Set Position ***/
              read(51) xp,yp,zp,vxp,vyp,vzp,mpart,rhop
              mbul = mbul+mpart
c *** initial positon and velocity shift
              xp=xp+xini(ig)            
              zp=zp+zini(ig)            
              vxp=vxp+vxini(ig)
              vyp=vyp+vyini(ig)
c *** DM bulge
c *** ascii output
              if(mod(i,npaskip).eq.0) then
                write(72,162) xp,yp,zp,vxp,vyp,vzp,mpart,rhop
              endif    
c *** binary output        
              write(64) xp,yp,zp,vxp,vyp,vzp,mpart,rhop
            enddo
          else
c *** stellar bulge ***
            do i=1,mwnb(ig)
c /*** Read and  Set Position ***/
              read(51) xp,yp,zp,vxp,vyp,vzp,mpart,rhop
              mbul = mbul+mpart
c *** initial positon and velocity shift
              xp=xp+xini(ig)            
              zp=zp+zini(ig)            
              vxp=vxp+vxini(ig)
              vyp=vyp+vyini(ig)
c *** set age ***
              ages=(tbgal(1,ig)+
     &            dble(ran1(idum))
     &            *(tbgal(2,ig)-tbgal(1,ig)))/TMUGYR
              fehp=metbgal(ig)
              fehs=fehp
c *** [alpha/Fe] as a function of [Fe/H]
              alfep=alfeb(ig)
c *** get H mass from scaling
              metp=10.0d0**fehp
              mzHs=MUSM*mpart*(0.76d0+(xsol(1)-0.76d0)*metp)
c *** total metal abundance 
              if(fehp.lt.fehp+alfep) then
                metp=10.0d0**(fehp+alfep)
              endif
c for C and/or N enhanced
              if(fehp.lt.fehp-alfep) then
                metp=10.0d0**(fehp-alfep)
              endif
              zmetp=metp
c *** check end time for group ***               
              if(zmetp.gt.0.0d0) then
                iz1=int((dlog10(zmetp)-llz_smgr)/dlz_smgr)+1
                if(iz1.lt.0) then
                  iz1=0
                else if(iz1.ge.nzmgr) then
                  iz1=nzmgr-1
                endif
                wz2=(zmetp-zmet_smgr(iz1))
     &           /(zmet_smgr(iz1+1)-zmet_smgr(iz1))
                wz1=1.0d0-wz2
                do im=1,nmgr
                  tepmgr(im)=wz1*te_smgr(im,iz1)
     &              +wz2*te_smgr(im,iz1+1)
                enddo
              else
                do im=1,nmgr
                  tepmgr(im)=te_smgr(im,0)
                enddo
              endif
c                flagfds=int(ran1(idum)*dble(nspgr-1))+2
 92           flagfds=int(ran1(idum)*dble(nmgr))+1
              if(flagfds.gt.nmgr) then
                 flagfds=nmgr
              endif
              if(ages.gt.tepmgr(flagfds)) then
                goto 92
              endif
c *** Z ***
c                mzZs=musm*mpart*XSZ*metp
              mzZs=(XSZ/xsol(1))*mzHs*metp
              mzHes=MUSM*mpart-mzHs-mzZs
              if(mzHes.lt.0.0d0) then
                write(6,*) ' He<0 for star=',i,ig,id
     &           ,fehp,alfep
              endif
c /* 1  2   3   4   5   6    7    8    9    10   11  12   13   14   15
c /* 1H,4He,12C,14N,16O,20Ne,23Na,24Mg,27Al,28Si,32S,36Ar.40Ca,56Fe,58Ni */
c *** C/H and N/H ***
              metpc=10.0d0**(fehp-alfep)
              MzCs=metpc*mzHs*(xsol(3)/xsol(1))
              mzNs=metpc*mzHs*(xsol(4)/xsol(1))
c *** O/H and etc
              metpo=10.0d0**(fehp+alfep)
              mzOs=metpo*mzHs*(xsol(5)/xsol(1))
              mzNes=metpo*mzHs*(xsol(6)/xsol(1))
              mzMgs=metpo*mzHs*(xsol(8)/xsol(1))
              mzSis=metpo*mzHs*(xsol(10)/xsol(1))
c *** Fe/H
              metpfe=10.0d0**(fehp)
              mzFes=metpfe*mzHs*(xsol(14)/xsol(1))
c *** check ***
              if(mzZs.lt.mzCs+mzNs+mzOs+mzNes
     &          +mzMgs+mzSis+mzFes) then
                write(6,*) ' Error: alfe is too big star ig,id='
     &           ,ig,id,alfep,fehp
                write(6,*) ' MZ, MZsum=',mzZs
     &           ,mzCs+mzNs+mzOs+mzNes
     &           +mzMgs+mzSis+mzFes
                stop
              endif
c *** ascii output
              if(mod(i,npaskip).eq.0) then
                vrotp=(vxp*yp-vyp*xp)/dsqrt(xp**2+yp**2)
                write(71,161) xp,yp,zp,vxp,vyp,vzp,mpart,mpart,rhop
     &           ,ages,mzHes,mzCs,mzNs
     &           ,mzOs,mzNes,mzMgs,mzSis,mzFes,mzZs
     &           ,mzZs/(mpart*MUSM*XSZ)
     &           ,dlog10((mzFes/mzHs)/(xsol(14)/xsol(1)))
     &           ,dlog10((mzOs/mzFes)/(xsol(5)/xsol(14)))
     &           ,fehs,rs,vrotp,flagfds
              endif
c *** binary output ***
              write(62) xp,yp,zp,vxp,vyp,vzp,mpart,mpart,rhop
     &          ,ages,flagfds
              write(63) mzHes,mzCs,mzNs,mzOs,mzNes,mzMgs
     &          ,mzSis,mzFes,mzZs
            enddo
          endif
        endif
c *** end of the bulge data

c *** halo hot gas ***
        if(nhgas(ig).gt.0) then
          write(filen,'(a8,i1,a4)') 'ini/hgas,',ig,'.dat'
          open(51,file=filen,status='old',form='unformatted')
          read(51) mwnhg(ig)
          do i=1,mwnhg(ig)
            read(51) xp,yp,zp,vxp,vyp,vzp,mpart,rhop,uthp
            mhgas=mhgas+mpart
c set the radius before inclination and the centre changed
            r2dg=dsqrt(xp**2+yp**2)
            rg=dsqrt(xp**2+yp**2+zp**2)
c *** inclination setting ***
            call rot(xp,yp,zp,vxp,vyp,vzp,ig,igal,omgal,flagrot(ig)
     &        ,xini,vxini,vyini)
c *** initial positon and velocity shift
            xp=xp+xini(ig)            
            zp=zp+zini(ig)            
            vxp=vxp+vxini(ig)
            vyp=vyp+vyini(ig)
c *** gas properties ***
            ug=uthp
            pg=(GAM-1.0d0)*rhop*ug
            hg=(mpart/rhop)**THIRD
            flagfdg=0
c *** hot gas metallicity constant ***
            fehp=fehhg(ig)
            fehg=fehp
c *** [alpha/Fe] constant
            alfep=alfehg(ig)
            metp=10.0d0**fehp
c *** get H mass from scaling
            mzHg=MUSM*mpart
     &        *(0.76d0+(xsol(1)-0.76d0)*metp)
c *** total metal abundance 
            if(fehp.lt.fehp+alfep) then
              metp=10.0d0**(fehp+alfep)
            endif
c for C and/or N enhanced
            if(fehp.lt.fehp-alfep) then
              metp=10.0d0**(fehp-alfep)
            endif
c *** Z ***
c                mzZg(i)=musm*mg(i)*XSZ*metp
            mzZg=metp*(XSZ/xsol(1))*mzHg
            mzHeg=MUSM*mpart-mzHg-mzZg
            if(mzHeg.lt.0.0d0) then
              write(6,*) ' He<0 for hot gas=',i,ig,id
     &         ,fehp,alfep
            endif
c /* 1  2   3   4   5   6    7    8    9    10   11  12   13   14   15
c /* 1H,4He,12C,14N,16O,20Ne,23Na,24Mg,27Al,28Si,32S,36Ar.40Ca,56Fe,58Ni */
c *** C/H and N/H ***
            metpc=10.0d0**(fehp-alfep)
            MzCg=metpc*mzHg*(xsol(3)/xsol(1))
            mzNg=metpc*mzHg*(xsol(4)/xsol(1))
c *** O/H and etc
            metpo=10.0d0**(fehp+alfep)
            mzOg=metpo*mzHg*(xsol(5)/xsol(1))
            mzNeg=metpo*mzHg*(xsol(6)/xsol(1))
            mzMgg=metpo*mzHg*(xsol(8)/xsol(1))
            mzSig=metpo*mzHg*(xsol(10)/xsol(1))
c *** Fe/H
            metpfe=10.0d0**(fehp)
            mzFeg=metpfe*mzHg*(xsol(14)/xsol(1))
c *** check ***
            if(mzZg.lt.mzCg+mzNg+mzOg+mzNeg
     &        +mzMgg+mzSig+mzFeg) then
              write(6,*) ' Error: alfe is too big ig,id=',ig,id
              write(6,*) ' metp,feh,alfe=',metp,fehp,alfep
              write(6,*) ' MZ, MZsum=',mzZg
     &         ,mzCg+mzNg+mzOg+mzNeg
     &         +mzMgg+mzSig+mzFeg
              stop
            endif
c *** Nh density
            nhp=(((mpart-(mzZg+mzHeg)/MUSM))/mpart)
     &        *rhop*(DU/MP)
c *** ascii output
            if(mod(i,npaskip).eq.0) then
              vrotp=(vxp*yp-vyp*xp)/dsqrt(xp**2+yp**2)
c              write(6,*) GAM,TPRHO,TUK,ug*((GAM-1.0d0)/TPRHO)*TUK
              write(70,163) xp,yp,zp,vxp,vyp,vzp
     &          ,mpart,rhop,ug,mzHeg,mzCg,mzNg
     &          ,mzOg,mzNeg,mzMgg,mzSig,mzFeg,mzZg
     &          ,ug*(GAM-1.0d0)*TUK/(TPRHO)
     &          ,mzZg/(mpart*MUSM*XSZ)
     &          ,dlog10((mzFeg/mzHg)/(xsol(14)/xsol(1)))
     &          ,dlog10((mzOg/mzFeg)/(xsol(5)/xsol(14)))
     &          ,fehg,rg,r2dg,vrotp,nhp
 163          format(27(1pE13.5))
            endif
c *** binary output ***
            write(60) xp,yp,zp,vxp,vyp,vzp,mpart
     &        ,rhop,ug,0.0d0,flagfdg
            write(61) mzHeg,mzCg,mzNg,mzOg,mzNeg,mzMgg,mzSig
     &        ,mzFeg,mzZg
          enddo
        endif

c *** DM halo ***
        if(nhalo(ig).gt.0) then
          write(filen,'(a8,i1,a4)') 'ini/halo',ig,'.dat'
          open(51,file=filen,status='old',form='unformatted')
          read(51) mwnh(ig)
          do i=1,mwnh(ig)
            read(51) xp,yp,zp,vxp,vyp,vzp,mpart,rhop
            mhalo = mhalo+mpart
            if(flagvch.ne.0) then
c adjust the velocity   
              vxp=vxp-vmxh(ig)
              vyp=vyp-vmyh(ig)
              vzp=vzp-vmzh(ig)
            endif
c *** initial positon and velocity shift
            xp=xp+xini(ig)            
            zp=zp+zini(ig)            
            vxp=vxp+vxini(ig)
            vyp=vyp+vyini(ig)
c *** ascii output
            if(mod(i,npaskip).eq.0) then
              write(73,162) xp,yp,zp,vxp,vyp,vzp,mpart,rhop
            endif    
c *** binary output        
            write(64) xp,yp,zp,vxp,vyp,vzp,mpart,rhop
          enddo
        endif
        mgalt(ig) = mhalo+mbul+mdisk+mhgas
        write(6,*) '*** Galaxy',ig,' ***'
        write(6,*) ' Mass halo, bulge disk halo gas total (1e12 Msun)='
     &    ,mhalo,mbul,mdisk,mhgas,mgalt(ig)
      enddo
c *** close files ***
      close(70)
      close(71)
      close(72)
      close(73)
c
      write(modname,'(a15)') ' diskm buildini'
      write(60) modname
      write(60) gam
      close(60)
      write(61) modname
      close(61)
      write(62) modname
      close(62)
      write(63) modname
      close(63)
      write(64) modname
      close(64)


      stop 
      end
