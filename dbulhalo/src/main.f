c /***************************************************
c   main.f  for bulhalo ver.37.1
c  15 Dec. 2015  written by D.Kawata
c ****************************************************/
c70*******************************************************************

      include 'define.f'
      integer MNJ
      parameter (MNJ=100000)
      double precision EPSJ
c      parameter (EPSJ=1.0e-1)
      parameter (EPSJ=1.5)
      integer i,j,k,ip,id,idj,it,ir,mnit,itd,ihg
      integer np,pn,nbp,pnpi,pnpo
      integer is,ie,nit
      character filen*60
      integer flagb,flagmexth,flag,flaglog
c /*** halo parameter NFW ***/
      double precision mvir,rvir,cnfw,rs,facnfw,tvir,tvir0
      double precision vvir
c /*** bulge parameter hernquist profile ***/
      double precision ah,re
      double precision xcut,refac,xrt,rt,rcutdmh
      double precision rcutbul
c *** halo hot gas parameters ***
      integer nhg
      double precision mtothg,lamjhg,ajzhg,xrhg,rcuthg
      double precision mhgp,mtothgrc
c *** disk exponential profile ***
      integer ndisk,idsds(MNC),ideds(MNC)
      integer idgdisk
      integer flaggd(MNC),flagvcf,flagflzd(MNC)
      double precision hd(MNC),zd(MNC),Rflzd(MNC)
      double precision rdlim(MNC),zdlim(MNC),frdlim(MNC),fzdlim(MNC)
      double precision maxrdlim,maxzdlim
      double precision frq(MNC)
      integer rlistd(MN)
      real r3dr(MN)
      double precision r3dd(MN)
      double precision xcutgd,frgd
      double precision dfeh(MNC),feh0(MNC),dalfe(MNC),alfe0(MNC)
      double precision nhhalo(MNC),fehhalo(MNC)
      double precision tdisk
      integer nr,nth,nz
      double precision eps,dsig
c *** jeans mass instability u threshold
      double precision fpjm,zred
c /*** stellar/DM ratio ***/
      double precision fs
      double precision mbul,bmtot,tmh,tmtot,tmhg
      double precision mdisk(MNC),mtdisk(MNC)
c /*** cosmology ***/ 
      double precision omg0,lam0,hub,zt,omz,xi,iai,omgb,fb
c /*** for velocity dispersion ***/
      integer npg,ngr,ngz
      double precision ri,ro,rr,dr,prr,dpdr
      double precision zi,zo,zz,dz,pzz
c /*** to generate the particles ***/
      integer nh,nb,nd(MNC)
      double precision mdmhp,mbp,mdp
c *** for external force ***
      double precision fxe(0:MN-1),fye(0:MN-1),fze(0:MN-1)
c *** for check ***
      integer flagexit
      double precision xij,yij,zij,rij,pot
      double precision errmax,errm
c *** external function ***
      double precision mrexpdf
      external mrexpdf

c ***** MPI Initialization *****
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      if(myrank.eq.0) then
        write(6,*) ' ver.37 15/12/2015'
        print *, "Process ",myrank, " of ", nprocs, " is alive"
      endif

c *** cosmological parameter ***
      omg0 = 0.24d0
      hub = 0.73d0
      zt = 1.5d0
c mean molecular weight 1.2 is better?
c      myug=0.6d0
      mnit=100

c /*** read input parameters ***/
      open(50,file='ini/input.dat',status='old')
      read(50,*) flagb
      read(50,*) nh
      read(50,*) flagmexth
      read(50,*) mvir
      read(50,*) cnfw
      read(50,*) xcut,xrt
      read(50,*) omg0,hub,zt
      read(50,*) omgb
      read(50,*) nb
      read(50,*) mbul,refac
      read(50,*) nhg
      read(50,*) mtothg,lamjhg
      read(50,*) xrhg
      if(xrhg.gt.xrt*xcut) then
        if(myrank.eq.0) then
          write(6,*) 'Error: xrhg must be smaller than xrt x xcut.'
          write(6,*) ' xrt,xcut,xrhg=',xrt,xcut,xrhg
        endif
        stop
      endif
      read(50,*) ndisk
      do i=1,ndisk
c flaggd: 0 star, otherwise gas
c <0: isothermal exp disk, 1: exp, 2: constant density, >2: exp(-R/h-h/R)
        read(50,*) nd(i),flaggd(i)
        read(50,*) mdisk(i)
        read(50,*) hd(i),zd(i),Rflzd(i)
        read(50,*) flagflzd(i)
        read(50,*) frdlim(i),fzdlim(i)
        rdlim(i) = frdlim(i)*hd(i)
        zdlim(i) = fzdlim(i)*zd(i)
        if(flaggd(i).eq.0) then
c *** stellar disk
          read(50,*) frq(i)
        else if(flaggd(i).gt.0) then
c gas disk
          read(50,*) fpjm
          read(50,*) dfeh(i),feh0(i),dalfe(i),alfe0(i)
          read(50,*) nhhalo(i),fehhalo(i)
        else
          read(50,*) tdisk
        endif
      enddo
      close(50)

      fb=omgb/omg0

c *** since ver. 6 ***
      fs=fb

c /*** cosmology ***/
      lam0 = 1.0d0-omg0
      iai = 1.0d0+zt
      xi = omg0*(iai**3)-(omg0+lam0-1.0d0)*(iai**2)+lam0
      omz = omg0*(iai**3)/xi

c /*** for the definition of NFW profile at z=0 ***/
      rvir = (1.63e-2*((mvir*hub)**(1.0d0/3.0d0))/hub)
     &  *((omg0/omz)**(-1.0d0/3.0d0))/(1.0d0+zt)     
      rt=xrt*rvir
c /*** Hernquist Profile with Boylan-Kolchin et al. 2005 parameters ***/
      re = refac*4.16d0*((mbul/1.0e11)**0.56d0)
      ah = re/1.8153d0

c *** set disksgen values
      eps=1.0e-3
      dsig=3.0
c necessary for thick disk model
      nr=1024
      nth=64
! nz should be high to sample steep decrease of density
      nz=512

c      nr=1024
c      nth=64
c      nz=64
c for test
      nr=256
      nth=32
      nz=32

c set total disk mass: mtdisk, mdisk is mass within rdlim
      flag=0
      do i=1,ndisk
        if(flaggd(i).le.1) then
c exp profile
c note radius unit is kpc, but normalised by hd
          mtdisk(i)=mdisk(i)*mrexpdf(rvir,hd(i),1.0d0)
     &      /mrexpdf(rdlim(i),hd(i),1.0d0)
        else if(flaggd(i).eq.2) then
c constant dnesity
          mtdisk(i)=mdisk(i)
        else
c check 
          if(flag.ne.0) then
            write(6,*) ' Error: flaggd>2 disk should be only one'
            call MPI_FINALIZE()
            stop
          endif
          flag=flag+1
c exp2d profile
c set total mass profile for exp2d 
          call setexp2dmr(mtdisk(i),mdisk(i),hd(i),rdlim(i))

        endif
      enddo

c output the input information
      if(myrank.eq.0) then
        open(60,file='inipara.dat',status='unknown')
        write(60,*) ' ver. 37 15/12/2015'
        write(6,*) ' ver. 37 15/12/2015'
        write(60,*) '###   NFW halo parameters  ###' 
        write(6,*) '###   NFW halo parameters  ###' 
        if(flagmexth.eq.0) then
          write(60,*) ' number of halo particles = ',nh
          write(6,*) ' number of halo particles = ',nh
          write(60,*) ' halo profile truncated with exp(r/rt)^2 '
          write(6,*) ' halo profile truncated with exp(r/rt)^2 '
        else if(flagmexth.eq.1) then
          write(60,*) ' halo is NFW halo'
          write(6,*) ' halo is NFW halo'
        else if(flagmexth.eq.-1) then
          write(60,*) ' + Bulge: Hernquist profile parameters'
          write(6,*) ' + Bulge: Hernquist profile parameters'
          write(60,*) ' Mbul (Msun)   =',mbul
          write(60,*) ' Re   (kpc)    =',re
          write(60,*) ' a    (kpc)    =',ah
          write(6,*) ' Mbul,Re,ah    =',mbul,re,ah
          write(6,*) ' nb set to be 0'
          nb=0
        else
          write(60,*) ' halo is constant density halo'
          write(60,*) ' radius -',cnfw
          write(6,*) ' halo is constant density halo'
          write(6,*) ' radius -',cnfw
        endif
        write(60,*) ' Mvir (Msun), cnfw = ',mvir,cnfw
        write(60,*) ' omg0,h,z = ',omg0,hub,zt,' lam0=1-omg0'
        write(60,*) ' omgb,fb=',omgb,fb
        write(60,*) ' maximum radius = rt x xcut, xcut,xrt='
     &   ,xcut,xrt
        write(60,*) ' truncation radius rt (kpc)=',rt
        write(60,*) ' rcut (kpc)=',rt*xcut
        write(60,*) ' r200 (kpc)  =',rvir
        write(60,*) ' r200 (kpc/h)=',rvir*hub
        write(60,*) '  c          =',cnfw
c
        write(6,*) ' Mvir (Msun), cnfw = ',mvir,cnfw
        write(6,*) ' omg0,h,z = ',omg0,hub,zt,' lam0=1-omg0'
        write(6,*) ' omgb,fb=',omgb,fb
        write(6,*) ' maximum radius = rt x xcut, xcut,xrt='
     &   ,xcut,xrt
        write(6,*) ' truncation radius rt (kpc)=',rt
        write(6,*) ' rcut (kpc)=',rt*xcut
        write(6,*) ' r200 (kpc)  =',rvir
        write(6,*) ' r200 (kpc/h)=',rvir*hub
        write(6,*) '  c          =',cnfw
        if(nb.gt.0) then
          write(60,*) '###   Bulge: Hernquist profile parameters  ###'
          write(60,*) ' Np bulge      =',nb
          write(60,*) ' Mbul (Msun)   =',mbul
          write(60,*) ' Re   (kpc)    =',re
          write(60,*) ' a    (kpc)    =',ah
c
          write(6,*) '###   Bulge: Hernquist profile parameters  ###'
          write(6,*) ' Np bulge      =',nb
          write(6,*) ' Mbul (Msun)   =',mbul
          write(6,*) ' Re   (kpc)    =',re
          write(6,*) ' a    (kpc)    =',ah
        endif
        if(nhg.gt.0) then
          write(60,*) '###   Halo Hot Gas  ###'
          write(60,*) ' Np hotgas   =',nhg
          write(60,*) ' Mhg  (Msun) < rvir =',mtothg
          write(60,*) ' Lam spin    =',lamjhg
          write(60,*) ' hot gas generated up to R =',rvir*xrhg
c       
          write(6,*) '###   Halo Hot Gas  ###'
          write(6,*) ' Np hotgas   =',nhg
          write(6,*) ' Mhg  (Msun) =',mtothg
          write(6,*) ' Lam spin    =',lamjhg
c convert to GCD+ unit
          mtothg=mtothg/MUSM    
        endif
        do i=1,ndisk
          if(nd(i).gt.0) then
            write(60,*) '###   Disk:',i,' Exp profile parameters  ###'
            write(60,*) ' Mdisk total (Msun)    =',mtdisk(i)
            write(60,*) ' Mdisk (<rdisk,Msun)    =',mdisk(i)
            write(60,*) ' h  (kpc)        =',hd(i)
            write(60,*) ' z0 (kpc)        =',zd(i)
            if(flagflzd(i).ne.0) then
              write(60,*) ' zd = z0 exp(R/',Rflzd(i),')'
            endif
            write(60,*) ' rdisk_max (kpc) =',frdlim(i),'x h =',rdlim(i)
            write(60,*) ' zdisk_max (kpc) =',fzdlim(i),'x z0=',zdlim(i)
c
            write(6,*) '###   Disk:',i,' Exp profile parameters  ###'
            write(6,*) ' Mdisk total (Msun)    =',mtdisk(i)
            write(6,*) ' Mdisk (<rdisk,Msun)    =',mdisk(i)
            write(6,*) ' h  (kpc)        =',hd(i)
            write(6,*) ' z0 (kpc)        =',zd(i)
            if(flagflzd(i).ne.0) then
              write(6,*) ' zd = z0 exp(R/',Rflzd(i),')'
            endif
            write(6,*) ' rdisk_max (kpc) =',frdlim(i),'x h =',rdlim(i)
            write(6,*) ' zdisk_max (kpc) =',fzdlim(i),'x z0=',zdlim(i)
            if(flagflzd(i).ne.0
     &         .and.zdlim(i).lt.zd(i)*exp(rdlim(i)/Rflzd(i))) then
              write(6,*) ' Error zd lim is too small.'
              write(6,*) ' zdlim should be bigger than zd at rdlim='
     &          ,zd(i)*exp(rdlim(i)/Rflzd(i))
              stop
            endif
            if(flaggd(i).eq.0) then
              write(60,*) ' stellar disk'
              write(60,*) ' fr=',frq(i)
c
              write(6,*) ' stellar disk'
              write(6,*) ' fr=',frq(i)
            else
              write(60,*) ' gas disk flaggd=',flaggd(i)
              write(6,*) ' gas disk flaggd=',flaggd(i)
              if(flaggd(i).gt.0) then
                write(60,*) ' factor for jeans mass limit =',fpjm
                write(60,*) ' ETAH =',ETAH
                write(6,*) ' factor for jeans mass limit =',fpjm
                write(6,*) ' ETAH =',ETAH
                if(flaggd(i).eq.2) then
                  write(60,*) ' constant density'
                  write(60,*) ' r and z initial limit=',hd(i),zd(i)
                  write(60,*) ' rhoc,nh='
     &           ,mdisk(i)/(M_PI*(hd(i)**2)*(2.0d0*zd(i)))
     &            *((LUKPC**3)/MUSM)
     &           ,0.76*mdisk(i)*(DU/MP)/(M_PI*(hd(i)**2)*(2.0d0*zd(i)))
     &            *((LUKPC**3)/MUSM)
                  write(6,*) ' constant density'
                  write(6,*) ' r and z initial limit=',hd(i),zd(i)
                  write(6,*) ' rhoc,nh='
     &           ,mdisk(i)/(M_PI*(hd(i)**2)*(2.0d0*zd(i)))
     &            *((LUKPC**3)/MUSM)
     &           ,0.76*mdisk(i)*(DU/MP)/(M_PI*(hd(i)**2)*(2.0d0*zd(i)))
     &            *((LUKPC**3)/MUSM)
                else
                  write(60,*) ' Exp(-R/h-h/R) profile'
                  write(6,*) ' Exp(-R/h-h/R) profile'
                endif
              else
                write(60,*) ' isothermal T=',tdisk
                write(6,*) ' isothermal T=',tdisk
                dfeh(i)=0.0d0
                feh0(i)=0.0d0
                dalfe(i)=0.0d0
                alfe0(i)=0.0d0         
                nhhalo(i)=-100.0
                fehhalo(i)=0.0d0
              endif
              write(60,*) ' [Fe/H](R)=Rx',dfeh(i),'+',feh0(i)
              write(60,*) ' [al/Fe]=[Fe/H]x',dalfe(i),'+',alfe0(i)
              write(60,*) ' nH density limit for disk nHlim=',nhhalo(i)
              write(60,*) ' if nH<nHlim, [Fe/H]=',fehhalo(i)
c
              write(6,*) ' [Fe/H](R)=Rx',dfeh(i),'+',feh0(i)
              write(6,*) ' [al/Fe]=[Fe/H]x',dalfe(i),'+',alfe0(i)
              write(6,*) ' nH density limit for disk nHlim=',nhhalo(i)
              write(6,*) ' if nH<nHlim, [Fe/H]=',fehhalo(i)
            endif
          endif
        enddo
        if(ndisk.gt.0) then
          write(60,*) ' in diskspgen'
          write(60,*) ' nr,nth,nz=',nr,nth,nz
          write(60,*) ' eps,dsig=',eps,dsig
          write(60,*) ' RHO_LIM=',RHO_LIM
c
          write(6,*) ' in diskspgen'
          write(6,*) ' nr,nth,nz=',nr,nth,nz
          write(6,*) ' eps,dsig=',eps,dsig
          write(6,*) ' RHO_LIM=',RHO_LIM
        endif
        close(60)
      endif
c /*** convert to my units ****/
      mvir = mvir/MUSM
      rvir = rvir/LUKPC
      rt=rt/LUKPC
      rcutdmh=rt*xcut
      rs = rvir/cnfw
      vvir=dsqrt(G*mvir/rvir)  

c *** for bulge ***
      if(nb.gt.0.or.flagmexth.eq.-1) then
        ah = ah/LUKPC   
        mbul = mbul/MUSM
      else
        mbul=0.0d0
      endif
c *** for disk ***
      do id=1,ndisk
        if(nd(id).gt.0) then
          hd(id)=hd(id)/LUKPC   
          zd(id)=zd(id)/LUKPC   
          Rflzd(id)=Rflzd(id)/LUKPC
          mdisk(id)=mdisk(id)/MUSM
          mtdisk(id)=mtdisk(id)/MUSM
          rdlim(id)=rdlim(id)/LUKPC
          zdlim(id)=zdlim(id)/LUKPC
        else
          hd(id)=0.0d0
          zd(id)=0.0d0
          Rflzd(id)=INF
          mdisk(id)=0.0d0
          mtdisk(id)=0.0d0
        endif
      enddo

c *** DM halo profile ***
      tmh=0.0d0
      if(abs(flagmexth).le.1) then
        call dmhaloprof(omg0,omz,hub,cnfw,rvir,zt,xcut,rt,fb
     &    ,flagmexth)
      else
        rvir=cnfw
        call dmhaloprof2(mvir,rvir,fb)
      endif
      tmh=mrdmh(ndr) 
      mdmhp=tmh/dble(nh)
      if(myrank.eq.0) then
        write(6,*) ' total mass of halo mass < rcut= ',tmh
        write(6,*) ' mass of particle, mdmhp= ',mdmhp
      endif
      if(nb.gt.0.or.flagmexth.eq.-1) then
c *** add hernquist bulge to an external potential ***
        call hpbulext(mbul,ah,xcut,rt,flagmexth)
        tmh=mrdmh(ndr) 
        if(myrank.eq.0) then
          write(6,*) ' after adding Hernquist bulge to extPot'
          write(6,*) ' total mass of halo+bulge mass < rcut= ',tmh
        endif
      endif
      if(nhg.gt.0) then
c *** add hot gas mass to an external potential ***
        call hgext(mtothg,cnfw,rvir,xcut,rt)
      endif
      if(flagmexth.ne.-1.and.myrank.eq.0) then
c use DM profile only for mext.dat
        call system('mv mext-dmh.dat mext.dat')
      endif

      tmtot=mvir*(1.0d0-fb)
      tmtot=tmtot+mbul
      bmtot=mbul
      do id=1,ndisk
        tmtot=tmtot+mdisk(id)
        bmtot=bmtot+mdisk(id)
      enddo
c      tmtot=tmtot+tmhg
c      bmtot=bmtot+tmhg
      if(myrank.eq.0) then
        write(6,*) ' initial and final baryon fraction ',fb,
     &    bmtot/tmtot
      endif
c *** generate sample disk particles ***
      idgdisk=-1
      np=0

      do id=1,ndisk
        idsds(id)=np
! not log grid for disks particles
        flaglog=0
        if(flaggd(id).ne.0.or.flagflzd(id).ne.0) then
! log grid for disks particles
          flaglog=1
        endif
        call diskspgen(np,rdlim(id),zdlim(id),flaglog,id
c        call diskspgen(np,rdlim(id),zdlim(id),0,id
     &   ,nr,nth,nz)
        ideds(id)=np-1
        call disksdenm(idsds(id),ideds(id),mtdisk(id)
     &   ,hd(id),zd(id),id,flaggd(id),flagflzd(id),Rflzd(id))
        if(flaggd(id).ne.0) then
c *** store gas disk id 
          idgdisk=id
        endif
      enddo
      if(myrank.eq.0) then
        write(6,*) myrank,'generated sample disk particle np=',np
      endif
      call setkernel()

c      do id=1,ndisk
c        write(filen,'(a5,i2.2,a1,i3.3,a4)') 'disks',id,'-',myrank,'.dat'
c        open(60,file=filen,status='unknown') 
c        do i=idsds(id),ideds(id)
c         write(60,'(7(1pE13.5))') x_p(i),y_p(i),z_p(i),m_p(i)
c     &    ,h_p(i),dvr_p(i),dsqrt(x_p(i)**2+y_p(i)**2)
c        enddo
c        close(60)
c      enddo

      if(myrank.eq.0) then
        if(idgdisk.ge.1) then
          write(6,*) ' gas disk id,h,z=',idgdisk,hd(idgdisk),zd(idgdisk)
          write(6,*) ' density limit=',RHO_LIM
        endif
      endif

      if(idgdisk.ge.0) then
        id=idgdisk
c *** calculate ptential structures ***
c *** generate reference grid for gas disk ***
        ngr=200
        ngz=200
        xcutgd=1.0d0
        frgd=frini0d
        ri = frgd*rdlim(id)
        ro = rdlim(id)*xcutgd
        zi = frini0dz*zdlim(id)
        zo = zdlim(id)*xcutgd
        call setgrid(npg,ngr,ngz,ri,ro,zi,zo)
        if(myrank.eq.0) then
           write(6,*) ' Ngrid for gas disk =',npg
        endif
c *** set initial density profile ***
        call diskdeng(npg,mtdisk(id),hd(id),zd(id),flaggd(id)
     &   ,flagflzd(id),Rflzd(id))
c
c        open(60,file='diskden.dat',status='unknown')
c        do i=0,npg-1
c          write(60,'(3(1pE13.5))') x_rzg(i),z_rzg(i),rho_rzg(i)
c        enddo
c
c *** set metal profile ***
        call diskmetg(ngr,ngz,dfeh(id),feh0(id)
     &   ,dalfe(id),alfe0(id),nhhalo(id),fehhalo(id))
        call seteos(zred)
        call setcool()
c *** set Jeans mass limit factor
        mdp=mdisk(id)/dble(nd(id))
        if(flaggd(id).le.0) then
c *** isothermal case if fpjm<0, then -fpjm=temperature of gas disk
          fpjm=-tdisk
        endif

c *** find gas density profile iteratively ***
        nit=10
        flagexit=0
        do it=1,nit
c *** get dphi/dz
          if(np.gt.0) then
            call treebuild(np)
            call treeforce(npg)
          endif
c disksf
c       write(filen,'(a4,i2.2,a1,i3.3,a4)') 'grzg',id,'-',myrank,'.dat'
c       open(60,file=filen,status='unknown') 
c       do i=0,npg-1
c         write(60,'(6(1pE13.5))') x_rzg(i),z_rzg(i),rho_rzg(i)
c     &    ,fx_rzg(i),fy_rzg(i),fz_rzg(i)
c       enddo
c       close(60)
c       call MPI_FINALIZE()

          call extf(npg,x_rzg,y_rzg,z_rzg,fxe,fye,fze)
          do i=0,npg-1
            fx_rzg(i)=fx_rzg(i)+fxe(i)        
c ignore position fxyz
            if(fx_rzg(i).gt.0d0) then
              fx_rzg(i)=0.0d0
            endif
            fy_rzg(i)=fy_rzg(i)+fye(i)        
            if(fy_rzg(i).gt.0d0) then
              fy_rzg(i)=0.0d0
            endif
            fz_rzg(i)=fz_rzg(i)+fze(i)        
            if(fz_rzg(i).gt.0d0) then
              fz_rzg(i)=0.0d0
            endif
          enddo
c          if(it.gt.1.and.errm.lt.EPSJ) then
          if(it.gt.1.and.errmax.lt.EPSJ) then
c            goto 72
            flagexit=1
          endif
c *** integrate hydrostatic eq.
          call hydstat(ngr,ngz,errmax,errm,it,fpjm,mdp,zred)
          if(myrank.eq.0) then
            write(6,*) ' errmax,mean=',errmax,errm,it
          endif
c *** adjust mass for gas disk
          call ndenggrid(ngr,ngz,mtdisk(id),hd(id)
     &      ,idsds(id),ideds(id),fpjm,mdp,zred,dsig
     &      ,flagexit,flaggd(id))
          if(flagexit.ne.0) then
             goto 72
          endif
c *** re-set metal profile ***
          call diskmetg(ngr,ngz,dfeh(id),feh0(id)
     &     ,dalfe(id),alfe0(id),nhhalo(id),fehhalo(id))
        enddo
        if(myrank.eq.0) then
          write(6,*) ' no convergent for gas density errmax=',errmax
          open(60,file='nit-gas.dat') 
          write(60,*) 'nit=',nit
          write(60,*) 'errmax,errm=',errmax,errm
          close(60)
        endif
 72     if(myrank.eq.0) then
          write(6,*) ' finish generating gas density map nit=',it  
        endif
c for test
        if(myrank.eq.0) then
    
          open(60,file='gasgrid.dat',status='unknown')
c          open(60,file=filen,status='unknown')
          do i=0,ngr
            do j=0,ngz
            if(i.eq.0) then
              pnpi=id_rzg(i,j)
              pnpo=id_rzg(i+1,j)
            else if(i.ne.ngr) then
              pnpi=id_rzg(i-1,j)
              pnpo=id_rzg(i+1,j)
            else 
              pnpi=id_rzg(i-1,j)
              pnpo=id_rzg(i,j)
            endif
            dr=x_rzg(pnpo)-x_rzg(pnpi)
            dpdr=(p_rzg(pnpo)-p_rzg(pnpi))/dr
            pn=id_rzg(i,j)
            write(60,'(11(1pE13.5))') x_rzg(pn),z_rzg(pn),rho_rzg(pn)
     &       ,dlog10(tm_rzg(pn)),lmet_rzg(pn),p_rzg(pn)
     &       ,dlog10(rho_rzg(pn)*fh_rzg(pn)*(DU/MP))
     &       ,p_rzg(pn)/((GAM-1.0d0)*rho_rzg(pn))
     &       ,myu_rzg(pn),dpdr,mz_rzg(pn)
            enddo
          enddo
          close(60)
        endif
c *** kinematics 
c *** for grid
        call vgasgrid(ngr,ngz)
        call gasdgen(nd(id),mtdisk(id),mdisk(id),hd(id)
     &    ,rdlim(id),zdlim(id),ngr,ngz,id,flagb
     &    ,dfeh(id),feh0(id),dalfe(id),alfe0(id),flaggd(id))
      endif

c *** calculate velocity structures of stellar disks ***
      flagvcf=0
      do id=1,ndisk
        if(flaggd(id).eq.0) then
          ngr=200
          ngz=200
c too high ngr increase noise. 200 is good enough after convergence test
c 03/07/2011
          xcutgd=1.0d0
          frgd=frini0d
          ri = frgd*rdlim(id)
          ro = rdlim(id)*xcutgd
          zi = frini0dz*zdlim(id)
          zo = zdlim(id)*xcutgd
c          ri = frini0*rt
c          ro = rvir*xcutgd
c          zi = frini0*rvir
c          zo = rvir*xcutgd

          call setgrid(npg,ngr,ngz,ri,ro,zi,zo)
c *** set initial density profile ***
          call diskdeng(npg,mtdisk(id),hd(id),zd(id),flaggd(id)
     &   ,flagflzd(id),Rflzd(id))
c *** get dphi/dz from sample disk particles (whose mass is adjusted to follow the calculated mass distribution of stars and the gas ***
          if(np.gt.0) then
            call treebuild(np)
            call treeforce(npg)
          endif

          call extf(npg,x_rzg,y_rzg,z_rzg,fxe,fye,fze)

! output the contribution from disk and DM and bulge
          if(myrank.eq.0.and.flagvcf.eq.0) then
            flagvcf=1
            open(60,file="vccont.dat",status='unknown') 
            j=0
            do i=0,ngr
              pn=id_rzg(i,j)
              write(60,'(4(1pE13.5))') x_rzg(pn)*LUKPC
     &         ,VUKMS*dsqrt(-x_rzg(pn)*(fx_rzg(pn)+fxe(pn)))
     &         ,VUKMS*dsqrt(-x_rzg(pn)*fx_rzg(pn))
     &         ,VUKMS*dsqrt(-x_rzg(pn)*fxe(pn))
            enddo
          endif

          do i=0,npg-1
            fx_rzg(i)=fx_rzg(i)+fxe(i)        
c ignore position fxyz
            if(fx_rzg(i).gt.0d0) then
              fx_rzg(i)=0.0d0
            endif
            fy_rzg(i)=fy_rzg(i)+fye(i)        
            if(fy_rzg(i).gt.0d0) then
              fy_rzg(i)=0.0d0
            endif
            fz_rzg(i)=fz_rzg(i)+fze(i)        
            if(fz_rzg(i).gt.0d0) then
              fz_rzg(i)=0.0d0
            endif
          enddo
c
          if(myrank.eq.0) then
          write(filen,'(a6,i3.3,a4)') 'sfgrid',id,'.dat' 
          open(60,file=filen,status='unknown')
          do i=0,npg-1
            write(60,'(4(1pE13.5))') x_rzg(i),z_rzg(i)
     &       ,-fx_rzg(i),-fz_rzg(i)
          enddo
          endif

          if(myrank.eq.0) then
c *** set velocity structure ***
          call vgrid(ngr,ngz,frq(id),mtdisk(id),hd(id),id)
c *** generate particles ***

          call expdgen(id,nd(id),mtdisk(id),mdisk(id),hd(id),zd(id)
     &      ,rdlim(id),zdlim(id),ngr,ngz,flagb
     &      ,omg0,omz,hub,cnfw,rvir,zt,fb,flagmexth
     &      ,flagflzd(id),Rflzd(id))
          write(6,*) ' generate disk',id,' nd=',nd(id)
          endif
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      enddo

c *** generating live bulge particles.
      if(nb.gt.0.and.flagmexth.ne.-1) then
c bulge particle mass
        mbp=mbul/dble(nb)
        if(myrank.eq.0) then
          write(6,*) ' generating bulge particles: mbp=',mbp
        endif
        ngr=600
        ngz=600
        xcutgd=1.0d0
        frgd=frini0d
        ri = frgd*rcutdmh
        ro = xcutgd*rcutdmh
        zi = frini0dz*rcutdmh
        zo = xcutgd*rcutdmh

        call setgrid(npg,ngr,ngz,ri,ro,zi,zo)
c *** set initial density profile ***
        call buldengen(npg,mbul,ah)
c *** get dphi/dz
        if(np.gt.0) then
          call treebuild(np)
          call treeforce(npg)
        endif
c        
        call extf(npg,x_rzg,y_rzg,z_rzg,fxe,fye,fze)
        do i=0,npg-1
          fx_rzg(i)=fx_rzg(i)+fxe(i)        
c ignore position fxyz
          if(fx_rzg(i).gt.0d0) then
            fx_rzg(i)=0.0d0
          endif
          fy_rzg(i)=fy_rzg(i)+fye(i)        
          if(fy_rzg(i).gt.0d0) then
            fy_rzg(i)=0.0d0
          endif
          fz_rzg(i)=fz_rzg(i)+fze(i)        
          if(fz_rzg(i).gt.0d0) then
            fz_rzg(i)=0.0d0
          endif
        enddo

c *** set velocity structure ***
        call bulvgrid(ngr,ngz)
c *** generate particles ***
        call bulpgen(nb,mbp,ngr,ngz,mbul,ah,rcutdmh,flagb)
      endif

c *** generating live halo hot gas particles.
      if(nhg.gt.0) then
        if(myrank.eq.0) then
          write(6,*) '### Generating Hot gas ###'
        endif
c        ngr=600
c        ngz=600
        ngr=200
        ngz=200
        xcutgd=1.0d0
        frgd=frini0d
        ri = frgd*rcutdmh
        ro = xcutgd*rcutdmh
        zi = frini0dz*rcutdmh
        zo = xcutgd*rcutdmh

        call setgrid(npg,ngr,ngz,ri,ro,zi,zo)
c *** set initial density profile ***

        call hgdengen(npg,mtothg,cnfw,rvir)

c *** get dphi/dz
        if(np.gt.0) then
          call treebuild(np)
          call treeforce(npg)
        endif
c        
        call extf(npg,x_rzg,y_rzg,z_rzg,fxe,fye,fze)
        do i=0,npg-1
          fx_rzg(i)=fx_rzg(i)+fxe(i)        
c ignore position fxyz
          if(fx_rzg(i).gt.0d0) then
            fx_rzg(i)=0.0d0
          endif
          fy_rzg(i)=fy_rzg(i)+fye(i)        
          if(fy_rzg(i).gt.0d0) then
            fy_rzg(i)=0.0d0
          endif
          fz_rzg(i)=fz_rzg(i)+fze(i)        
          if(fz_rzg(i).gt.0d0) then
            fz_rzg(i)=0.0d0
          endif
        enddo

c *** set velocity structure ***
c        call hgvgrid(ngr,ngz)
c set temperature profile, rmhg and etc. set in constant ln(r) bin. 
c and get mtothgrc: hot gas mass within rcuthg
        rcuthg=xrhg*rvir
        call hgtemprof(mtothg,cnfw,rvir,rcuthg,mtothgrc)
c *** get angular momentum ***
        call hgsetang(mtothg,cnfw,mvir,rvir,lamjhg,ajzhg)
c *** generate particles ***
c hot gas particle mass
        mhgp=mtothgrc/dble(nhg)
        if(myrank.eq.0) then
          write(6,*) ' generating hot gas particles: mhgp=',mhgp
        endif
        call hgpgen(nhg,mhgp,ngr,ngz,mtothg,cnfw,rvir,ajzhg 
     &    ,rcuthg,mtothgrc,flagb)
      endif

c *** generating DM halo particles.
      if(flagmexth.eq.0) then
        if(myrank.eq.0) then
          write(6,*) ' generating DM halo'          
        endif
        ngr=600
        ngz=600
        xcutgd=1.0d0
        frgd=frini0d
        ri = frgd*rcutdmh
        ro = xcutgd*rcutdmh
        zi = frini0dz*rcutdmh
        zo = xcutgd*rcutdmh

        call setgrid(npg,ngr,ngz,ri,ro,zi,zo)
c *** set initial density profile ***
        call dmhdengen(npg,omg0,omz,hub,cnfw,rvir,zt,xcut,rt,fb
     &    ,flagmexth)
c *** get dphi/dz
        if(np.gt.0) then
          call treebuild(np)
          call treeforce(npg)
        endif
c        
        call extf(npg,x_rzg,y_rzg,z_rzg,fxe,fye,fze)
        do i=0,npg-1
          fx_rzg(i)=fx_rzg(i)+fxe(i)        
c ignore position fxyz
          if(fx_rzg(i).gt.0d0) then
            fx_rzg(i)=0.0d0
          endif
          fy_rzg(i)=fy_rzg(i)+fye(i)        
          if(fy_rzg(i).gt.0d0) then
            fy_rzg(i)=0.0d0
          endif
          fz_rzg(i)=fz_rzg(i)+fze(i)        
          if(fz_rzg(i).gt.0d0) then
            fz_rzg(i)=0.0d0
          endif
        enddo
c
c        open(60,file='dmfgrid.dat',status='unknown')
c        do i=0,npg-1
c          write(60,'(4(1pE13.5))') x_rzg(i),z_rzg(i)
c     &     ,-fx_rzg(i),-fz_rzg(i)
c        enddo

c *** set velocity structure ***
        call dmhvgrid(ngr,ngz)
c reset mdmhr
        call dmhaloprof(omg0,omz,hub,cnfw,rvir,zt,xcut,rt,fb
     &    ,flagmexth)
c *** generate particles ***
        call dmhpgen(nh,mdmhp,ngr,ngz,omg0,omz,hub,cnfw
     &    ,rvir,zt,rcutdmh,rt,fb,flagb,flagmexth)

        if(myrank.eq.0) then
          write(6,*) ' finish generating the particles, nh=',nh
        endif
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE()
      stop
      end

