c ******************************************
c    cool.F  cooling  copy for gcd+ pv.33.6
c  8 Feb. 2011    written by D. Kawata
c ******************************************
c ***********************************************************
c  set net cooling and myu
c ************************************************************

      SUBROUTINE cool(rhop,pp,lnhp,lmetp,current_z,myup)
      include 'define.f'
      integer np
      DOUBLE PRECISION current_z

      integer i
      double precision cradp,hradp,radp
      double precision lowmet
      double precision rhop,pp,lnhp,lmetp,myup
      double precision lognhp,logTp,metp,logmetp
      double precision dnh,dmet,dt
      double precision wt1,wt2,wd1,wd2,wz1,wz2,dz
      double precision wm1,wm2
      integer pn,it,id,iz,im

c lowest metallicity
      lowmet=10.0d0**met_crtb(0)
c set delta nH, dmet, dT
      dnh=nh_crtb(2)-nh_crtb(1)
      dmet=met_crtb(2)-met_crtb(1)
      dt=t_crtb(2)-t_crtb(1)
c

c find redshift 
      if(current_z.lt.z_crtb(nz_crtb-1)) then
        do i=1,nz_crtb-1
          if(z_crtb(i).gt.current_z) then
            iz = i
            goto 70
          endif
        enddo
        if(myrank.eq.0) then
          write(6,*) ' cool(): cannot find z =',current_z
          write(6,*) ' it must be smaller than ',z_crtb(nz_crtb-1)
        endif
        iz = nz_crtb-1
   70   dz = z_crtb(iz)-z_crtb(iz-1)
c /*** weight for redshift ***/
        wz1 = (z_crtb(iz)-current_z)/dz
        wz2 = (current_z-z_crtb(iz-1))/dz
      else
        iz = nz_crtb-1
        wz1 = 0.0d0
        wz2 = 1.0d0
c *** just for test ***
c        do i =0,np-1
c          pn = plist(i)
c          ram(pn)=0.0d0
c        enddo
c        return
      endif

c for check
c      if(np.gt.5000) then
c        open(60,file='cool.dat',status='unknown')
c      endif
      
!ocl novrec
!cdir nodep
c *** nH ***
      Lognhp=lnhp
c *** find metallicity ***
      if(lmetp.lt.met_crtb(0)) then
        logmetp=met_crtb(0)
      else
        logmetp=lmetp
      endif
      im = int((logmetp-met_crtb(0))/dmet)
      if(im.lt.0) then
c *** use the vale at met_crtb(0) ***
        im = 0
        wm1=1.0d0
        wm2=0.0d0
        else if(im.ge.nmet_crtb-1) then
c *** use the vale at met_crtb(nmet_crtb-1) ***
        im = nmet_crtb-2
        wm1=0.0d0
        wm2=1.0d0     
      else
        wm1 = (met_crtb(im+1)-logmetp)/dmet
        wm2 = (logmetp-met_crtb(im))/dmet
c          if(dabs(wm1).gt.1.05d0.or.dabs(wm2).gt.1.05d0) then
c            write(6,*) ' Warning: extrapolation in cool()'
c            write(6,*) ' Zmet id,wm1,wd2,lognh=',id,wm1,wm2,logmetp
c            write(6,*) ' logZ im+1,im=',met_crtb(im+1),met_crtb(im)
c          endif
      endif
      id=int((lognhp-nh_crtb(0))/dnh)
      if(id.lt.0) then
c *** use the vale at nh_crtb(0) ***
        id = 0
c          wd1=1.0d0
c          wd2=0.0d0
      else if(id.ge.nnh_crtb-1) then
c *** use the vale at met_crtb(nmet_crtb-1) ***
        id = nnh_crtb-2
c          wd1=0.0d0
c          wd2=1.0d0
      endif
      wd1 = (nh_crtb(id+1)-lognhp)/dnh
      wd2 = (lognhp-nh_crtb(id))/dnh

c          if(dabs(wd1).gt.1.05d0.or.dabs(wd2).gt.1.05d0) then
c            write(6,*) ' Warning: extrapolation in cool()'
c            write(6,*) ' nh id,wd1,wd2,lognh=',id,wd1,wd2,lognhp
c            write(6,*) ' lognh id+1,id=',nh_crtb(id+1),nh_crtb(id)
c          endif
c        endif
c /*** temperature ***/
      logTp = dlog10((pp*myup
     &    /(rhop*TPRHO*MYU))*TUK)
C /*** find temperature ***/
      it = int((logTp-t_crtb(0))/dt)
      if(it.lt.0) then
        it = 0
      else if(it.ge.nt_crtb-1) then
        it = nt_crtb-2
      endif
      wt1 = (t_crtb(it+1)-logTp)/dt
      wt2 = (logTp-t_crtb(it))/dt
c        if(dabs(wt1).gt.1.0d0.or.dabs(wt2).gt.1.0d0) then
c          write(6,*) ' Warning: extrapolation in cool()'
c          write(6,*) ' T it,w1,w2,logT=',it,wt1,wt2,logTp
c          write(6,*) ' logT it+1,it=',t_crtb(it+1),t_crtb(it)
c        endif

c *** get myu ***
c /*** if low and high temperature, no extrapolation ***/
      if(logTp.lt.t_crtb(0)) then
        it = 0
        wt1=1.0d0
        wt2=0.0d0
      else if(logTp.gt.t_crtb(nt_crtb-1)) then
        it = nt_crtb-2
        wt1=0.0d0
        wt2=1.0d0
      endif
c /*** if low and high density, no extrapolation ***/
      if(lognhp.lt.nh_crtb(0)) then
        id=0
        wd1=1.0d0
        wd2=0.0d0
      else if(lognhp.gt.nh_crtb(nnh_crtb-1)) then
        id = nnh_crtb-2
        wd1=0.0d0
        wd2=1.0d0
      endif

      myup
c iz-1
c   id
c     im
c       it
     &    =wt1*wm1*wd1*wz1*(myu_crtb(it,im,id,iz-1))
c       it+1
     &    +wt2*wm1*wd1*wz1*(myu_crtb(it+1,im,id,iz-1))
c     im+1
c       it
     &    +wt1*wm2*wd1*wz1*(myu_crtb(it,im+1,id,iz-1))
c       it+1
     &    +wt2*wm2*wd1*wz1*(myu_crtb(it+1,im+1,id,iz-1))
c   id+1
c     im
c       it
     &    +wt1*wm1*wd2*wz1*(myu_crtb(it,im,id+1,iz-1))
c       it+1
     &    +wt2*wm1*wd2*wz1*(myu_crtb(it+1,im,id+1,iz-1))
c     im+1
c       it
     &    +wt1*wm2*wd2*wz1*(myu_crtb(it,im+1,id+1,iz-1))
c       it+1
     &    +wt2*wm2*wd2*wz1*(myu_crtb(it+1,im+1,id+1,iz-1))
c iz
c   id
c     im
c       it
     &    +wt1*wm1*wd1*wz2*(myu_crtb(it,im,id,iz))
c       it+1
     &    +wt2*wm1*wd1*wz2*(myu_crtb(it+1,im,id,iz))
c     im+1
c       it
     &    +wt1*wm2*wd1*wz2*(myu_crtb(it,im+1,id,iz))
c       it+1
     &    +wt2*wm2*wd1*wz2*(myu_crtb(it+1,im+1,id,iz))
c   id+1
c     im
c       it
     &    +wt1*wm1*wd2*wz2*(myu_crtb(it,im,id+1,iz))
c       it+1
     &    +wt2*wm1*wd2*wz2*(myu_crtb(it+1,im,id+1,iz))
c     im+1
c       it
     &    +wt1*wm2*wd2*wz2*(myu_crtb(it,im+1,id+1,iz))
c       it+1
     &    +wt2*wm2*wd2*wz2*(myu_crtb(it+1,im+1,id+1,iz))

      return
      end

