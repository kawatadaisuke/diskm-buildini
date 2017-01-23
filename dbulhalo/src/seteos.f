c *********************************************
c  seteos.f for plotting cooling rate
c  25 Jan. 2011   written by D. Kawata
c *********************************************
c70*******************************************************************

      subroutine  seteos()
      include 'define.f'
      integer i,j
      character rch*6
      double precision zred

      open(50,file='./ini/eos.dat',status='old')
      read(50,'(a3,1pE13.5)') rch,zred
      read(50,*)
      read(50,'(a6,i5)') rch,nmet_eos
      read(50,'(a5,i5)') rch,nnh_eos

      if(myrank.eq.0) then
        write(6,*) ' redshift from eos.dat=',zred
      endif

      do j=0,nmet_eos-1
        do i=0,nnh_eos-1
          read(50,150) lnh_eos(i),ltm_eos(i,j),myu_eos(i,j)
     &      ,lmet_eos(j)
 150      format(4(1pE13.5))
        enddo
      enddo
      close(50)

      return
      end

c *** calculate p and gam ***
      subroutine pgameos(rhop,lmetp,pp,gamp,temp,lnhp,myupi)
      include 'define.f'
      double precision rhop,lmetp,pp,gamp,temp,myupi
      integer i,j,inh,im,it,id
      double precision lnhp,ltmp,myup,hfrac
      double precision dlmet,dlnh
      double precision lnrhol,lnrhoh,lnpl,lnph,rhoi
      double precision wm1,wm2,wd1,wd2
c /*** Solar Abundance (meteorites from WW95) in terms of mass ***/
c /* 1  2   3   4   5   6    7    8    9    10   11  12   13   14   15
c /* 1H,4He,12C,14N,16O,20Ne,23Na,24Mg,27Al,28Si,32S,36Ar.40Ca,56Fe,58Ni */
      double precision xsol(15),XSZ
      data xsol/0.706d0,0.275d0,3.03d-3,1.11d-3,9.59d-3,1.62d-3,
     & 3.34d-5,5.15d-4,5.81d-5,6.53d-4,3.96d-4,7.74d-5,5.99d-5,
     & 1.17d-3,4.94d-5/
      parameter (XSZ=0.019d0)

c *** set delta
      dlnh=lnh_eos(2)-lnh_eos(1)
      dlmet=lmet_eos(2)-lmet_eos(1)

c *** N_H 
c      lnhp=dlog10(rhop*(0.76d0+(xsol(1)-0.76d0)*(10.0d0**lmetp))
c     &  *(DU/MP))
c      write(6,*) rhop,lnhp,lmetp
c *** find metallicity
      im=int((lmetp-lmet_eos(0))/dlmet)
      if(im.lt.0) then
c *** use lowest metallicity values ***
        im=0
        wm1=1.0d0
        wm2=0.0d0
      else if(im.ge.nmet_eos-1) then
c *** use the value at nment_eos-1
        im=nmet_eos-2
        wm1=0.0d0
        wm2=1.0d0
      else
        wm1=(lmet_eos(im+1)-lmetp)/dlmet
        wm2=(lmetp-lmet_eos(im))/dlmet
      endif
c *** find density ***
      id=int((lnhp-lnh_eos(0))/dlnh)
      if(id.lt.0) then
        id=0
        wd1=(lnh_eos(id+1)-lnhp)/dlnh
        wd2=(lnhp-lnh_eos(id))/dlnh
c *** no extrapolation
c        wd1=1.0d0
c        wd2=0.0d0
      else if(id.ge.nnh_eos-1) then
        id=nnh_eos-2
        wd1=0.0d0
        wd2=1.0d0
      else
        wd1=(lnh_eos(id+1)-lnhp)/dlnh
        wd2=(lnhp-lnh_eos(id))/dlnh
      endif
c temperature and myu
      ltmp=wd1*wm1*ltm_eos(id,im)
     &  +wd1*wm2*ltm_eos(id,im+1)
     &  +wd2*wm1*ltm_eos(id+1,im)
     &  +wd2*wm2*ltm_eos(id+1,im+1)

      myup=wd1*wm1*myu_eos(id,im)
     &  +wd1*wm2*myu_eos(id,im+1)
     &  +wd2*wm1*myu_eos(id+1,im)
     &  +wd2*wm2*myu_eos(id+1,im+1)
      temp=10.0d0**ltmp

c      if(lnhp.gt.-1.0d0.and.lnhp.lt.-0.5d0) then
c      if(lmetp.gt.-1.5d0.and.lmetp.lt.-1.0d0) then
c      open(90,file='testeos.dat',access='append')
c      write(90,'(4(1pE13.5))') lnhp,lmetp,ltmp,myup
c      close(90)
c      endif
c      endif

      myupi=myup
c *** pressure
      pp=(10.0d0**ltmp)*rhop*TPRHO*MYU/(myup*TUK)
c *** get gam 
c *** low density one
      rhoi=(10.0d0**lnh_eos(id))*(MP/DU)
     &  /(0.76d0+(xsol(1)-0.76d0)*(10.0d0**lmetp))
      lnrhol=dlog(rhoi)
      ltmp=wm1*ltm_eos(id,im)+wm2*ltm_eos(id,im+1)
      myup=wm1*myu_eos(id,im)+wm2*myu_eos(id,im+1)
      lnpl=dlog((10.0d0**ltmp)*rhoi*TPRHO*MYU/(myup*TUK))
c *** high density one
      rhoi=(10.0d0**lnh_eos(id+1))*(MP/DU)
     &  /(0.76d0+(xsol(1)-0.76d0)*(10.0d0**lmetp))
      lnrhoh=dlog(rhoi)
      ltmp=wm1*ltm_eos(id+1,im)+wm2*ltm_eos(id+1,im+1)
      myup=wm1*myu_eos(id+1,im)+wm2*myu_eos(id+1,im+1)
      lnph=dlog((10.0d0**ltmp)*rhoi*TPRHO*MYU/(myup*TUK))
      gamp=(lnph-lnpl)/(lnrhoh-lnrhol)

      return
      end

