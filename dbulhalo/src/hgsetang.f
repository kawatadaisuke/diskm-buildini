c *****************************************************
c    hgsetang.f
c 23 Mar., 2015   written by D. Kawata
c ***************************************************** 
c70*******************************************************************

c calculating the constant, A, for specific angular momentum A R
c Use spin parameter in Bullock et al. (2001) 
c lambda = J/sqrt(2) Mvir Vvir Rvir
c for gas assume DM and gas has the same mass and j profile (which is not true)
c lambda = Jgas (M/Mgas) / sqrt(2) Mvir Vvir Rvir (eq. 1)
c Assume specific angular momentum of the shape of j = A R
c Jgas = Integrate[ Pi R Sig_gas(R) A R , {R, 0, rvir}]
c Then, A can be calculate from (eq. 1)

      subroutine hgsetang(mtothg,cnfw,mvir,rvir,lamjhg,ajzhg)
      include './define.f'
      integer NDRZ
      parameter (NDRZ=1000)
      double precision mtothg,cnfw,mvir,rvir,lamjhg,ajzhg 
      double precision vvir

      integer ngr,ngz,npg
      integer i,j,k,pn,pnpi,pnpo,ii,io
      double precision fin,dz,fp,fpp,dr,rp,r3dp,x3dp
      double precision gcv,rho0hg,rs
      double precision mgas2dr(0:NDRZ),mtot
      double precision ri,ro,lri,lro,dlr
      double precision zi,zo,lzi,lzo,dlz
      double precision rhorizi,rhorizo,rhorozi,rhorozo,rhop
      double precision jtotgas
      double precision drhovr
      double precision fr,qval,sdend
c external function
      double precision rhombhg
      external rhombhg

c constants for gas density profile
      rs=rvir/cnfw
      gcv=9.0d0*dlog(1.0d0+(4.0d0/3.0d0)*cnfw)
     &   -8.0d0*dlog(1.0d0+cnfw)-(4.0d0*cnfw/(1.0d0+cnfw))
      rho0hg=mtothg/(4.0d0*M_PI*(rs**3)*gcv)
c virial velocity 
      vvir=dsqrt(G*mvir/rvir)

      if(myrank.eq.0) then
        open(60,file='hgm2dprof.dat',status='unknown')
c        open(61,file='hgm3dprof.dat',status='unknown')
      endif

      ngz=NDRZ
      ngr=NDRZ
      ri=frini0*rvir
      ro=rvir
      lri=dlog10(ri)
      lro=dlog10(ro)
      dlr=(lro-lri)/dble(NDRZ)
c      write(6,*) ' ri,ro,rend=',ri,ro,10.0d0**(lri+dlr*dble(NDRZ))
c      write(6,*) ' lri,lro,dlr=',lri,lro,dlr
      ro=0.0d0
c angular momentum
      jtotgas=0.0d0
      mtot=0.0d0
      do i=0,ngr
        ri=ro
        ro=10.0d0**(lri+dlr*dble(i)) 
        if(i.gt.0) then
          rp=10.0d0**(lri+dlr*(dble(i)-0.5d0))
        else
          rp=0.50d0*(ri+ro)
        endif
        dr=ro-ri
c set z integration range
        zo=dsqrt(rvir**2-rp**2)
        zi=frini0*zo
        lzo=dlog10(zo)
        lzi=dlog10(zi)
        if(lzo.gt.lzi) then
          ngz=1+(lzo-lzi)/dlr
        else
          if(myrank.eq.0) then
            write(6,*) ' lzo>lzi at ir,rp,rvir=',i,rp,rvir
            write(6,*) ' zi,zo=',zi,zo
          endif
          stop
        endif
        dlz=(lzo-lzi)/dble(ngz)
        zo=0.0d0
        mgas2dr(i)=0.0d0
        do j=0,ngz
          zi=zo
          zo=10.0d0**(lzi+dlz*dble(j))
          dz=zo-zi
c get density
c rizi
          r3dp=dsqrt(ri**2+zi**2)
          x3dp=r3dp/rs
          rhorizi=rho0hg*rhombhg(x3dp)
c rizo
          r3dp=dsqrt(ri**2+zo**2)
          x3dp=r3dp/rs
          rhorizo=rho0hg*rhombhg(x3dp)
c rozi
          r3dp=dsqrt(ro**2+zi**2)
          x3dp=r3dp/rs
          rhorozi=rho0hg*rhombhg(x3dp)
c rozo
          r3dp=dsqrt(ro**2+zo**2)
          x3dp=r3dp/rs
          rhorozo=rho0hg*rhombhg(x3dp)
c mean rho
          rhop=0.25d0*(rhorizi+rhorizo+rhorozi+rhorozo)
c          if(myrank.eq.0) then
c            write(61,'(3(1pE13.5))') rp
c     &       ,10.0d0**(lzi+dlz*(dble(j)-0.5d0)),rhop
c          endif
c integration
          mgas2dr(i)=mgas2dr(i)+rhop*dz
        enddo
c surface mass density
        mgas2dr(i)=2.0d0*mgas2dr(i)
c summing up angular momentum A*rp 
        jtotgas=jtotgas+M_PI*rp*mgas2dr(i)*(ro**2-ri**2)
c total mass integration for check
        mtot=mtot+mgas2dr(i)*M_PI*(ro**2-ri**2)
c output
        if(myrank.eq.0) then
          write(60,'(5(1pE13.5),I8)') rp,mgas2dr(i),jtotgas,mtot
     &     ,10.0d0**lzo,ngz
        endif
      enddo
      ajzhg=lamjhg*dsqrt(2.0d0)*mtothg*rvir*vvir/jtotgas
      if(myrank.eq.0) then
        write(6,*) ' Total gas mass after Integration =',mtot
        write(6,*) ' Jtot_gas=',jtotgas
        write(6,*) ' Angular momentum profile is A R, A=',ajzhg
      endif
      close(60)

      return
      end
