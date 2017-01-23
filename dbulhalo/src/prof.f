c /***************************************************
c   prof.f  for bulhalo
c  22 Aug. 2006  written by D.Kawata
c ****************************************************/
c70*******************************************************************

c /*****   NFW density profile ****/
      function rhonfwsw(r,c,rvir,rt,flag)
      include 'define.f'
      integer flag
      double precision rhonfwsw
      double precision r,c,rvir,rcut
      double precision x,t,a,n,ct,rt
 
      x=r/rvir
      t=rt*c/rvir

      if(flag.eq.0) then
        rhonfwsw=dexp(-r**2/rt**2)/(c*x*((1.0d0+c*x)**2))
      else
        rhonfwsw=1.0d0/(c*x*((1.0d0+c*x)**2))
      endif

      return
      end

      function mrnfw(r,c,rvir)
      include 'define.f'
      double precision mrnfw
      double precision r,c,rvir
 
      mrnfw=(rvir**3)*(4.0d0*M_PI)*((c*r+(c*r+rvir)*dlog(rvir))
     &   /(-c*r-rvir)+dlog(c*r+rvir))/(c**3)

      return
      end

c /*****   Hernquist mass profile  *****/
c /*****   M(r)  *****/
      function mrhp(r,ms,ah)
      include 'define.f'
      double precision mrhp
      double precision r,ms,ah

      mrhp = ms*(r**2)/((r+ah)**2)

      return
      end

c /*****   rho(r)  *****/
      function rhohp(r,ms,ah)
      include 'define.f'
      double precision rhohp
      double precision r,ms,ah

      rhohp = (ms/(2.0d0*M_PI))*(ah/r)/((r+ah)**3)

      return
      end

c /*****   r mrhp(r)  *****/
      function rmrhp(f,ms,ah)
      include 'define.f'
      double precision rmrhp
      double precision f,ms,ah

      rmrhp = ah*dsqrt(f)/(1.0d0-dsqrt(f))

      return
      end

c modified NFW by Maller & Bullock (2004) for halo hot gas
      function rhombhg(x)
      double precision rhombhg,x

      rhombhg=1.0d0/((x+0.75d0)*((1.0d0+x)**2))

      return 
      end

      function mrmbhg(x)
      double precision mrmbhg,x

      mrmbhg=9.0d0*dlog(1.0d0+(4.0d0/3.0d0)*x)-8.0d0*dlog(1.0d0+x)
     &  -4.0d0*x/(1.0d0+x)

      return 
      end


c /*****   Exponential disk mass profile  *****/
c /*****   M(r)  from particle *****/
      function mr3dexpd(r,nd,rlistd,r3dd,massp)
      include 'define.f'
      double precision mr3dexpd
      integer i,ip,ipp
      integer nd,rlistd(MN)
      double precision r,r3dd(MN),mrp,massp(MN)

      mrp = 0.0d0
      ipp = rlistd(1)
      do i=1,nd
        ip=rlistd(i)
        if(r3dd(ip).lt.r) then
          mrp=mrp+massp(ip)
        else
          if(ipp.eq.rlistd(1)) then
            ip = rlistd(2)
            mrp=massp(ip)+massp(ip)*(r-r3dd(ipp))/(r3dd(ip)-r3dd(ipp))
          else
            mrp=mrp+massp(ip)*(r-r3dd(ipp))/(r3dd(ip)-r3dd(ipp))
          endif
          goto 70
         endif
         ipp = ip
      enddo
 70   mr3dexpd = mrp
      if(mr3dexpd.lt.0.0d0) then
        mr3dexpd = 0.0d0
      endif

      return
      end

c /*****   M(<r)  assuming thin disk *****/
      function mrexpdf(r,hd,mdisk)
      include 'define.f'
      double precision mrexpdf
      double precision r,hd,mdisk

c      if(flagdp.eq.0) then
        mrexpdf = mdisk*(1.0d0-(1.0d0+r/hd)*dexp(-r/hd))
c      else
c        mrexpdf = mdisk*(1.0d0-(1.0d0+(r/hd)+0.5d0*((r/hd)**2)
c     &    +(1.0d0/6.0d0)*((r/hd)**3))*dexp(-r/hd))
c      endif

      return
      end
      
