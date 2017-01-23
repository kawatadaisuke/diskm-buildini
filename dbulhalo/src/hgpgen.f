c /***************************************************
c   hgpgen.f  for dbulhalo
c   20 Apr. 2015  written by D.Kawata
c ****************************************************/
c70*******************************************************************

c /*****   generating halo hot gas *****/
      subroutine hgpgen(nhg,mhgp,ngr,ngz,mtothg,cnfw,rvir
     &  ,ajzhg,rcuthg,mtothgrc,flagb)
      include 'define.f'

      integer nhg,ngr,ngz
      integer flagb
      double precision mhgp,mtothg,cnfw,rvir,ajzhg
      double precision rcuthg,mtothgrc
      character filen*60
      integer i,j
      integer nrp,ndrm
      double precision rcut
      double precision ri,ro,lri,lro,dlr,rp,rr,mrp,rp2d
      double precision pfr,ph,th
c for density profile
      double precision gcv,rho0hg,rs
c *** for particles ***
      double precision xbp,ybp,zbp,vxbp,vybp,vzbp,rhobp,ubp
c *** for velocity ***
      integer ir,iz
      integer pnrizi,pnrizo,pnrozi,pnrozo
      double precision wrizi,wrizo,wrozi,wrozo
      double precision dr,dz
      double precision vsigz,vsigr,vsigph,vphm,vsig2p
      double precision vrp,vphp

c to measure profiles 
      integer ndpr,ipr
      double precision lripr,lropr,mdrpr(0:MNR),dlrpr,dvdr

c *** external function 
      double precision rhombhg
      external rhombhg
      integer idum
      real ran1,gasdev
      external ran1,gasdev

      idum = -111

      rs=rvir/cnfw
      gcv=9.0d0*dlog(1.0d0+(4.0d0/3.0d0)*cnfw)
     &   -8.0d0*dlog(1.0d0+cnfw)-(4.0d0*cnfw/(1.0d0+cnfw))
      rho0hg=mtothg/(4.0d0*M_PI*(rs**3)*gcv)

c *** generating particles ***
      if(myrank.eq.0) then
      
        write(6,*) ' hgpgen nhg=',nhg

        write(filen,'(a15)') 'output/hgas.dat'
        if(flagb.eq.0) then
          open(60,file=filen,status='unknown',form='unformatted')
          write(60) nhg
        else
          open(60,file=filen,status='unknown')
          write(60,'(a3,I10)') '#n=',nhg
        endif
      endif
c for checking profile
      lripr=dlog10(0.01d0)
      lropr=dlog10(rcuthg)
      ndpr=100
      dlrpr=(lropr-lripr)/dble(ndpr)
      do i=0,ndpr
        mdrpr(i)=0.0d0
      enddo

      do i=0,nhg-1
 67     pfr = dble(ran1(idum))
c *** search radius ***
        do j=1,ndr
          if(mrhg(j)/mtothgrc.gt.pfr) then
            goto 72
          endif
        enddo
c        write(6,*) ' Warning: in finding radius for disk'
c        write(6,*) '  pfr =',pfr
        goto 67
 72     if(j.eq.1) then
c          write(6,*) ' Warning: in finding radius for disk: j,mr,prf='
c     &      ,j,mrexpd(j),pfr
          j = 2
        endif
        pfr=pfr*mtothgrc
c bin size is set in ln(r) bin in hgtemprof.f
        rp=dexp(dlog(rmhg(j-1))
     &    +(dlog(rmhg(j))-dlog(rmhg(j-1)))
     &    *(dlog(pfr)-dlog(mrhg(j-1)))
     &    /(dlog(mrhg(j))-dlog(mrhg(j-1))))
        if(rp.gt.rcuthg) then
          goto 67
        endif
c to measure profile
        ir=int((dlog10(rp)-lripr)/dlrpr)
        if(ir.ge.0.and.ir.le.ndpr) then
          mdrpr(ir)=mdrpr(ir)+mhgp
        endif
        zbp=-rp+2.0d0*rp*dble(ran1(idum))
        ph =  2.0d0*M_PI*dble(ran1(idum))
        rp2d=dsqrt(rp**2-zbp**2)
        xbp = rp2d*dcos(ph)
        ybp = rp2d*dsin(ph)
c *** density ***
        rhobp=rho0hg*rhombhg(rp/rs)
c velocity
c no velocity dispersion for gas
c *** vphm from angular momentum
        if(rp2d.lt.rvir) then
          vphm=rp2d*ajzhg
        else
! the maximum angular momentum is rvir x ajzhg
          vphm=rvir*ajzhg
        endif
        vphp=vphm
        vrp=0.0d0
        vzbp=0.0d0
c *** convert vr,vph -> vx,vy ***
        if(rp2d.gt.0.0d0) then
          vxbp=vrp*xbp/rp2d-vphp*ybp/rp2d
          vybp=vrp*ybp/rp2d+vphp*xbp/rp2d
        else
          vxbp=0.0d0
          vybp=0.0d0
        endif
c thermal energy
        ubp=dexp(dlog(uhg(j-1))
     &    +(dlog(uhg(j))-dlog(uhg(j-1)))
     &    *(dlog(rp)-dlog(rmhg(j-1)))
     &    /(dlog(rmhg(j))-dlog(rmhg(j-1))))
c *** output 
        if(myrank.eq.0) then
          if(flagb.eq.0) then
            write(60) xbp,ybp,zbp,vxbp,vybp,vzbp,mhgp,rhobp,ubp
          else
            write(60,160) xbp,ybp,zbp,vxbp,vybp,vzbp
     &       ,mhgp,rhobp,ubp,dsqrt(xbp**2+ybp**2+zbp**2)*LUKPC
     &       ,vsigph,vsigr,vsigz,rp2d,ubp*((GAM-1.0d0)/TPRHO)*TUK
 160        format(15(1pE13.5))
          endif
        endif
      enddo
      close(60)

      if(myrank.eq.0) then
        open(60,file='hgpprof.dat',status='unknown')
        do i=0,ndpr 
          rp=10.0d0**(lripr+dlrpr*(dble(i)+0.5d0))
          dvdr=(4.0d0*M_PI/3.0d0)
     &      *((10.0d0**(lripr+dlrpr*dble(i+1)))**3
     &       -(10.0d0**(lripr+dlrpr*dble(i)))**3)
          write(60,'(3(1pE13.5))') rp*LUKPC,mdrpr(i)/dvdr
     &     ,rho0hg*rhombhg(rp/rs)
        enddo
        close(60)
      endif

      return
      end
