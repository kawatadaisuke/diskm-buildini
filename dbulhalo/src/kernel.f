c *****************************************************
c    kernel.F for GCD+ pv.30.7
c 20 Aug., 2008   written by D. Kawata
c ***************************************************** 
c70*******************************************************************
      subroutine setkernel()

      include './define.f'
      integer i,j
      double precision f,cv,s

      ds_tb=1.0d0/dble(NKTAB+1)

c *** W: SPH kernel ***
c *** for constant ***
      cv=8.0d0/M_PI
      do i=0,NKTAB
        s_tb(i)=ds_tb*dble(i)
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=(1.0d0-6.0d0*(s**2)+6.0d0*(s**3))
        else if(s.le.1.0d0) then
          f=(2.0d0*((1.0d0-s)**3))
        else 
          f=0.0d0
        endif
        w_tb(i)=cv*f
      enddo
c *** dW(s)/ds/s: SPH kernel ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=-12.0d0+18.0d0*s
        else if(s.le.1.0d0) then
          f=-6.0d0*((1.0d0-s)**2)/s
        else 
          f=0.0d0
        endif
        dwds_s_tb(i)=cv*f
c *** core dw/ds ***
        if(s.le.THIRD) then
          f=-2.0d0
        else if(s.le.0.5d0) then
          f=-12.0d0*s+18.0d0*s*s
        else if(s.le.1.0d0) then
          f=-6.0d0*((1.0d0-s)**2)
        else 
          f=0.0d0
        endif
        dwdsc_tb(i)=cv*f
      enddo

c *** dphi/dr/r: phi kernel ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=(32.0d0/15.0d0)
     &     *(5.0d0-18.0d0*(s**2)+15.0d0*(s**3))
        else if(s.le.1.0d0) then
          f=(4.0d0/15.0d0)
     &     *(80.0d0*s-180.0d0*(s**2)+144.0d0*(s**3)
     &     -40.0d0*(s**4)-0.25d0/(s**2))/s
        endif
        dphidr_r_tb(i)=f
      enddo
c *** dphi/dh: phi kernel ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=-40.0d0*(s**2)+120.0d0*(s**4)-96.0d0*(s**5)+7.0d0
        else if(s.le.1.0d0) then
          f=-80.0d0*(s**2)+160.0d0*(s**3)-120.0d0*(s**4)
     &      +32.0d0*(s**5)+8.0d0
        endif
        dphidh_tb(i)=0.4d0*f
      enddo
c *** for multipole ***
c *** d2phi/dr2 ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=(32.0d0/15.0d0)
     &     *(5.0d0-54.0d0*(s**2)+60.0d0*(s**3))
        else if(s.le.1.0d0) then
          f=(4.0d0/15.0d0)
     &     *(80.0d0-360.0d0*s+432.0d0*(s**2)
     &     -160.0d0*(s**3)+2.0d0/(4.0d0*(s**3)))
        endif
        d2phidr2_tb(i)=f
      enddo
c *** d3phi/dr3 ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=(32.0d0/15.0d0)
     &     *(-108.0d0*s+180.0d0*(s**2))
        else if(s.le.1.0d0) then
          f=(4.0d0/15.0d0)
     &     *(-360.0d0+864.0d0*s
     &     -480.0d0*(s**2)-6.0d0/(4.0d0*(s**4)))
        endif
        d3phidr3_tb(i)=f
      enddo

      return
      end
