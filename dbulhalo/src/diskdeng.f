c *****************************************************
c    diskden.f
c 15 Sep, 2012   written by D. Kawata
c ***************************************************** 
c70*******************************************************************

      subroutine diskdeng(npg,mdisk,hd,zd0,flaggd
     &  ,flagflzd,Rflzd)
      include './define.f'
      integer npg,flaggd,flagflzd
      double precision rp,zp
      double precision mdisk,hd,zd0,rhoc,Rflzd,zdr
      integer i

      if(flaggd.le.1) then
c exp density profile
        do i=0,npg-1
          if(x_rzg(i).gt.0.0d0) then
            rp=x_rzg(i)
          else
            rp=10.0d0**lri_rzg
          endif
          if(z_rzg(i).gt.0.0d0) then
            zp=z_rzg(i)
          else
            zp=10.0d0**lzi_rzg
          endif
          if(flagflzd.eq.0) then
            rho_rzg(i)=(mdisk/(4.0d0*M_PI*zd0*(hd**2)))
     &        *(1.0d0/(dcosh(z_rzg(i)/zd0)**2))
     &        *dexp(-x_rzg(i)/hd)
          else
c flaring
            zdr=zd0*dexp(x_rzg(i)/Rflzd)
            rho_rzg(i)=(mdisk/(4.0d0*M_PI*zdr*(hd**2)))
     &        *(1.0d0/(dcosh(z_rzg(i)/zdr)**2))
     &        *dexp(-x_rzg(i)/hd)
          endif
          if(rho_rzg(i).lt.0.0d0) then
            write(6,*) x_rzg(i),z_rzg(i),rho_rzg(i)
          endif
        enddo
      else if(flaggd.eq.2) then
c constant density
        rhoc=mdisk/(M_PI*(hd**2)*(2.0d0*zd0))
        do i=0,npg-1
          if(x_rzg(i).gt.0.0d0) then
            rp=x_rzg(i)
          else
            rp=10.0d0**lri_rzg
          endif
          if(z_rzg(i).gt.0.0d0) then
            zp=z_rzg(i)
          else
            zp=10.0d0**lzi_rzg
          endif
          if(rp.lt.hd.and.zp.lt.zd0) then
            rho_rzg(i)=rhoc
          else
            rho_rzg(i)=0.0d0
          endif
        enddo
      else 
c exp(-x/h-h/x) profile
        do i=0,npg-1
          if(x_rzg(i).gt.0.0d0) then
            rp=x_rzg(i)
          else
            rp=10.0d0**lri_rzg
          endif
          if(z_rzg(i).gt.0.0d0) then
            zp=z_rzg(i)
          else
            zp=10.0d0**lzi_rzg
          endif
          rho_rzg(i)=(mdisk/(3.18884*2.0d0*zd0*(hd**2)))
     &      *(1.0d0/(dcosh(z_rzg(i)/zd0)**2))
     &      *dexp(-x_rzg(i)/hd-hd/x_rzg(i))
          if(rho_rzg(i).lt.0.0d0) then
            write(6,*) x_rzg(i),z_rzg(i),rho_rzg(i)
          endif
        enddo
      endif

      return
      end
