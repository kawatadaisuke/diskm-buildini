c *****************************************************
c    dmhdengen.f
c 16 June., 2011   written by D. Kawata
c ***************************************************** 
c70*******************************************************************

      subroutine dmhdengen(npg,omg0,omz,hub,cnfw,rvir,zt,xcut,rt,fb
     &  ,flag)
      include './define.f'

      integer npg
      integer i,flag
      double precision omg0,omz,hub,cnfw,rvir,zt,xcut,rt,fb
      double precision rhoc,delc,facnfw
      double precision rp,zp
      double precision rhonfwsw
      external rhonfwsw

      rhoc = (3.0d0*HUB0*HUB0*hub*hub/(8.0d0*M_PI*G))
     &  *(omg0/omz)*((1.0d0+zt)**3)
      delc = (200.0d0/3.0d0)*(cnfw**3/(dlog(1.0d0+cnfw)
     &  -cnfw/(1.0d0+cnfw)))
      facnfw = delc*rhoc*(1.0d0-fb)

      do i=0,npg-1
        if(x_rzg(i).gt.0.0d0.or.z_rzg(i).gt.0.0d0) then
          rp=dsqrt(x_rzg(i)**2+z_rzg(i)**2)
        else
          rp=10.0d0**lri_rzg
        endif
        if(z_rzg(i).gt.0.0d0) then
          zp=z_rzg(i)
        else
          zp=10.0d0**lzi_rzg
        endif
        rho_rzg(i)=rhonfwsw(rp,cnfw,rvir,rt,flag)*facnfw
      enddo

      return
      end
