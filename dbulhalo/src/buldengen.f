c *****************************************************
c    buldengen.f
c 16 June., 2011   written by D. Kawata
c ***************************************************** 
c70*******************************************************************

      subroutine buldengen(npg,mbul,ah)
      include './define.f'

      integer npg
      integer i
      double precision mbul,ah
      double precision rp,zp
      double precision rhohp
      external rhohp

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
        rho_rzg(i)=rhohp(rp,mbul,ah)
      enddo

      return
      end
