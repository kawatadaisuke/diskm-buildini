c *****************************************************
c    hgdengen.f
c 22 Mar., 2015   written by D. Kawata
c ***************************************************** 
c70*******************************************************************

      subroutine hgdengen(npg,mtothg,cnfw,rvir)
      include './define.f'

      integer npg
      integer i
      double precision mtothg,cnfw,rvir
      double precision rs,gcv,rho0hg
      double precision rp,xp
      double precision rhombhg
      external rhombhg

      rs=rvir/cnfw
      gcv=9.0d0*dlog(1.0d0+(4.0d0/3.0d0)*cnfw)
     &   -8.0d0*dlog(1.0d0+cnfw)-(4.0d0*cnfw/(1.0d0+cnfw))
      rho0hg=mtothg/(4.0d0*M_PI*(rs**3)*gcv)

      do i=0,npg-1
        if(x_rzg(i).gt.0.0d0.or.z_rzg(i).gt.0.0d0) then
          rp=dsqrt(x_rzg(i)**2+z_rzg(i)**2)
        else
          rp=10.0d0**lri_rzg
        endif
        xp=rp/rs
        rho_rzg(i)=rho0hg*rhombhg(xp)
      enddo

      return
      end
