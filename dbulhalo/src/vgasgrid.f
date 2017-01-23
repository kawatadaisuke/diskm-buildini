c *****************************************************
c    vgrid.f
c 19 Jan., 2011   written by D. Kawata
c ***************************************************** 
c70*******************************************************************

      subroutine vgasgrid(ngr,ngz)
      include './define.f'
      integer ngr,ngz,npg
      integer i,j,k,pn,pnpi,pnpo,ii,io
      double precision rp,ri,ro,dpdr,dr

      do j=0,ngz
        do i=0,ngr
          pn=id_rzg(i,j)
          rp=x_rzg(pn)
c *** set vsig=0
          vsigr_rzg(pn)=0.0d0
          vsigz_rzg(pn)=0.0d0
          vsigph_rzg(pn)=0.0d0
          if(i.eq.0) then
            vph_rzg(pn)=0.0d0
            vc2_rzg(pn)=0.0d0
          else
c *** get rotation velocity 
c ***dP/dR
c            if(i.ne.ngr) then
c              pnpi=id_rzg(i-1,j)
c              pnpo=id_rzg(i+1,j)
c            else 
              pnpi=id_rzg(i-1,j)
              pnpo=id_rzg(i,j)
c            endif
            dr=x_rzg(pnpo)-x_rzg(pnpi)
            dpdr=(p_rzg(pnpo)-p_rzg(pnpi))/dr
            if(fx_rzg(pn).lt.0.0d0) then
c since ver. 25
              if(dpdr.lt.0.0d0) then
                vc2_rzg(pn)=x_rzg(pn)*(-fx_rzg(pn)
     &            +(dpdr/rho_rzg(pn)))
              else
                vc2_rzg(pn)=x_rzg(pn)*(-fx_rzg(pn))
              endif
              if(vc2_rzg(pn).gt.0.0d0) then
                vph_rzg(pn)=dsqrt(vc2_rzg(pn))
              else
                vph_rzg(pn)=0.0d0
                vc2_rzg(pn)=0.0d0
              endif
            else
              vph_rzg(pn)=0.0d0
              vc2_rzg(pn)=0.0d0
            endif
          endif
        enddo
      enddo

c *** get sigma phi
      if(myrank.eq.0) then
      open(60,file='vgasgridz0.dat',status='unknown')
      open(61,file='vgasgrid.dat',status='unknown')
      npg=(ngr+1)*(ngz+1)
      do j=0,ngz
c       j=0
        do i=0,ngr
        pn=id_rzg(i,j)
        if(j.eq.0) then
          write(60,'(8(1pE13.5))') x_rzg(pn),y_rzg(pn),z_rzg(pn)
     &       ,dsqrt(vc2_rzg(pn))*VUKMS,vph_rzg(pn),-fx_rzg(pn)
     &       ,p_rzg(pn),rho_rzg(pn)
        endif
            write(61,'(7(1pE13.5))') x_rzg(pn),y_rzg(pn),z_rzg(pn)
     &       ,dsqrt(vc2_rzg(pn))*VUKMS,vph_rzg(pn),-fx_rzg(pn)
     &       ,p_rzg(pn)
        enddo
      enddo
      close(60)
      close(61)
      endif

      return
      end
