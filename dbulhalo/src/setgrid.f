c ****************************************************
c  setgrid.f
c  25 jan. 2011    written by D.KAWATA
c ****************************************************

      subroutine setgrid(npg,ngr,ngz,ri,ro,zi,zo)
      include 'define.f'
      integer npg,ngr,ngz
      double precision ri,ro,zi,zo
      integer i,j
      double precision zz,rr

      npg=0
      lri_rzg = dlog10(ri)
      lro_rzg = dlog10(ro)
      dlr_rzg = (lro_rzg-lri_rzg)/dble(ngr-1)
      lzi_rzg = dlog10(zi)
      lzo_rzg = dlog10(zo)
      dlz_rzg = (lzo_rzg-lzi_rzg)/dble(ngz-1)
c *** generate particles ***
      zz=0.0d0
      npg=0
      do j=0,ngz
        if(j.ne.0) then
          zz=10.0d0**(lzi_rzg+dble(j-1)*dlz_rzg)
        else
          zz=0.0d0
        endif
        rr=0.0d0
        do i=0,ngr
          if(i.ne.0) then
            rr=10.0d0**(lri_rzg+dble(i-1)*dlr_rzg)
          else
            rr=0.0d0
          endif
          x_rzg(npg)=rr
          y_rzg(npg)=0.0d0
          z_rzg(npg)=zz
          id_rzg(i,j)=npg
          npg=npg+1
        enddo
      enddo

      return
      end
