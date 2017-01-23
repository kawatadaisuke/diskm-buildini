c *****************************************************
c    hgvgrid.f
c 22 Mar., 2015   written by D. Kawata
c ***************************************************** 
c70*******************************************************************

      subroutine hgvgrid(ngr,ngz)
      include './define.f'
      integer ngr,ngz,npg
      integer i,j,k,pn,pnpi,pnpo,ii,io
      double precision fin,dz,fp,fpp,dr,rp
      double precision vz2(0:MN2D-1,0:MN2D-1),vr2(0:MN2D-1,0:MN2D-1)
     &  ,vph2(0:MN2D-1,0:MN2D-1)
      double precision drhovr
      double precision fr,qval,sdend

      do j=0,ngz
        do i=0,ngr
          pn=id_rzg(i,j)
          rp=x_rzg(pn)
c *** set vsigz ***
c *** integration ***
          fin=0.0d0
          fp=-rho_rzg(pn)*fz_rzg(pn)
          do k=j+1,ngz
            fpp=fp
            pnpi=id_rzg(i,k-1)
            pnpo=id_rzg(i,k)
            dz=z_rzg(pnpo)-z_rzg(pnpi)
            fp=-rho_rzg(pnpo)*fz_rzg(pnpo)
            fin=fin+dz*(fp+fpp)*0.5d0
c            if(i.eq.1.and.j.eq.1) then
c              write(6,*), fin,dz,fp,fpp
c            endif
          enddo
          if(rho_rzg(pn).gt.0.0d0.and.fin.gt.0.0d0) then
            vz2(i,j)=fin/rho_rzg(pn)
          else
            vz2(i,j)=0.0d0
          endif
c *** set vsigz ***
          if(vz2(i,j).gt.0.0d0) then
            vsigz_rzg(pn)=dsqrt(vz2(i,j))
          else
            vsigz_rzg(pn)=0.0d0
          endif
c *** set vsigr2 ***
          vr2(i,j)=vz2(i,j)
          vsigr_rzg(pn)=dsqrt(vr2(i,j))

          if(i.gt.0) then
c *** circular velocity ***
            if(fx_rzg(pn).lt.0.0d0) then
              vc2_rzg(pn)=-x_rzg(pn)*fx_rzg(pn)
            else
              vc2_rzg(pn)=0.0d0
            endif
          else
             vc2_rzg(pn)=0.0d0
          endif
        enddo
      enddo

c *** get sigma phi
      open(60,file='vph-hg.dat',status='unknown') 

      do j=0,ngz
        do i=0,ngr
          pn=id_rzg(i,j)
          rp=x_rzg(pn)
          if(i.eq.0) then
            ii=i
            io=i+1
            pnpi=id_rzg(i,j)
            pnpo=id_rzg(i+1,j)
          else if(i.eq.ngr) then
            ii=i-1
            io=i
            pnpi=id_rzg(i-1,j)
            pnpo=id_rzg(i,j)
          else
            ii=i-1
            io=i+1
            pnpi=id_rzg(i-1,j)
            pnpo=id_rzg(i+1,j)
          endif
          dr=x_rzg(pnpo)-x_rzg(pnpi)
          drhovr=(rho_rzg(pnpo)*vr2(io,j)-rho_rzg(pnpi)*vr2(ii,j))/dr
c <vph^2>
          vph2(i,j)=vr2(i,j)+(rp/rho_rzg(pn))*drhovr
     &      +vc2_rzg(pn)
          if(vph2(i,j).lt.0.0d0) then
c            write(6,*) ' vph2 < 0 i,j=',vph2(i,j),i,j
            vph2(i,j)=0.0d0
          endif
c vsig_ph
          vsigph_rzg(pn)=dsqrt(vph2(i,j))
c <vph>^2=0
          vph_rzg(pn)=0.0d0
c
c          if(j.eq.0) then
          write(60,'(13(1pE13.5))') x_rzg(pn),z_rzg(pn)
     &      ,dabs(drhovr),vc2_rzg(pn),rho_rzg(pn),vr2(i,j)
     &      ,vph2(i,j)
     &      ,rho_rzg(pn),(vr2(io,j)-vr2(ii,j))/dr
     &      ,(rho_rzg(pnpo)-rho_rzg(pnpi))/dr,dr
     &      ,(rp/rho_rzg(pn))*drhovr
     &      ,(rp/rho_rzg(pn))
c          endif
c
        enddo
c *** r=0, vph component is zero ***
c        i=0
c        pn=id_rzg(i,j)
c        vph2(i,j)=0.0d0
c        vph_rzg(pn)=0.0d0
c        vsigph_rzg(pn)=0.0d0
      enddo
c      close(60)

      if(myrank.eq.0) then
      open(60,file='hgvgridz0.dat',status='unknown')
      open(61,file='hgvgrid.dat',status='unknown')
      npg=(ngr+1)*(ngz+1)
c      j=0
      do j=0,ngz
        do i=0,ngr
          pn=id_rzg(i,j)
          if(j.eq.0) then
          write(60,'(6(1pE13.5))') x_rzg(pn)
     $         ,VUKMS*dsqrt(vc2_rzg(pn)),vph_rzg(pn),vsigr_rzg(pn)
     $         ,vsigph_rzg(pn),vsigz_rzg(pn)
         endif
          write(61,'(10(1pE13.5))') x_rzg(pn),z_rzg(pn)
     &         ,VUKMS*dsqrt(vc2_rzg(pn)),vph_rzg(pn),vsigr_rzg(pn)
     &         ,vsigph_rzg(pn),vsigz_rzg(pn),rho_rzg(pn)
     &         ,dsqrt(x_rzg(pn)**2+z_rzg(pn)**2)*LUKPC
     &         ,(vsigr_rzg(pn)*x_rzg(pn)+vsigz_rzg(pn)*z_rzg(pn))
     &          /dsqrt(x_rzg(pn)**2+z_rzg(pn)**2)
        enddo
      enddo
      close(60)
      close(61)
      endif

      return
      end
