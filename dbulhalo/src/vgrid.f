c *****************************************************
c    vgrid.f
c 10 Jun., 2014   written by D. Kawata
c ***************************************************** 
c70*******************************************************************

      subroutine vgrid(ngr,ngz,fr,mdisk,hd,id)
      include './define.f'
      integer ngr,ngz,npg,id
      double precision mdisk,hd
      integer i,j,k,pn,pnpi,pnpo,ii,io
      double precision fin,dz,fp,fpp,dr,rp
      double precision vz2(0:MN2D-1,0:MN2D-1),vr2(0:MN2D-1,0:MN2D-1)
     &  ,vph2(0:MN2D-1,0:MN2D-1),eta2(0:MN2D-1,0:MN2D-1)
     &  ,kapa2(0:MN2D-1,0:MN2D-1)
      double precision dph2dr2(0:MN2D-1,0:MN2D-1)
     &  ,drhovr(0:MN2D-1,0:MN2D-1)
      double precision fr,qval,sdend
      character filen*60

c *** sample eta
c      open(60,file='eta.dat',status='unknown') 

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
          if(vz2(i,j).gt.0.0d0) then
            vsigz_rzg(pn)=dsqrt(vz2(i,j))
          else
            vsigz_rzg(pn)=0.0d0
          endif
c *** set vsigr2 ***
          vr2(i,j)=fr*vz2(i,j)
          vsigr_rzg(pn)=dsqrt(vr2(i,j))
c          if(i.gt.0) then
c *** circular velocity ***
            if(fx_rzg(pn).lt.0.0d0) then
              vc2_rzg(pn)=-x_rzg(pn)*fx_rzg(pn)
            else
              vc2_rzg(pn)=0.0d0
            endif
c *** get eta ***
            if(i.le.1) then
              pnpi=id_rzg(i,j)
              pnpo=id_rzg(i+1,j)
! d dln A/d log10 r
              if((fx_rzg(pnpi).lt.0.0d0.and.fx_rzg(pnpo).lt.0.0d0)
     &          .and.fx_rzg(id_rzg(i+2,j)).lt.0.0d0) then
                dph2dr2(i,j)=(-dlog(-fx_rzg(id_rzg(i+2,j)))
     &           +4.0d0*dlog(-fx_rzg(pnpo))-3.0d0*dlog(-fx_rzg(pnpi)))
     &           /(2.0d0*dlr_rzg)
              else
                dph2dr2(i,j)=0.0d0
              endif
            else if(i.eq.ngr) then
              pnpi=id_rzg(i-1,j)
              pnpo=id_rzg(i,j)
              if((fx_rzg(pnpi).lt.0.0d0.and.fx_rzg(pnpo).lt.0.0d0)
     &          .and.fx_rzg(id_rzg(i-2,j)).lt.0.0d0) then
                dph2dr2(i,j)=(3.0d0*dlog(-fx_rzg(pnpo))
     &           -4.0d0*dlog(-fx_rzg(pnpi))
     &           +dlog(-fx_rzg(id_rzg(i-2,j))))
     &           /(2.0d0*dlr_rzg) 
              else
                dph2dr2(i,j)=0.0d0
              endif
            else
              pnpi=id_rzg(i-1,j)
              pnpo=id_rzg(i+1,j)
              if(fx_rzg(pnpi).lt.0.0d0.and.fx_rzg(pnpo).lt.0.0d0) then
! d dln A/d log10 r
                dph2dr2(i,j)=(dlog(-fx_rzg(pnpo))-dlog(-fx_rzg(pnpi)))
     &           /(2.0d0*dlr_rzg)
              else
                dph2dr2(i,j)=0.0d0
              endif
            endif
!            dr=x_rzg(pnpo)-x_rzg(pnpi)
!            dph2dr2(i,j)=(-fx_rzg(pnpo)+fx_rzg(pnpi))/dr
!
            if(fx_rzg(pn).lt.0.0d0.and.rp.gt.0.0d0) then
              dph2dr2(i,j)=-fx_rzg(pn)*(1.0d0/(rp*dlog(10.0d0)))
     &         *dph2dr2(i,j)
            else
              dph2dr2(i,j)=0.0d0
            endif
            if(fx_rzg(pn).lt.0.0d0.and.rp.gt.0.0d0) then
              eta2(i,j)=(4.0d0/rp)*(-fx_rzg(pn))
     &          /((3.0d0/rp)*(-fx_rzg(pn))+dph2dr2(i,j))
              kapa2(i,j)=((3.0d0/rp)*(-fx_rzg(pn))+dph2dr2(i,j))
            else
              eta2(i,j)=0.0d0
              kapa2(i,j)=0.0d0
            endif
            if(eta2(i,j).lt.0.0d0) then
              eta2(i,j)=0.0d0
              kapa2(i,j)=0.0d0
            endif
c          endif
        enddo
      enddo

c *** get sigma phi
c      open(60,file='vph.dat',status='unknown') 

      do j=0,ngz
        do i=0,ngr
          pn=id_rzg(i,j)
          rp=x_rzg(pn)
          drhovr(i,j)=0.0d0
c *** r=0, vph component is zero ***
          if(i.eq.0) then
            vph2(i,j)=0.0d0
            vph_rzg(pn)=0.0d0
            vsigph_rzg(pn)=0.0d0
          else
            if(i.eq.1) then
              ii=i
              io=i+1
              pnpi=id_rzg(i,j)
              pnpo=id_rzg(i+1,j)

c! d dln A/d log10 r
              if((rho_rzg(pnpo)*vr2(io,j).gt.0.0d0
     &          .and.rho_rzg(pnpi)*vr2(ii,j).gt.0.0d0)
     &          .and.rho_rzg(id_rzg(i+2,j))*vr2(i+2,j).gt.0.0d0) then
                drhovr(i,j)=(-dlog(rho_rzg(id_rzg(i+2,j))*vr2(i+2,j))
     &           +4.0d0*dlog(rho_rzg(pnpo)*vr2(io,j))
     &           -3.0d0*dlog(rho_rzg(pnpi)*vr2(ii,j)))/(2.0d0*dlr_rzg)
              else
                drhovr(i,j)=0.0d0
              endif
            else if(i.eq.ngr) then
              ii=i-1
              io=i
              pnpi=id_rzg(i-1,j)
              pnpo=id_rzg(i,j)
              if((rho_rzg(pnpo)*vr2(io,j).gt.0.0d0
     &          .and.rho_rzg(pnpi)*vr2(ii,j).gt.0.0d0)
     &          .and.rho_rzg(id_rzg(i-2,j))*vr2(i-2,j).gt.0.0d0) then
                drhovr(i,j)=(3.0d0*dlog(rho_rzg(pnpo)*vr2(io,j))
     &            -4.0d0*dlog(rho_rzg(pnpi)*vr2(ii,j))
     &            +dlog(rho_rzg(id_rzg(i-2,j))*vr2(i-2,j)))
     &            /(2.0d0*dlr_rzg)
              else
                drhovr(i,j)=0.0d0
              endif
            else
              ii=i-1
              io=i+1
              pnpi=id_rzg(i-1,j)
              pnpo=id_rzg(i+1,j)
              if((rho_rzg(pnpo)*vr2(io,j).gt.0.0d0
     &          .and.rho_rzg(pnpi)*vr2(ii,j).gt.0.0d0)) then
                drhovr(i,j)=(dlog(rho_rzg(pnpo)*vr2(io,j))
     &            -dlog(rho_rzg(pnpi)*vr2(ii,j)))/(2.0d0*dlr_rzg)
              else
                drhovr(i,j)=0.0d0
              endif
            endif
! dA/dR
            if(rp.gt.0.0d0) then
              drhovr(i,j)=rho_rzg(id_rzg(i,j))*vr2(i,j)*drhovr(i,j)
     &          /(rp*dlog(10.0d0)) 
            else
              drhovr(i,j)=0.0d0
            endif
c <vph^2>
            vph2(i,j)=vr2(i,j)+(rp/rho_rzg(pn))*drhovr(i,j)
     &        +vc2_rzg(pn)
c not allow to vph>vc
            if(vph2(i,j).gt.vc2_rzg(pn)) then
              vph2(i,j)=vc2_rzg(pn)
            endif

           if(i.eq.1.and.j.eq.100) then
             write(6,*) vr2(i,j),(rp/rho_rzg(pn)),drhovr(i,j)
     &         ,vc2_rzg(pn),vph2(i,j)
           endif


            if(vph2(i,j).lt.0.0d0) then
c            write(6,*) ' vph2 < 0 i,j=',vph2(i,j),i,j
              vph2(i,j)=0.0d0
            endif
c vsig_ph
c          if(vph2(i,j).lt.vph_rzg(pn)) then
c            write(6,*) ' Vph error <vph2>,<vph>2,r='
c     &       ,vph2(i,j),vph_rzg(pn),x_rzg(pn)
c          endif

            if(eta2(i,j).lt.1.0d0) then
              eta2(i,j)=1.0d0
            endif
            if(eta2(i,j).gt.0.0d0) then
              vsigph_rzg(pn)=dsqrt(vr2(i,j)/eta2(i,j))
c <vph>^2
              vph_rzg(pn)=vph2(i,j)-vr2(i,j)/eta2(i,j)
            else
              vsigph_rzg(pn)=0.0d0
              vph_rzg(pn)=vph2(i,j)
            endif
            if(vph_rzg(pn).gt.0.0d0) then
              vph_rzg(pn)=dsqrt(vph_rzg(pn))
            else
               vph_rzg(pn)=0.0d0
            endif
          endif
c
c          if(j.eq.0) then
c          write(60,'(13(1pE13.5))') x_rzg(pn),z_rzg(pn)
c     &      ,dabs(drhovr(i,j)),vc2_rzg(pn),rho_rzg(pn),vr2(i,j)
c     &      ,vph2(i,j)
c     &      ,rho_rzg(pnpo),vr2(io,j),rho_rzg(pnpi),vr2(ii,j),dr
c          endif
c
        enddo
      enddo
c      close(60)
c set i=0 values, except vsigz, vsigr
c      i=0
c      do j=0,ngz
c        pn=id_rzg(i,j)
c        vc2_rzg(pn)=0.0d0
c        vph_rzg(pn)=0.0d0
c      enddo


      if(myrank.eq.0) then
      write(filen,'(a8,i3.3,a4)') 'vgridz0c',id,'.dat'
      open(60,file=filen,status='unknown')
      write(filen,'(a6,i3.3,a4)') 'vgridc',id,'.dat'
      open(61,file=filen,status='unknown')
      npg=(ngr+1)*(ngz+1)
c      j=0
      do j=0,ngz
        do i=0,ngr
          pn=id_rzg(i,j)
          if(j.eq.0) then
          sdend=(mdisk/(2.0d0*M_PI*(hd**2)))
     &      *dexp(-x_rzg(pn)/hd)
          qval=vsigr_rzg(pn)*dsqrt(kapa2(i,j))
     &     /(3.36d0*G*sdend)
          write(60,'(17(1pE13.5))') x_rzg(pn),y_rzg(pn),z_rzg(pn)
     $         ,VUKMS*dsqrt(vc2_rzg(pn)),vph_rzg(pn),vsigr_rzg(pn)
     $         ,vsigph_rzg(pn),vsigz_rzg(pn),eta2(i,j) ,vph2(i,j)
     $         ,dph2dr2(i,j),fx_rzg(pn),rho_rzg(pn) ,vr2(i,j)
     $         ,qval,-fx_rzg(pn),drhovr(i,j)
         endif
          write(61,'(17(1pE13.5))') x_rzg(pn),y_rzg(pn),z_rzg(pn)
     &         ,VUKMS*dsqrt(vc2_rzg(pn)),vph_rzg(pn),vsigr_rzg(pn)
     &         ,vsigph_rzg(pn),vsigz_rzg(pn),eta2(i,j) ,vph2(i,j)
     &         ,dph2dr2(i,j),fx_rzg(pn),rho_rzg(pn) ,vr2(i,j)
     &         ,qval,rho_rzg(pn),drhovr(i,j)
        enddo
      enddo
      close(60)
      close(61)
      endif

      return
      end

