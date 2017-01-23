c ***********************************************
c      treef.f copy from pv32.17 and modified
c  03 Feb. 2011    written by D.KAWATA
c ***********************************************
c70*******************************************************************
c *** Definition of treeforce() ***
c * calculate acceraration *

      subroutine treeforce(npg)
      include 'define.f'
      integer npg
      integer i,ip,level,is
      integer isend,npj
c * Number of Notfinished particle, temp *      
      integer nlist,tnlist
c * Node searching for each Particle *      
      integer nd,pn,pnj
      double precision r2,ir2,ir3,ir4,ir5,ir7,rij,ir1
c * for calculate gradients *
      double precision dwi,dwj
c * for test far away *            
      double precision theta2
      double precision l_r2
      double precision xij,yij,zij
      double precision dphiir,dphijr,si,sj,hsi,hsj
      double precision d2phii,d2phij,d3phii,d3phij
      double precision dpijr,dp_ij,d2p_ij,d3p_ij
      double precision x2ij,y2ij,z2ij,x3ij,y3ij,z3ij
c *** potential and correction adaptive softening ***
      double precision pot,dhc
c *** for check
      integer ncell,nnode

      ncell=0
      nnode=0

c *** Basic Value ***
      theta2 = THETA*THETA

c *** initialization for a?_p ***
      do i=0,npg-1
        fx_rzg(i)=0.0d0
        fy_rzg(i)=0.0d0
        fz_rzg(i)=0.0d0
      enddo

      npj=npg
c *** initialization ***
      do i=0,npj-1
        list(i)=i
        node(i)=0 
      enddo
      nlist=npj
      level=0
c *** calculate contribution from SPH Particle ***
   77 if(nlist.le.0) then
        goto 93
      endif
      level = level+1
      if(level.gt.MAXNODE) then
        write(6,*) ' Error in treef(): failure in walking tree'
        stop
      endif
      do i=0,nlist-1
        pn=list(i)         
        nd = node(pn)
c *** distance from the center of the mass ***
        xij=x_rzg(pn)-cmx_tr(nd)
        yij=y_rzg(pn)-cmy_tr(nd)
        zij=z_rzg(pn)-cmz_tr(nd)
        r2=xij*xij+yij*yij+zij*zij
        if(r2.gt.0.0d0) then
          ir2 = 1.0d0/r2
          l_r2=(l_tr(nd)+delta_tr(nd))
     &      *(l_tr(nd)+delta_tr(nd))*ir2
        else
          l_r2 = INF
          ir2 = 0.0d0
        endif
c * preriminary task for softened gravity *
c * check well separated *
        if(l_r2.eq.INF.and.np_tr(nd).eq.1) then
          node(pn)=next_tr(nd)
c * calculate softened gravity *
        else if(np_tr(nd).eq.1) then
          pnj = pn_tr(nd)
          rij=dsqrt(r2)
c *** dphi/dr(hj) ***
          hsj=h_p(pnj)
c            hsj=rij
          sj=rij/hsj
          if(sj.lt.1.0d0) then
            is=int(sj/ds_tb)
            if(is.lt.0) then
              is=0
            else if(is.ge.NKTAB) then
              is=NKTAB-1
            endif
c dphi/dr/r
            dphijr=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is))
     &        *(sj-s_tb(is))/ds_tb
            dphijr=dphijr/(hsj**3)
c dw/dr/r
            dwj=dwds_s_tb(is)+(dwds_s_tb(is+1)-dwds_s_tb(is))
     &        *(sj-s_tb(is))/ds_tb
            dwj=dwj/(hsj**5)
          else
            dphijr=ir2*dsqrt(ir2)
c            dphijr=1.0d0/(1.0d0/ir2+h_p(pnj)**2)**(3.0d0/2.0d0)
            dwj=0.0d0
          endif
c *** use hj
          hsi=hsj
          dphiir=dphijr
          dwi=dwj
c *** force calculation ***
c *** dvx_p,dvy_p,dvz_p ***
          pot=m_p(pnj)*G*0.5d0*(dphiir+dphijr)
          fx_rzg(pn)=fx_rzg(pn)-xij*pot
          fy_rzg(pn)=fy_rzg(pn)-yij*pot
          fz_rzg(pn)=fz_rzg(pn)-zij*pot
c * update node *
          node(pn)=next_tr(nd)
        else if(l_r2.le.theta2) then
          ir1=dsqrt(ir2)
          ir3=ir2*ir1
c *** dphi/dr(hi) ***
c use maximum h for node
          hsi=hm_tr(nd)
          rij=dsqrt(r2)
          si=rij/hsi
          if(si.lt.1.0d0) then
            is=int(si/ds_tb)
            if(is.lt.0) then
              is=0
            else if(is.ge.NKTAB) then
              is=NKTAB-1
            endif
c dphi/dr/r
            dphiir=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is))
     &        *(si-s_tb(is))/ds_tb
            dphiir=dphiir/(hsi**3)
c d2phi/dr2(hi)
            d2phii=d2phidr2_tb(is)+(d2phidr2_tb(is+1)-d2phidr2_tb(is))
     &        *(si-s_tb(is))/ds_tb
            d2phii=d2phii/(hsi**3)
c d3phi/dr3(hi)
            d3phii=d3phidr3_tb(is)+(d3phidr3_tb(is+1)-d3phidr3_tb(is))
     &        *(si-s_tb(is))/ds_tb
            d3phii=d3phii/(hsi**4)
          else
            dphiir=ir3
c d2phi/dr2(hi)
            d2phii=-2.0d0*ir3
c d3phi/dr3(hi)
            d3phii=6.0d0*ir2*ir2
          endif
c *** dphi/dr(hj) ***
          hsj=hm_tr(nd)
          sj=rij/hsj
          if(sj.lt.1.0d0) then
            is=int(sj/ds_tb)
            if(is.lt.0) then
              is=0
            else if(is.ge.NKTAB) then
              is=NKTAB-1
            endif
c dphi/dr/r
            dphijr=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is))
     &        *(sj-s_tb(is))/ds_tb
            dphijr=dphijr/(hsj**3)
c d2phi/dr2(hj)
            d2phij=d2phidr2_tb(is)+(d2phidr2_tb(is+1)-d2phidr2_tb(is))
     &        *(sj-s_tb(is))/ds_tb
            d2phij=d2phij/(hsj**3)
c d3phi/dr3(hj)
            d3phij=d3phidr3_tb(is)+(d3phidr3_tb(is+1)-d3phidr3_tb(is))
     &        *(sj-s_tb(is))/ds_tb
            d3phij=d3phij/(hsj**4)
          else
            dphijr=ir3
c d2phj/dr2(hj)
            d2phij=-2.0d0*ir3
c d3phj/dr3(hj)
            d3phij=6.0d0*ir2*ir2
          endif
c *** force calculation ***
c *** dvx_p,dvy_p,dvz_p ***
          dpijr=0.5d0*(dphiir+dphijr)
          pot=mass_tr(nd)*G*dpijr
          dp_ij=dpijr*rij
          d2p_ij=0.5d0*(d2phii+d2phij)
          d3p_ij=0.5d0*(d3phii+d3phij)
          ir4=ir2**2
          ir5=ir3*ir2
          ir7=ir5*ir2
          x2ij=xij**2
          y2ij=yij**2
          z2ij=zij**2
          x3ij=xij*x2ij
          y3ij=yij*y2ij
          z3ij=zij*z2ij
c *** ax ***
          fx_rzg(pn)=fx_rzg(pn)-(xij*pot
c dipole
     &       -mx_tr(nd)*((ir1-x2ij*ir3)*dp_ij
     &        +x2ij*ir2*d2p_ij)
     &       -my_tr(nd)*(-xij*yij*ir3*dp_ij+xij*yij*ir2*d2p_ij)
     &       -mz_tr(nd)*(-xij*zij*ir3*dp_ij+xij*zij*ir2*d2p_ij)
c quadrupole
     &       +0.5d0*mxx_tr(nd)
     &        *((-3.0d0*xij*ir3+3.0d0*x3ij*ir5)*dp_ij
     &         +(3.0d0*xij*ir2-3.0d0*x3ij*ir4)*d2p_ij
     &         +x3ij*ir3*d3p_ij)
c
     &       +0.5d0*myy_tr(nd)
     &        *((-xij*ir3+3.0d0*xij*y2ij*ir5)*dp_ij
     &         +(xij*ir2-3.0d0*xij*y2ij*ir4)*d2p_ij
     &          +xij*y2ij*ir3*d3p_ij)
c
     &       +0.5d0*mzz_tr(nd)
     &        *((-xij*ir3+3.0d0*xij*z2ij*ir5)*dp_ij
     &         +(xij*ir2-3.0d0*xij*z2ij*ir4)*d2p_ij
     &         +xij*z2ij*ir3*d3p_ij)
c
     &       +mxy_tr(nd)
     &        *((-yij*ir3+3.0d0*x2ij*yij*ir5)*dp_ij
     &         +(yij*ir2-3.0d0*x2ij*yij*ir4)*d2p_ij
     &          +x2ij*yij*ir3*d3p_ij)
c
     &       +mzx_tr(nd)
     &        *((-zij*ir3+3.0d0*x2ij*zij*ir5)*dp_ij
     &         +(zij*ir2-3.0d0*x2ij*zij*ir4)*d2p_ij
     &          +x2ij*zij*ir3*d3p_ij)
c
     &       +myz_tr(nd)
     &        *(3.0d0*xij*yij*zij*ir5*dp_ij
     &         +3.0d0*xij*yij*zij*ir4*d2p_ij
     &         +xij*yij*zij*ir3*d3p_ij))
c *** ay ***
          fy_rzg(pn)=fy_rzg(pn)-(yij*pot
c dipole
     &       -mx_tr(nd)*(-xij*yij*ir3*dp_ij+xij*yij*ir2*d2p_ij)
     &       -my_tr(nd)*((ir1-y2ij*ir3)*dp_ij
     &        +y2ij*ir2*d2p_ij)
     &       -mz_tr(nd)*(-yij*zij*ir3*dp_ij+yij*zij*ir2*d2p_ij)
c quadrupole
     &       +0.5d0*mxx_tr(nd)
     &        *((-yij*ir3+3.0d0*yij*x2ij*ir5)*dp_ij
     &         +(yij*ir2-3.0d0*yij*x2ij*ir4)*d2p_ij
     &          +yij*x2ij*ir3*d3p_ij)
c
     &       +0.5d0*myy_tr(nd)
     &        *((-3.0d0*yij*ir3+3.0d0*y3ij*ir5)*dp_ij
     &         +(3.0d0*yij*ir2-3.0d0*y3ij*ir4)*d2p_ij
     &         +y3ij*ir3*d3p_ij)
c
     &       +0.5d0*mzz_tr(nd)
     &        *((-yij*ir3+3.0d0*yij*z2ij*ir5)*dp_ij
     &         +(yij*ir2-3.0d0*yij*z2ij*ir4)*d2p_ij
     &          +yij*z2ij*ir3*d3p_ij)
c
     &       +mxy_tr(nd)
     &        *((-xij*ir3+3.0d0*y2ij*xij*ir5)*dp_ij
     &         +(xij*ir2-3.0d0*y2ij*xij*ir4)*d2p_ij
     &          +y2ij*xij*ir3*d3p_ij)
c
     &       +mzx_tr(nd)
     &        *(3.0d0*xij*yij*zij*ir5*dp_ij
     &         +3.0d0*xij*yij*zij*ir4*d2p_ij
     &         +xij*yij*zij*ir3*d3p_ij))
c
     &       +myz_tr(nd)
     &        *((-zij*ir3+3.0d0*y2ij*zij*ir5)*dp_ij
     &         +(zij*ir2-3.0d0*y2ij*zij*ir4)*d2p_ij
     &          +y2ij*zij*ir3*d3p_ij)
*** az ***
          fz_rzg(pn)=fz_rzg(pn)-(zij*pot
c dipole
     &       -mx_tr(nd)*(-xij*zij*ir3*dp_ij+xij*zij*ir2*d2p_ij)
     &       -my_tr(nd)*(-yij*zij*ir3*dp_ij+yij*zij*ir2*d2p_ij)
     &       -mz_tr(nd)*((ir1-z2ij*ir3)*dp_ij
     &        +z2ij*ir2*d2p_ij)
c quadrupole
     &       +0.5d0*mxx_tr(nd)
     &        *((-zij*ir3+3.0d0*zij*x2ij*ir5)*dp_ij
     &         +(zij*ir2-3.0d0*zij*x2ij*ir4)*d2p_ij
     &          +zij*x2ij*ir3*d3p_ij)
c
     &       +0.5d0*myy_tr(nd)
     &        *((-zij*ir3+3.0d0*zij*y2ij*ir5)*dp_ij
     &         +(zij*ir2-3.0d0*zij*y2ij*ir4)*d2p_ij
     &          +zij*y2ij*ir3*d3p_ij)
c
     &       +0.5d0*mzz_tr(nd)
     &        *((-3.0d0*zij*ir3+3.0d0*z3ij*ir5)*dp_ij
     &         +(3.0d0*zij*ir2-3.0d0*z3ij*ir4)*d2p_ij
     &         +z3ij*ir3*d3p_ij)
c
     &       +mxy_tr(nd)
     &        *(3.0d0*xij*yij*zij*ir5*dp_ij
     &         +3.0d0*xij*yij*zij*ir4*d2p_ij
     &         +xij*yij*zij*ir3*d3p_ij))
c
     &       +mzx_tr(nd)
     &        *((-xij*ir3+3.0d0*z2ij*xij*ir5)*dp_ij
     &         +(xij*ir2-3.0d0*z2ij*xij*ir4)*d2p_ij
     &          +z2ij*xij*ir3*d3p_ij)
c
     &       +myz_tr(nd)
     &        *((-yij*ir3+3.0d0*z2ij*yij*ir5)*dp_ij
     &         +(yij*ir2-3.0d0*z2ij*yij*ir4)*d2p_ij
     &          +z2ij*yij*ir3*d3p_ij)
c *** update node
          node(pn)=next_tr(nd)
        else
          node(pn) = daughter_tr(nd)
        endif
      enddo
c * update not-finished particle list *
      tnlist = nlist
      nlist = 0
!ocl novrec (list)
!cdir nodep (list)
      do i=0,tnlist-1
        if(node(list(i)).ne.0) then
          list(nlist)=list(i)
          nlist=nlist+1
        endif
      enddo
      goto 77      
c *** sum up the force
   93 do i=0,npg-1
        tx(i)=0.0d0
      enddo

      call MPI_ALLREDUCE(fx_rzg,tx,npg,MPI_DOUBLE_PRECISION
     &  ,MPI_SUM,MPI_COMM_WORLD,ierr)
      do i=0,npg-1
        fx_rzg(i)=tx(i)
      enddo
      do i=0,npg-1
        tx(i)=0.0d0
      enddo
      call MPI_ALLREDUCE(fy_rzg,tx,npg,MPI_DOUBLE_PRECISION
     &  ,MPI_SUM,MPI_COMM_WORLD,ierr)
      do i=0,npg-1
        fy_rzg(i)=tx(i)
      enddo
      do i=0,npg-1
        tx(i)=0.0d0
      enddo
      call MPI_ALLREDUCE(fz_rzg,tx,npg,MPI_DOUBLE_PRECISION
     &  ,MPI_SUM,MPI_COMM_WORLD,ierr)
      do i=0,npg-1
        fz_rzg(i)=tx(i)
      enddo

c comparison with direct summation
c *** output force field ***
c      open(60,file='fgrid.dat',status='unknown')
c      do i=0,npg-1
c        write(60,'(6(1pE13.5))') x_rzg(i),y_rzg(i),z_rzg(i)
c     &   ,fx_rzg(i),fy_rzg(i),fz_rzg(i)
c      enddo
c      close(60)
c
c      do i=0,npg-1
c        fx_rzg(i)=0.0d0
c        fy_rzg(i)=0.0d0
c        fz_rzg(i)=0.0d0
c        do j=0,np-1
c          xij=x_p(j)-x_rzg(i)
c          yij=y_p(j)-y_rzg(i)
c          zij=z_p(j)-z_rzg(i)
c          rij=xij**2+yij**2+zij**2
c          if(dsqrt(rij).gt.h_p(j)) then
c            pot=m_p(j)/((rij)**(3.0d0/2.0d0))
c          else
c            pot=m_p(j)/((rij+(h_p(j)/4.0d0)**2)**(3.0d0/2.0d0))
c          endif
c          fx_rzg(i)=fx_rzg(i)+pot*xij
c          fy_rzg(i)=fy_rzg(i)+pot*yij
c          fz_rzg(i)=fz_rzg(i)+pot*zij
c        enddo
c      enddo
c      open(60,file='dfgrid.dat',status='unknown')
c      do i=0,npg-1
c        write(60,'(6(1pE13.5))') x_rzg(i),y_rzg(i),z_rzg(i)
c     &   ,fx_rzg(i),fy_rzg(i),fz_rzg(i)
c      enddo
c      close(60)

      return
      end

