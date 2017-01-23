c *********************************************
c  setcool.f copy for GCD+ pv.33.6
c  8 Feb. 2011   written by D. Kawata
c *********************************************
c70*******************************************************************

      subroutine  setcool()
      include 'define.f'
      integer MNVALP
      parameter (MNVALP=1000000)
      integer i,j,k
      integer iz,inh,imet,it,nval,nc
      integer lmet10,dmet10,lnh10,dnh10
      integer z100,met10,nh10
      character czred*5,cmet*5,cnh*5,filen*60
      double precision crad,hrad
      integer tivr(0:10)
      double precision tdvr(0:MNVALP)

      if(myrank.eq.0) then

        open(50,file='./cool/setcool.dat',status='old')
        read(50,'(I6)') nz_crtb
        read(50,150) nnh_crtb,lnh10,dnh10
        read(50,150) nmet_crtb,lmet10,dmet10
        read(50,'(I6)') nt_crtb
c        read(50,'(F5.2)') SI_zeor
  150   format(3I6)
        close(50)
      
        open(51,file='./cool/setcoolzred.dat',status='old')
        do iz=0,nz_crtb-1
          read(51,*) z100
          z_crtb(iz)=dble(z100)/100.0d0
          write(czred,'(a1,i4.4)') 'z',z100
          do inh=0,nnh_crtb-1
            nh10=lnh10+inh*dnh10
            nh_crtb(inh)=dble(nh10)/10.0d0
            if(nh10.ge.0) then
              write(cnh,'(a2,i3.3)') 'hp',nh10
            else
              write(cnh,'(a2,i3.3)') 'hm',-nh10
            endif
            do imet=0,nmet_crtb-1
              met10=lmet10+imet*dmet10
              met_crtb(imet)=dble(met10)/10.0d0
              if(met10.ge.0) then
                write(cmet,'(a2,i3.3)') 'Zp',met10
              else
                write(cmet,'(a2,i3.3)') 'Zm',-met10
              endif
c *** set filename ***
              write(filen,'(a13,a5,a5,a5,a4)')
     &          './cool/chrate',czred,cmet,cnh,'.dat'
              open(50,file=filen,status='old')
              do i=1,4
                read(50,*)
              enddo
              do i=0,nt_crtb-1
                read(50,151) t_crtb(i),cr_crtb(i,imet,inh,iz)
     &           ,hr_crtb(i,imet,inh,iz)
     &           ,myu_crtb(i,imet,inh,iz)
  151           format(f5.3,f8.3,f8.3,f6.3) 
c
c                write(60,'(7(1pE13.5))') t_crtb(i)
c     &            ,cr_crtb(i,imet,inh,iz)
c     &            ,hr_crtb(i,imet,inh,iz)
c     &            ,myu_crtb(i,imet,inh,iz)
c     &            ,met_crtb(imet),nh_crtb(inh),z_crtb(iz)
c
              enddo
              close(50)
c              write(60,*)
            enddo
c            write(60,*)
          enddo
c          write(60,*)
        enddo
        close(51)
c        close(60)
c *** sending numbers to the other procs ***
        tivr(0)=nz_crtb
        tivr(1)=nnh_crtb
        tivr(2)=nmet_crtb
        tivr(3)=nt_crtb
        write(6,*) ' cooling table highest z=',z_crtb(nz_crtb-1)
c        write(6,*) ' EoR z=',SI_zeor
      endif
      call MPI_BCAST(tivr,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        nz_crtb=tivr(0)
        nnh_crtb=tivr(1)
        nmet_crtb=tivr(2)
        nt_crtb=tivr(3)
      endif
c *** send zeor ***
c      call MPI_BCAST(SI_zeor,1,MPI_DOUBLE_PRECISION,0
c     &  ,MPI_COMM_WORLD,ierr)
c *** send zred grid ***
      nval=nz_crtb
      if(myrank.eq.0) then
        do iz=0,nz_crtb-1
          tdvr(iz)=z_crtb(iz)
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0
     &  ,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        do iz=0,nz_crtb-1
          z_crtb(iz)=tdvr(iz)
        enddo
      endif
c *** send nh grid ***
      nval=nnh_crtb
      if(myrank.eq.0) then
        do inh=0,nnh_crtb-1
          tdvr(inh)=nh_crtb(inh)
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0
     &  ,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        do inh=0,nnh_crtb-1
          nh_crtb(inh)=tdvr(inh)
        enddo
      endif
c *** send metal grid ***
      nval=nmet_crtb
      if(myrank.eq.0) then
        do imet=0,nmet_crtb-1
          tdvr(imet)=met_crtb(imet)
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0
     &  ,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        do imet=0,nmet_crtb-1
          met_crtb(imet)=tdvr(imet)
        enddo
      endif
c *** send T grid ***
      nval=nt_crtb
      if(myrank.eq.0) then
        do i=0,nt_crtb-1
          tdvr(i)=t_crtb(i)
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0
     &  ,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        do i=0,nt_crtb-1
          t_crtb(i)=tdvr(i)
        enddo
      endif
c *** send cooling rate ***
      nval=nt_crtb*nmet_crtb*nnh_crtb*nz_crtb
      if(myrank.eq.0) then
        if(MNVALP.lt.nval) then
          write(6,*) ' Error in setcool(): too big cooling array.'
          write(6,*) ' MNVALP,nval=',MNVALP,nval
          stop
        endif
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                tdvr(nc)=cr_crtb(i,imet,inh,iz)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0
     &  ,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                cr_crtb(i,imet,inh,iz)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif
c *** send heating rate ***
      if(myrank.eq.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                tdvr(nc)=hr_crtb(i,imet,inh,iz)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0
     &  ,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                hr_crtb(i,imet,inh,iz)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif
c *** send myu ***
      if(myrank.eq.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                tdvr(nc)=myu_crtb(i,imet,inh,iz)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0
     &  ,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                myu_crtb(i,imet,inh,iz)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif
c for test
c      if(myrank.eq.1) then
c        write(filen,'(a7,i3.3,a4)') 'setcool',myrank,'.dat'
c        open(60,file=filen,status='unknown')
c        do iz=0,nz_crtb-1
c          do inh=0,nnh_crtb-1
c            do imet=0,nmet_crtb-1
c              do i=0,nt_crtb-1
c                write(60,'(7(1pE13.5))') t_crtb(i)
c     &            ,cr_crtb(i,imet,inh,iz)
c     &            ,hr_crtb(i,imet,inh,iz)
c     &            ,myu_crtb(i,imet,inh,iz)
c     &            ,met_crtb(imet),nh_crtb(inh),z_crtb(iz)
c              enddo
c              write(60,*)
c            enddo
c          enddo
c        enddo
c        close(60)
c      endif
      return
      end

