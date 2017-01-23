c ***********************************************
c    tree.f copy from pv.32.17
c  18 Jan. 2011   written by D. Kawata
c ***********************************************
c70*******************************************************************

c ***********************************************
c    Definition of function about octtree
c  This program build octtree and
c    compute the center of mass and total mass and 
c ***********************************************

      subroutine treebuild(np)
      include 'define.f'
      integer np	  
      integer i,pn,nd,level,npn
c * for not-finished particle (nfp) *
c * number of not-finished particle *	  
      integer numnfp
c * for subcell (sc) *
c * start subcell in tree *	  
      integer stsc
c * number of subcells *	  
      integer numsc,clist(0:MNS-1)
c * diameter of subcell *	  
      double precision l_sc
c * for work *
      integer numtmp
c *** for calculating hm_tr ***
      integer d,num,npare,pnd,j

c *** for Tree ***
      tcx(0) = 1.0d0
      tcx(1) = 1.0d0
      tcx(2) = 1.0d0
      tcx(3) = 1.0d0
      tcx(4) = -1.0d0
      tcx(5) = -1.0d0
      tcx(6) = -1.0d0
      tcx(7) = -1.0d0
      tcy(0) = 1.0d0
      tcy(1) = 1.0d0
      tcy(4) = 1.0d0
      tcy(5) = 1.0d0
      tcy(2) = -1.0d0
      tcy(3) = -1.0d0
      tcy(6) = -1.0d0
      tcy(7) = -1.0d0
      tcz(0) = 1.0d0
      tcz(2) = 1.0d0
      tcz(4) = 1.0d0
      tcz(6) = 1.0d0
      tcz(1) = -1.0d0
      tcz(3) = -1.0d0
      tcz(5) = -1.0d0
      tcz(7) = -1.0d0
	  
c *** Make root ***	
      call makeroot(np)
      next_tr(0)=0
      pare_tr(0)=-1
      num_tr=1	
      if(np.eq.0) then
        np_tr(0)=0
        num_tr=0
        return
      endif

c * initialization *
      numnfp = np

      do i=0,np-1
        pn_nfp(i)=i
        nd_nfp(i)=0
        label_nfp(i)=0
      enddo
      numsc=8
      stsc=1
      l_sc=l_tr(0)*0.5d0	
      do i=0,7
        np_sc(i)=0
        pare_sc(i)=0
        c_sc(i)=i
      enddo
      flag_pc(0)=0
c * in subcell *	  
      daughter_tr(0) = 0
      level = 0
      list(0) = 0
c *****   start iteration in level *****
   77 if(numnfp.le.1) then
        goto 99
      endif
c      write(6,*) numnfp,level,l_sc
      if(level.gt.MAXNODE) then
        write(6,*) ' Error in tree(): failure in treebuild'
        stop
      endif        
c * initialization for the case np_tr < 8 *
      do i=list(level),num_tr-1
        nd_sc(i)=0
      enddo
      level=level+1	  
      list(level) = num_tr
      l_sc = 0.5d0*l_sc
c * find out which subcell it is in *
      do i=0,numnfp-1
        pn=pn_nfp(i)
        nd=nd_nfp(i)
        tx(i)=x_p(pn)-cx_tr(nd)
        ty(i)=y_p(pn)-cy_tr(nd)
        tz(i)=z_p(pn)-cz_tr(nd)
      enddo
      do i=0,numnfp-1
        c_nfp(i)=0
        if(tx(i).lt.0.0d0) then
          c_nfp(i)=c_nfp(i)+4
        endif		  
        if(ty(i).lt.0.0d0) then
          c_nfp(i)=c_nfp(i)+2
        endif		  
        if(tz(i).lt.0.0d0) then
          c_nfp(i)=c_nfp(i)+1
        endif
      enddo
      npn = 0
      do i=0,numnfp-1
        if(np_tr(nd_nfp(i)).gt.8) then
          nd_nfp(i) = daughter_tr(nd_nfp(i))+c_nfp(i)
        else
          nalist(npn)=nd_nfp(i)
          talist(npn)=i
          npn = npn+1
        endif
      enddo                

!ocl novrec (nd_nfp)
!cdir nodep (nd_nfp)
      do i=0,npn-1
        nd_nfp(talist(i))=daughter_tr(nalist(i))+nd_sc(nalist(i))
        nd_sc(nalist(i))=nd_sc(nalist(i))+1
      enddo

      do i=0,npn-1
        c_sc(nd_nfp(talist(i)))=c_nfp(talist(i))  
      enddo

c * update info of subcell *
      do i=0,numnfp-1
        np_sc(nd_nfp(i))=np_sc(nd_nfp(i))+1
      enddo

      if(num_tr.gt.MAXNODE-numnfp) then
        write(6,*) ' Error in buildtree():Node is overflow!'
        write(6,*) ' This level is ',level,'numnfp=',numnfp
        stop
      endif

c *** subcell is conected to tree ***
      npn = 0
      do i=0,numsc-1
        if(np_sc(i).ge.1) then
          nalist(npn)=i
          talist(npn)=num_tr
          node(npn)=pare_sc(i)
          clist(npn)=c_sc(i)
          num_tr=num_tr+1
          npn=npn+1
        endif
      enddo

      do i=0,npn-1
        if(flag_pc(node(i)).eq.0) then
          daughter_tr(node(i))=talist(i)
          flag_pc(node(i))=1
        endif
      enddo
!ocl novrec (cx_tr,cy_tr,cz_tr)
!cdir nodep (cx_tr,cy_tr,cz_tr)
      do i=0,npn-1
        nd_sc(nalist(i))=talist(i)
c *****  Set subcell in tree   *****
        np_tr(talist(i))=np_sc(nalist(i))
        l_tr(talist(i))=2.0d0*l_sc
        pare_tr(talist(i))=node(i)
        cx_tr(talist(i))=tcx(clist(i))*l_sc+cx_tr(node(i))
        cy_tr(talist(i))=tcy(clist(i))*l_sc+cy_tr(node(i))
        cz_tr(talist(i))=tcz(clist(i))*l_sc+cz_tr(node(i))
      enddo		

c *** Set label not-finished particle ***
      do i=0,numnfp-1
        nd_nfp(i)=nd_sc(nd_nfp(i))
        pn_tr(nd_nfp(i)) = pn_nfp(i)
c *  this node is leaf *		  
        if(np_tr(nd_nfp(i)).eq.1) then
          label_nfp(i)=1
        endif			
      enddo
c *** rebuild not finished particle list ***
      numtmp = 0
!ocl novrec (label_nfp,pn_nfp,nd_nfp)
!cdir nodep (label_nfp,pn_nfp,nd_nfp)
      do i=0,numnfp-1
        if(label_nfp(i).eq.0) then
          pn_nfp(numtmp)=pn_nfp(i)
          nd_nfp(numtmp)=nd_nfp(i)
          label_nfp(numtmp)=0
          numtmp=numtmp+1
        endif
      enddo
      numnfp=numtmp
c *** rebuild subcell ***
      numtmp = num_tr-stsc
      numsc = 0
!ocl novrec (np_sc,nd_sc,pare_sc,c_sc)
!cdir nodep (np_sc,nd_sc,pare_sc,c_sc)
      do i=0,numtmp-1
        if(np_tr(stsc).ge.2) then
          daughter_tr(stsc)=numsc
          flag_pc(stsc) = 0
          np_sc(numsc)=0
          nd_sc(numsc)=0
          pare_sc(numsc)=stsc
          c_sc(numsc) = 0
          np_sc(numsc+1)=0
          nd_sc(numsc+1)=0
          pare_sc(numsc+1)=stsc
          c_sc(numsc+1) = 1
          np_sc(numsc+2)=0
          nd_sc(numsc+2)=0
          pare_sc(numsc+2)=stsc
          c_sc(numsc+2) = 2
          np_sc(numsc+3)=0
          nd_sc(numsc+3)=0
          pare_sc(numsc+3)=stsc
          c_sc(numsc+3) = 3
          np_sc(numsc+4)=0
          nd_sc(numsc+4)=0
          pare_sc(numsc+4)=stsc
          c_sc(numsc+4) = 4
          np_sc(numsc+5)=0
          nd_sc(numsc+5)=0
          pare_sc(numsc+5)=stsc
          c_sc(numsc+5) = 5
          np_sc(numsc+6)=0
          nd_sc(numsc+6)=0
          pare_sc(numsc+6)=stsc
          c_sc(numsc+6) = 6
          np_sc(numsc+7)=0
          nd_sc(numsc+7)=0
          pare_sc(numsc+7)=stsc
          c_sc(numsc+7) = 7
          numsc=numsc+8
        else if(np_tr(stsc).eq.1) then
          daughter_tr(stsc)=-1
        endif
        stsc=stsc+1
      enddo
      stsc = num_tr
      if(numsc.gt.MAXNODE) then
        write(6,*) ' Error in buildtree() : '
        write(6,*) '  Subcell is overflow	!'
        write(6,*) '  This level is ',level
        stop
      endif
      goto 77
c *** set next node ***
   99 next_tr(0)=0
      do i=0,num_tr-2
        pare_sc(i)=pare_tr(i+1)
      enddo		
      do i=1,num_tr-2
        if(pare_sc(i).eq.pare_tr(i)) then 
          next_tr(i)=i+1
        else
          next_tr(i) = next_tr(pare_tr(i))
        endif
      enddo
      next_tr(num_tr-1)=next_tr(pare_tr(num_tr-1))

c *** compute mass and center of mass ***
      call compute_mass(level)
      return
      end

c *** Definition of makeroot() ***
      subroutine makeroot(np)
      include 'define.f'
      integer i,np
c * max coordinate *	  
      double precision max_x,max_y,max_z
c * min coordinate <0 *	  
      double precision min_x,min_y,min_z
c * max,temp length *
      double precision maxl,tl

c *** Define root node ***
      max_x=-INF
      max_y=-INF
      max_z=-INF
      min_x=INF
      min_y=INF
      min_z=INF
      do i=0,np-1
        if(x_p(i).lt.min_x) then
          min_x = x_p(i)
        endif		  
        if(y_p(i).lt.min_y) then
          min_y = y_p(i)
        endif		  
        if(z_p(i).lt.min_z) then
          min_z = z_p(i)
        endif
        if(x_p(i).gt.max_x) then
          max_x = x_p(i)
        endif
        if(y_p(i).gt.max_y) then
          max_y = y_p(i)
        endif
        if(z_p(i).gt.max_z) then
          max_z = z_p(i)
        endif
      enddo
c *** get the maximum and minimum for all the particles ***
c *** maximum ***
c      tdvs(0)=max_x
c      tdvs(1)=max_y
c      tdvs(2)=max_z
c      call MPI_ALLREDUCE(tdvs,tdvr,3,MPI_DOUBLE_PRECISION
c     &   ,MPI_MAX,MPI_COMM_WORLD,ierr)
c      max_x=tdvr(0)
c      max_y=tdvr(1)
c      max_z=tdvr(2)
c *** minimum ***
c      tdvs(0)=min_x
c      tdvs(1)=min_y
c      tdvs(2)=min_z
c      call MPI_ALLREDUCE(tdvs,tdvr,3,MPI_DOUBLE_PRECISION
c     &   ,MPI_MIN,MPI_COMM_WORLD,ierr)
c      min_x=tdvr(0)
c      min_y=tdvr(1)
c      min_z=tdvr(2)
c *** check x range ***
      tl=max_x - min_x
      if(tl.lt.0.0d0) then
        tl = -tl
      endif
      maxl=tl
c *** check y range ***
      tl=max_y - min_y
      if(tl.lt.0.0d0) then
        tl = -tl
      endif
      if(tl.gt.maxl) then
        maxl = tl
      endif		
c *** check z range ***
      tl=max_z-min_z
      if(tl.lt.0.0d0) then
        tl = -tl
      endif
      if(tl.gt.maxl) then
        maxl = tl
      endif
c *** Set root node ***
      l_tr(0)=MGROOT*maxl
      cx_tr(0)=(max_x+min_x)*0.5d0
      cy_tr(0)=(max_y+min_y)*0.5d0
      cz_tr(0)=(max_z+min_z)*0.5d0
      np_tr(0)=np
      
c      write(6,*) ' cx=',cx_tr(0),cy_tr(0),cz_tr(0)
c      write(6,*) ' l =',l_tr(0)
      return
      end

c *** Definition of compute_mass() ***
c * This function compute the center of mass and total mass *
c * tree, center of mass, total mass, node No. *	 	  
      subroutine compute_mass(level)
      include 'define.f'
      integer i,level
      integer j,d,npare,num,nd,pnd
      double precision xsd,ysd,zsd

c * compute total mass in node *
      do i=0,num_tr-1
        if(np_tr(i).eq.1) then
          mass_tr(i)=m_p(pn_tr(i))
          cmx_tr(i)=m_p(pn_tr(i))*x_p(pn_tr(i))
          cmy_tr(i)=m_p(pn_tr(i))*y_p(pn_tr(i))
          cmz_tr(i)=m_p(pn_tr(i))*z_p(pn_tr(i))
        else 
          mass_tr(i)=0.0d0
          cmx_tr(i)=0.0d0
          cmy_tr(i)=0.0d0
          cmz_tr(i)=0.0d0
        endif
      enddo

c *** initialization for set hm_tr ***
      do i=0,num_tr-1
        if(np_tr(i).eq.1) then
          hm_tr(i)=h_p(pn_tr(i))
        else
          hm_tr(i)=0.0d0
        endif
      enddo

      num=num_tr-1	
      do j=level,1,-1
        npare = 0
        do i=num,list(j),-1
          if(pare_tr(i).ne.pare_tr(i-1)) then
            pare_sc(npare) = pare_tr(i)
            nd_sc(npare) = i		  
            npare=npare+1
          endif
        enddo
        do d=0,7
!ocl novrec (mass_tr,cmx_tr,cmy_tr,cmz_tr)
!cdir nodep (mass_tr,cmx_tr,cmy_tr,cmz_tr)
          do i=0,npare-1
            nd=nd_sc(i)+d
            pnd=pare_sc(i)		  
            if(pnd.eq.pare_tr(nd).and.nd.le.num) then
              mass_tr(pnd)=mass_tr(pnd)+mass_tr(nd)
              cmx_tr(pnd)=cmx_tr(pnd)+cmx_tr(nd)
              cmy_tr(pnd)=cmy_tr(pnd)+cmy_tr(nd)
              cmz_tr(pnd)=cmz_tr(pnd)+cmz_tr(nd)
              if(hm_tr(nd).gt.hm_tr(pnd)) then
                hm_tr(pnd)=hm_tr(nd)
              endif
            endif
          enddo
        enddo		  
        num=list(j)-1
      enddo		
c * compute center of mass *
      do i=0,num_tr-1
        cmx_tr(i)=cmx_tr(i)/mass_tr(i)
        cmy_tr(i)=cmy_tr(i)/mass_tr(i)
        cmz_tr(i)=cmz_tr(i)/mass_tr(i)
      enddo
c *** compute delta ***
!cdir nodep
      do i=0,num_tr-1
        delta_tr(i)=dsqrt((cmx_tr(i)-cx_tr(i))*(cmx_tr(i)-cx_tr(i))
     &    +(cmy_tr(i)-cy_tr(i))*(cmy_tr(i)-cy_tr(i))
     &    +(cmz_tr(i)-cz_tr(i))*(cmz_tr(i)-cz_tr(i)))
      enddo		
c * compute Multipole Momentum *
c * initialization *
!cdir nodep
      do i=0,num_tr-1
        mx_tr(i)=0.0d0
        my_tr(i)=0.0d0
        mz_tr(i)=0.0d0
        mxx_tr(i)=0.0d0
        myy_tr(i)=0.0d0
        mzz_tr(i)=0.0d0
        mxy_tr(i)=0.0d0
        myz_tr(i)=0.0d0
        mzx_tr(i)=0.0d0 
      enddo
      num = num_tr-1	
      do j=level,1,-1
        npare = 0
        do i=num,list(j),-1
          if(pare_tr(i).ne.pare_tr(i-1)) then
            pare_sc(npare) = pare_tr(i)
            nd_sc(npare) = i
            npare=npare+1
          endif
        enddo
        do d=0,7
!ocl novrec(mxx_tr,myy_tr,mzz_tr,mxy_tr,myz_tr,mzx_tr)
!cdir nodep (mxx_tr,myy_tr,mzz_tr,mxy_tr,myz_tr,mzx_tr)
          do i=0,npare-1
            pnd=pare_sc(i)
            nd=nd_sc(i)+d
            if(pnd.eq.pare_tr(nd).and.nd.le.num) then
              xsd=cmx_tr(nd)-cmx_tr(pnd)
              ysd=cmy_tr(nd)-cmy_tr(pnd)
              zsd=cmz_tr(nd)-cmz_tr(pnd)
              mx_tr(pnd)=mx_tr(pnd)
     &          +(mx_tr(nd)-xsd*mass_tr(nd))
              my_tr(pnd)=my_tr(pnd)
     &          +(my_tr(nd)-ysd*mass_tr(nd))
              mz_tr(pnd)=mz_tr(pnd)
     &          +(mz_tr(nd)-zsd*mass_tr(nd))
              mxx_tr(pnd)=mxx_tr(pnd)
     &          +(mxx_tr(nd)-2.0d0*xsd*mx_tr(nd)
     &          +xsd*xsd*mass_tr(nd))
              myy_tr(pnd)=myy_tr(pnd)
     &          +(myy_tr(nd)-2.0d0*ysd*my_tr(nd)
     &           +ysd*ysd*mass_tr(nd))
              mzz_tr(pnd)=mzz_tr(pnd)
     &          +(mzz_tr(nd)-2.0d0*zsd*mz_tr(nd)
     &           +zsd*zsd*mass_tr(nd))
              mxy_tr(pnd)=mxy_tr(pnd)
     &          +(mxy_tr(nd)-xsd*my_tr(nd)-ysd*mx_tr(nd)
     &           +xsd*ysd*mass_tr(nd))
              myz_tr(pnd)=myz_tr(pnd)
     &          +(myz_tr(nd)-ysd*mz_tr(nd)-zsd*my_tr(nd)
     &           +ysd*zsd*mass_tr(nd))
              mzx_tr(pnd)=mzx_tr(pnd)
     &          +(mzx_tr(nd)-zsd*mx_tr(nd)-xsd*mz_tr(nd)
     &           +zsd*xsd*mass_tr(nd))
            endif
          enddo
        enddo
        num = list(j)-1
      enddo
c      write(6,*) ' num_tr = ',num_tr
c      write(6,*) ' cmx',cmx_tr(120)
c      write(6,*) ' mx,my,mz',mx_tr(120),my_tr(120),mz_tr(120)
c      write(6,*) ' mxx,myy,mzz',mxx_tr(120),myy_tr(120),
c     &  mzz_tr(120)
c      write(6,*) ' mxy,myz,mzx',mxy_tr(120),myz_tr(120),
c     &  mzx_tr(120)
      return
      end


