c ****************************************************
c  copy from  set_value.F for gcd+ pver.32.17
c  23 jan. 2011    written by D.KAWATA
c ****************************************************

      subroutine set_den(npg)
      include 'define.f'
      integer npg,np
      integer i,j,nit,mnit,mnitb,npj,srank,rrank,ip,nval,isend
     &  ,irecv,nival,is
c * Number of Notfinished particle, temp *      
      integer nlist,tnlist
c * Particle Number in this node *      
      integer pni,pn,nd
      double precision r,s,rc,rh
c * Smoothing length *      
      double precision hsi,ihsi,wij
      double precision xij,yij,zij,tij
      double precision crit,ipi
c *** for work from common ***
c tx,ty,tz : temporary rotx,roty,rotz
      character fileo*60

c **** Basic Value ***
      ipi = 1.0d0/M_PI
      crit=dsqrt(3.0d0)*0.5d0

      nval = 4
c *** for communication ***
      np=npg
      npj=np
c
c
c reset flagc_p for gas
      do i=0,npg-1
        pn_nfp(i)=i
      enddo
      do i=0,np-1
        pni=pn_nfp(i)
        tdvr(i)=x_rzg(pni)
        tdvr(i+np)=y_rzg(pni)
        tdvr(i+np*2)=z_rzg(pni)
      enddo
c rho_p,omgh,div,rotx,roty,rotz,zetah
      do i=np*3,np*nval-1
        tdvr(i)=0.0d0
      enddo 

      srank = myrank+1
      if(srank.ge.nprocs) then
        srank = srank-nprocs
      endif
      rrank = myrank-1            
      if(rrank.lt.0) then
        rrank = rrank+nprocs
      endif  
c *** tree walk for all the proc ***
      npj = np
   70 do ip=0,nprocs-1
        if(np_tr(0).eq.0) then
          goto 99
        endif
        do i=0,npj-1
          list(i)=i
          node(i)=0
        enddo 
        nlist=npj
c **********    start tree walk   *********
   77   if(nlist.eq.0) then
          goto 99
        endif

        do i=0,nlist-1
          pni = list(i)
          nd = node(pni)         
          xij=tdvr(pni)-cx_tr(nd)
          yij=tdvr(pni+npj)-cy_tr(nd)
          zij=tdvr(pni+npj*2)-cz_tr(nd)
          rc=xij*xij+yij*yij+zij*zij
c *** only take into account gather condition ***
          rh=hm_tr(nd)+crit*l_tr(nd)
          rh=rh*rh
          if(rc.lt.rh.and.np_tr(nd).eq.1) then
            pn = pn_tr(nd)
c *** need to recalculate rij, because cx_tr is the center of the node. ***
            xij=tdvr(pni)-x_p(pn)
            yij=tdvr(pni+npj)-y_p(pn)
            zij=tdvr(pni+npj*2)-z_p(pn)
            r=xij*xij+yij*yij+zij*zij
            if(r.ne.0.0d0) then
              r=dsqrt(r)
            endif            
c *** calculate density for both gas and star,
c     in case if some neighbour are changed to stars ****
c *** use h_p(pn)
            hsi=h_p(pn)
            ihsi = 1.0d0/hsi
            s = r*ihsi
            if(s.lt.1.0d0) then
c * satisfy the criterion for neigbour particles *
c * calculate density if flag!=1 *
              is=int(s/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
              wij=w_tb(is)+(w_tb(is+1)-w_tb(is))*(s-s_tb(is))/ds_tb
              wij=wij*(ihsi**3)
c *** update density ***
              tdvr(pni+npj*3)=tdvr(pni+npj*3)+m_p(pn)*wij
            endif
          endif
c * update node *		
          if(np_tr(nd).eq.1.or.rc.ge.rh) then
            node(pni)=next_tr(nd)
          else
            node(pni)=daughter_tr(nd)
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
c        write(6,*) ' tnlist,nlist=',nlist,tnlist
        goto 77
c *** end itteration within the proc ***
c   99   if(nprocs.gt.1) then
cc *** communicating with the other proc ***
cc *** send and recieve the data ***
c          isend=npj*nval
c          do i=0,npj*nval-1
c            tdvs(i)=tdvr(i)
c          enddo
c          call MPI_ISEND(tdvs,isend,MPI_DOUBLE_PRECISION,srank,3,
c     &       MPI_COMM_WORLD,ireqs(2),ierr)
c          call MPI_IRECV(tdvr,MNVALP-1,MPI_DOUBLE_PRECISION,rrank,3,
c     &       MPI_COMM_WORLD,ireqr(2),ierr)
c          call MPI_WAIT(ireqs(2),istatus,ierr)
c          call MPI_WAIT(ireqr(2),istatusr,ierr)
c          call MPI_GET_COUNT(istatusr,MPI_DOUBLE_PRECISION,npj,ierr)
c          npj=npj/nval
cc *** tivr=flagc_p is also sent to track gather contribution ***
c          isend=isend/nval
c          isend=isend*nival
c          irecv=npj*nival
c          do i=0,isend-1
c            tivs(i)=tivr(i)
c          enddo
c          call MPI_ISEND(tivs,isend,MPI_INTEGER,srank,3,
c     &     MPI_COMM_WORLD,ireqs(1),ierr)
c          call MPI_IRECV(tivr,irecv,MPI_INTEGER,rrank,3,
c     &      MPI_COMM_WORLD,ireqr(1),ierr)
c          call MPI_WAIT(ireqs(1),istatus,ierr)
c            call MPI_WAIT(ireqr(1),istatus,ierr)
c        endif
      enddo

c *** update variables ***
 99   do i=0,npj-1
        pni = pn_nfp(i)
c        rho_rzg(pni)=tdvr(i+npj*3)
      enddo

      return
      end
