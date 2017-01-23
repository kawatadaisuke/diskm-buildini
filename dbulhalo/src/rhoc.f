c /***************************************************
c   tempc.f  for dbulhalo ver.6
c  20 Jun. 2009  written by D.Kawata
c ****************************************************/
c70*******************************************************************

      subroutine tempc(tvir,rvir,zred,met,nrhg,rhg,rhohg,tecdhg)
      include 'define.f'
      double precision tvir,rvir
      double precision zred,met
      integer nrhg
      double precision rhg(MNR),rhohg(MNR),tecdhg(MNR)
      double precision tdyn,tc,dte,temp,ptemp
      double precision fx,hx,xx,dtc,pdtc
      integer i,j,it,mnit
      integer listb(0:MNB-1)
      double precision epsc
      parameter (epsc=1.0e-3)
      double precision tcool
      external tcool

      mnit=1000

c *** set tc=hage pressure line ***
      open(60,file='tempc.dat',status='unknown') 
      dte=0.05d0
      write(6,*) ' tvir=',tvir,' in tempc'
      do i=1,nrhg
        temp=tvir
        tdyn=dsqrt(3.0d0*M_PI/(16.0d0*G*rhohg(i)))
        do j=1,mnit
          ptemp=temp
c *** get temperature more than tvir ***
          temp=10.0d0**(dlog10(tvir)-dte*dble(j-1))
          tc=tcool(temp,zred,met,rhohg(i))

          pdtc=dtc
          dtc=tdyn-tc
          if(tc.gt.tdyn) then
            goto 96
          endif
        enddo
        write(6,*) ' Error cannot find T at r=',rhg(i)
        write(6,*) ' T,rho=',temp,rho_p(0)
        stop
 96     if(j.gt.1) then
          hx=ptemp
          xx=temp
          fx=dtc
          it=0
 97       it=it+1
          if(it.gt.mnit) then
            write(6,*) ' Error in tempc(): at r=',rhg(i)
            write(6,*) ' fx,hx,temp,dtc=',fx,hx,temp,dtc,dabs(dtc/tdyn)
            stop
          endif
          temp=0.5d0*(hx+xx)
          tc=tcool(temp,zred,met,rhohg(i))
          dtc=tdyn-tc
          if(dabs(dtc/tc).gt.epsc) then
            if(dtc*fx.gt.0.0d0) then
              xx=temp            
              fx=dtc
            else
              hx=temp
            endif
            goto 97
          endif
        endif
        write(60,*) rhg(i),rhohg(i),temp,rhohg(i)*temp*TPRHO,tc
        tecdhg(i)=temp
      enddo
      close(60)

      return
      end

      
