c *****************************************************
c    diskmetg.f
c 26 Jan., 2011   written by D. Kawata
c ***************************************************** 
c70*******************************************************************

      subroutine diskmetg(ngr,ngz,dfeh,feh0,dalfe,alfe0
     &  ,nhhalo,fehhalo)
      include './define.f'
      integer ngr,ngz
      double precision dfeh,feh0,dalfe,alfe0
      integer i,j,pn
      double precision rp,fehp,alfep,lmetp,metp
      double precision xsolh,xsolz,xh0,fehhalo,nhp,nhhalo

      xh0=0.76d0
      xsolh=0.706d0
      xsolz=0.019d0

      if(myrank.eq.0) then
        write(6,*) 'nhhalo,fehhalo=',nhhalo,fehhalo
      endif

      do i=0,ngr
        do j=0,ngz
          pn=id_rzg(i,j)
c check density
          fh_rzg(pn)=(xh0+(xsolh-xh0)*(10.0d0**fehhalo))        
          nhp=rho_rzg(pn)*fh_rzg(pn)*(DU/MP)
          if(nhp.lt.nhhalo) then
            fehp=fehhalo
          else
c *** radius 
            rp=x_rzg(pn)
            fehp=dfeh*(rp*LUKPC)+feh0
          endif
          alfep=dalfe*fehp+alfe0
c *** set [Z/H] from [Fe/H] and [al/Fe]
          lmetp=fehp
          if(fehp.lt.fehp+alfep) then
c if alpha overabundant 
            lmetp=fehp+alfep
          endif
          if(fehp.lt.fehp-alfep) then
c for the case of C or N enhanced 
            lmetp=fehp-alfep
          endif
          metp=10.0d0**lmetp
          fh_rzg(pn)=(xh0+(xsolh-xh0)*(10.0d0**fehp))
          lmet_rzg(pn)=dlog10((1.0d0/xsolh)
     &        *(xh0+(xsolh-xh0)*metp)*metp)
        enddo
      enddo

      return
      end
