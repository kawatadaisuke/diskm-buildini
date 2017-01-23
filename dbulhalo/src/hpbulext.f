c /***************************************************
c   hpbulext.f  for dbulhalo
c  9 Feb. 2014  written by D.Kawata
c ****************************************************/
c70*******************************************************************


c /*****   generating NFW DM halo or mext.dat *****/
      subroutine hpbulext(mbul,ah,xcut,rt,flag)
      include 'define.f'

      double precision mbul,ah
c for work
      integer i,flag
      character filen*60
c *** for radial profile ***
      double precision ri,ro,xcut,lri,lro,dlr,dr,rt
      double precision mr,pmr
      double precision rhoi,rhoo
c *** external function 
      double precision rhohp,mrhp
      external rhohp,mrhp

      ri = frini0*rt
      ro = xcut*rt
      lri = dlog10(ri)
      lro = dlog10(ro)
      dlr = (lro-lri)/dble(ndr)
      ro = 0.0d0
      mr = 0.0d0
      do i=0,ndr
        ri=ro
        ro=10.0d0**(lri+dlr*dble(i))
        dr=ro-ri
        rhoo=rhohp(ro,mbul,ah)
        mrdmh(i)=mrdmh(i)+mrhp(ro,mbul,ah)
        rdmh(i) = ro
        rhodmh(i)=rhodmh(i)+rhoo
      enddo

c *** set common values
      lrodmh=lro
      lridmh=lri
      dlrdmh=dlr

      if(myrank.eq.0.and.flag.eq.-1) then
        open(62,file='mext.dat',status='unknown'
     &   ,form='unformatted')
        open(63,file='mexta.dat',status='unknown')
        write(62) ndr+1
        write(62) dlr,lri,lro
        do i=0,ndr
          write(63,'(3(1pE13.5))') rdmh(i)*LUKPC,mrdmh(i)
     &     ,rhodmh(i)
          write(62) rdmh(i),mrdmh(i)
        enddo
        close(62)
        close(63)
      endif

      return
      end
