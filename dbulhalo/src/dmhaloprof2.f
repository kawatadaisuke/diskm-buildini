c /***************************************************
c   dmhaloprof2.f  for dbulhalo
c  1 Mar. 2011  written by D.Kawata
c ****************************************************/
c70*******************************************************************


c /*****   generating constant density DM halo or mext.dat *****/
      subroutine dmhaloprof2(mvir,rvir,fb)
      include 'define.f'

      double precision mvir,rvir
      double precision fb
c for work
      integer i
      character filen*60
c *** for radial profile ***
      double precision ri,ro,xcut,lri,lro,dlr,dr
      double precision mr,pmr,facnfw
      double precision rhodm
      double precision mhalo

      mhalo=mvir*(1.0d0-fb)
      rhodm=mhalo/(4.0d0*M_PI*(rvir**3)/3.0d0)

c /*** integral of NFW profile ****/
      ri = frini0*rvir
      ro = rvir
      lri = dlog10(ri)
      lro = dlog10(ro)
      dlr = (lro-lri)/dble(ndr-1)
c *** set common values
      lrodmh=lro
      lridmh=lri
      dlrdmh=dlr

      if(myrank.eq.0) then
        open(62,file='mext.dat',status='unknown'
     &   ,form='unformatted')
        open(63,file='mexta.dat',status='unknown')
        write(62) ndr+1
        write(62) dlr,lri,lro
        do i=0,ndr
          ro=10.0d0**(lri+dlr*dble(i))
          rdmh(i)=ro
          mrdmh(i)=mhalo*((ro/rvir)**3)
          rhodmh(i)=rhodm
          write(63,'(3(1pE13.5))') rdmh(i)*LUKPC,mrdmh(i)
     &     ,rhodmh(i)
          write(62) rdmh(i),mrdmh(i)
        enddo
        close(62)
        close(63)
      endif

      return
      end
