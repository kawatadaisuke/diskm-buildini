c *****************************************************
c    extf.f
c 19 Jan., 2011   written by D. Kawata
c ***************************************************** 
c70*******************************************************************

      subroutine extf(np,xp,yp,zp,fxe,fye,fze)
      include './define.f'
      integer np
      double precision xp(0:MN-1),yp(0:MN-1),zp(0:MN-1)
     &  ,fxe(0:MN-1),fye(0:MN-1),fze(0:MN-1)
      integer i,j,ir
      double precision pot,rp

      do i=0,np-1
        rp=dsqrt(xp(i)**2+yp(i)**2+zp(i)**2)
        if(rp.gt.10.0d0**lrodmh) then
          pot=mrdmh(ndr)/(rp**3)
        else if(rp.gt.0.0d0) then
          ir=int((dlog10(rp)-lridmh)/dlrdmh)
          if(ir.lt.0) then
            ir=0
          else if(ir.ge.ndr) then
            ir=ndr-1
          endif
          pot=(mrdmh(ir)+(rp-rdmh(ir))
     &     *(mrdmh(ir+1)-mrdmh(ir))/(rdmh(ir+1)-rdmh(ir)))
     &     /(rp**3)
        else
          pot=0.0d0
        endif
        fxe(i)=-xp(i)*pot
        fye(i)=-yp(i)*pot
        fze(i)=-zp(i)*pot
      enddo

      return
      end
